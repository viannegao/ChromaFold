import argparse
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import torch.backends.cudnn as cudnn

import numpy as np
from numpy import savez_compressed

from tqdm import tqdm
import h5py

from utils import *
from dataloader import *
from model import *
from testing import *
from util_chrom_starts import *

# DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

parser = argparse.ArgumentParser(description="PyTorch chromaFold")

parser.add_argument(
    "--data-path", metavar="DIR", default="./datasets", help="path to dataset"
)
parser.add_argument(
    "--model-path", default="./checkpoints/chromafold_CTCFmotif.pth.tar", help="path to model checkpoint"
)
parser.add_argument(
    "--save-path", default="./predictions", help="path to save model predictions"
)

parser.add_argument("-ct", default="umb_endo")
parser.add_argument("-chrom", "--test-chrom", nargs="+", default=[])
parser.add_argument("--genome", default="hg38", help="genome of training cell type")
parser.add_argument(
    "-j",
    "--workers",
    default=8,
    type=int,
    help="number of data loading workers (default: 8)",
)
parser.add_argument(
    "-b", "--batch-size", default=64, type=int, help="batch size (default: 64)"
)
parser.add_argument("-offset", default=0, type = int, help="Chromosome start index offset (bp). For training, we use 4,000,000 to skip over the low-signal regions. For testing, use 0 to start at 2Mb and ensure all regions have full input context for prediction. we can use -2,000,000 to make prediction starting from the beginning; input will be zero-padded in this case.")

parser.add_argument("--use-ctcf-chip", action="store_true", help="Use CTCF ChIP seq data instead of motif score")
parser.add_argument("--disable-cuda", action="store_true", help="Disable CUDA")
parser.add_argument("--deterministic", action="store_true", help="cudnn deterministic")


def main():
    args = parser.parse_args()
    # check if gpu training is available
    if not args.disable_cuda and torch.cuda.is_available():
        args.device = torch.device("cuda")
        if args.deterministic:
            cudnn.deterministic = True
            logging.warning(
                "You have chosen to seed training. This will turn on the CUDNN deterministic setting, which can slow down "
                "your training considerably! You may see unexpected behavior when restarting from checkpoints."
            )
    #             cudnn.benchmark = True
    else:
        args.device = torch.device("cpu")

    data_path = args.data_path
    save_path = args.save_path
    model_path = args.model_path
    BATCH_SIZE = args.batch_size

    ct = args.ct
    genome = args.genome

    assert os.path.isdir(
        data_path
    ), "data_path does not exist. Please create directory first."
    assert os.path.exists(
        model_path
    ), "model_path does not exist. Please double check path."

    print("============================ Running ChromaFold ============================")
    print("\n-> Predicting Hi-C interactions in {}".format(ct))
    print("\n-> Using checkpoint {}".format(model_path))
    if args.use_ctcf_chip:
        print('\n-> Using CTCF ChIP seq data. Please make sure you are using the CTCF ChIP seq model.')

    test_chrom = args.test_chrom

    input_size = 4010000
    test_step = 1e04

    start_dict, end_dict = get_chrom_starts(genome)
    
    # LOAD DATA
    print("\n-> Params:")
    print(args)

    pbulk_dict = {}
    scatac_dict = {}
    n_cells = {}
    n_metacells = {}
    libsize_cell = {}

    pbulk_dict, scatac_dict = load_atac_eval(data_path, ct, pbulk_only=False)
    n_cells = get_num_cells(pbulk_dict)
    
    effective_genome_size = get_effective_genome_size(genome)
    normalized_pbulk_dict = normalize_bulk_dict(pbulk_dict, effective_genome_size)
    
    del pbulk_dict

    ctcf_dict = load_ctcf_motif(data_path, genome, args.use_ctcf_chip)

    test_dataset = test_Dataset(
        input_size,
        normalized_pbulk_dict,
        scatac_dict,
        ctcf_dict,
        *get_starts(
            chrom_list=test_chrom,
            start_dict=start_dict,
            end_dict=end_dict,
            chrom_start_offset=args.offset,
            chrom_end_offset=5000000,
            multi_ct=False,
            ct_list=[ct],
            step=test_step,
        ),
        transform=normalize
    )

    test_loader = DataLoader(
        dataset=test_dataset,
        batch_size=BATCH_SIZE,
        shuffle=False,
        num_workers=args.workers,
    )

    # Run inference
    print("Running Inference...")
    
    mod_branch_pbulk = nn.DataParallel(branch_pbulk(), device_ids=[0])
    mod_branch_cov = nn.DataParallel(branch_cov(), device_ids=[0])
    model = nn.DataParallel(trunk(mod_branch_pbulk, mod_branch_cov), device_ids=[0])
    model.to(args.device)
    model.load_state_dict(torch.load(model_path))
    model.eval()
    
    y_hat_list = test_model(model, test_loader, args.device)
    
    y_hat_list = [x.detach().cpu() for x in y_hat_list]
    y_hat_list = np.concatenate([x for x in 
                    y_hat_list], axis = 0)
    
    # Save predictions
    test_chrom_df = get_starts(
            chrom_list=test_chrom,
            start_dict=start_dict,
            end_dict=end_dict,
            chrom_start_offset=args.offset,
            chrom_end_offset=5000000,
            multi_ct=False,
            ct_list=[ct],
            step=test_step,
        )
    
    test_chrom_df = pd.DataFrame(test_chrom_df[2],
            test_chrom_df[2])

    for i in range(len(test_chrom)):
        start_end = np.where(test_chrom_df[0]=='chr{}'.format(test_chrom[i]))[0]
        print(start_end[0], start_end[-1]+1, test_chrom[i])
        
        savez_compressed('{}/prediction_{}_chr{}.npz'.format(save_path, ct, test_chrom[i]), 
                         y_hat_list[start_end[0] : start_end[-1] + 1, :])
        
    print("Inference done and predictions saved to {}.".format(save_path))


if __name__ == "__main__":
    main()
