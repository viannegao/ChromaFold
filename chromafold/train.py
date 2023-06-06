import argparse
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import torch.backends.cudnn as cudnn

import numpy as np
from tqdm import tqdm
import h5py

from utils import *
from dataloader import *
from model import *
from training_utils import *

# DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'

parser = argparse.ArgumentParser(description="PyTorch chromaFold")

parser.add_argument(
    "--data-path", metavar="DIR", default="./datasets", help="path to dataset"
)
parser.add_argument(
    "--save-path", default="./checkpoints", help="path to save model chekpoints"
)
parser.add_argument(
    "--bulk-checkpoint",
    default="./checkpoints/deterministic/gm12878_hg38_umb_endo_imr90_CTCFmotifScore_seed1229.pth.tar",
    help="pbulk chekpoint",
)
parser.add_argument("-ct", "--ct-list", nargs="+", default=[])
parser.add_argument("--genome", default="hg38", help="genome of training cell type")
parser.add_argument(
    "--min-stripe-signal",
    default=-170,
    type=int,
    help="Stripe signal cutoff; avoid training on low-signal stripes",
)
parser.add_argument(
    "-j",
    "--workers",
    default=8,
    type=int,
    help="number of data loading workers (default: 8)",
)
parser.add_argument(
    "--epochs", default=30, type=int, help="number of total epochs to run"
)
parser.add_argument(
    "-b", "--batch-size", default=64, type=int, help="batch size (default: 64)"
)
parser.add_argument(
    "--patience", default=10, type=int, help="training patience in epochs (default: 10)"
)
parser.add_argument(
    "--lr",
    "--learning-rate",
    default=1e-6,
    type=float,
    metavar="LR",
    help="initial learning rate",
    dest="lr",
)
parser.add_argument(
    "--wd",
    "--weight-decay",
    default=1e-3,
    type=float,
    metavar="W",
    help="weight decay (default: 1e-3)",
    dest="weight_decay",
)
parser.add_argument(
    "--seed", default=1229, type=int, help="seed for initializing training. "
)
parser.add_argument("--disable-cuda", action="store_true", help="Disable CUDA")
parser.add_argument("--deterministic", action="store_true", help="cudnn deterministic")
parser.add_argument("--mod-name", default="", help="model name")


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

    torch.manual_seed(args.seed)

    data_path = args.data_path
    save_path = args.save_path
    BATCH_SIZE = args.batch_size
    N_EPOCHS = args.epochs
    initial_rate = args.lr
    wd = args.weight_decay

    ct_list = args.ct_list
    genome = args.genome

    # Verify validity of input and output paths.
    assert os.path.isdir(
        data_path
    ), "data_path does not exist. Please create directory first."
    assert os.path.isdir(
        save_path
    ), "save_path does not exist. Please create directory first."

    if args.mod_name:
        save_path = os.path.join(save_path, args.mod_name)
        if not os.path.isdir(save_path):
            os.mkdir(save_path)

    print("Saving to {}".format(save_path))

    # Training specifications.
    train_chrom = [1, 2, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19]
    val_chrom = [3, 15]
    test_chrom = [5, 18, 20, 21]

    # The following should not be changed unless the user modifies the dimensions in the model accordingly.
    input_size = 4010000
    train_step = 5e04
    val_step = 5e04

    #     start_dict = {}
    #     end_dict = {}
    #     for i in hic_dict[ct_list[1]].keys():
    #         idx = np.where(hic_dict[ct_list[0]][i].toarray().sum(1)>-600)[0]
    #         start_dict[i] = idx[0] * 10000
    #         end_dict[i] = idx[-1] * 10000
    #     for ct in ct_list[1:]:
    #         for i in hic_dict[ct_list[1]].keys():
    #             idx = np.where(hic_dict[ct][i].toarray().sum(1)>-600)[0]
    #             start_dict[i] = max(start_dict[i], idx[0] * 10000)
    #             end_dict[i] = min(end_dict[i], idx[-1] * 10000)

    # Start and end indices of each chromosome
    start_dict = {
        "chr13": 0,
        "chr21": 0,
        "chr3": 0,
        "chr20": 0,
        "chrX": 0,
        "chr8": 0,
        "chr10": 0,
        "chr15": 0,
        "chr5": 0,
        "chr2": 0,
        "chr14": 0,
        "chr6": 0,
        "chr7": 0,
        "chr17": 0,
        "chr18": 0,
        "chr16": 0,
        "chr1": 0,
        "chr11": 0,
        "chr4": 0,
        "chr9": 0,
        "chr12": 0,
        "chr19": 0,
    }

    end_dict = {
        "chr13": 114360000,
        "chr21": 46700000,
        "chr3": 198020000,
        "chr20": 63020000,
        "chrX": 155270000,
        "chr8": 145130000,
        "chr10": 133790000,
        "chr15": 101990000,
        "chr5": 180910000,
        "chr2": 242190000,
        "chr14": 107040000,
        "chr6": 170800000,
        "chr7": 159130000,
        "chr17": 81190000,
        "chr18": 78070000,
        "chr16": 90330000,
        "chr1": 248950000,
        "chr11": 135000000,
        "chr4": 190210000,
        "chr9": 138390000,
        "chr12": 133270000,
        "chr19": 58610000,
    }

    # Specify efective genome size for normalization
    if "hg" in genome:
        effective_genome_size = 2913022398
    elif "mm" in genome:
        effective_genome_size = 2652783500
    else:
        raise ValueError(
            "Please compute effective_genome_size manuallly and modify this function accordingly."
        )

    # LOAD DATA
    print("Loading training data for these cell types:")
    print(ct_list)
    print(args)

    hic_dict = {}
    hic_qval_dict = {}
    pbulk_dict = {}
    scatac_dict = {}
    n_cells = {}
    n_metacells = {}
    libsize_cell = {}

    for ct in ct_list:
        hic_dict[ct], hic_qval_dict[ct] = load_hic(data_path, ct)
        pbulk_dict[ct], scatac_dict[ct] = load_atac(data_path, ct, pbulk_only=False)
        n_cells[ct] = get_num_cells(pbulk_dict[ct])
        n_metacells[ct] = get_num_cells(scatac_dict[ct], dim=1)
        libsize_cell[ct] = get_libsize(pbulk_dict[ct])

    ctcf_dict = load_ctcf_motif(data_path, genome)

    # Create dataset objects for training and validation
    train_dataset = Dataset(
        input_size,
        effective_genome_size,
        hic_dict,
        pbulk_dict,
        scatac_dict,
        ctcf_dict,
        n_cells,
        n_metacells,
        libsize_cell,
        *get_starts(
            chrom_list=train_chrom,
            start_dict=start_dict,
            end_dict=end_dict,
            chrom_start_offset=4000000,
            chrom_end_offset=5000000,
            multi_ct=True,
            ct_list=ct_list,
            step=train_step,
        ),
        transform=normalize,
        is_training=True
    )

    val_dataset = Dataset(
        input_size,
        effective_genome_size,
        hic_dict,
        pbulk_dict,
        scatac_dict,
        ctcf_dict,
        n_cells,
        n_metacells,
        libsize_cell,
        *get_starts(
            chrom_list=val_chrom,
            start_dict=start_dict,
            end_dict=end_dict,
            chrom_start_offset=4000000,
            chrom_end_offset=5000000,
            multi_ct=True,
            ct_list=ct_list,
            step=val_step,
        ),
        transform=normalize
    )

    # Create dataloader objects for training and validation
    train_loader = DataLoader(
        dataset=train_dataset,
        batch_size=BATCH_SIZE,
        shuffle=True,
        num_workers=args.workers,
    )
    val_loader = DataLoader(
        dataset=val_dataset,
        batch_size=BATCH_SIZE,
        shuffle=False,
        num_workers=args.workers,
    )

    # Define Model
    mod_branch_pbulk = nn.DataParallel(branch_pbulk(), device_ids=[0])
    PATH_branch_pbulk = args.bulk_checkpoint
    mod_branch_pbulk.load_state_dict(torch.load(PATH_branch_pbulk), strict=True)

    mod_branch_cov = nn.DataParallel(branch_cov(), device_ids=[0])
    model = nn.DataParallel(trunk(mod_branch_pbulk, mod_branch_cov), device_ids=[0])

    optimizer = torch.optim.SGD(
        model.parameters(), lr=initial_rate, momentum=0.9, weight_decay=wd
    )
    scheduler = torch.optim.lr_scheduler.MultiStepLR(
        optimizer, milestones=[N_EPOCHS / 10 * 4, N_EPOCHS / 10 * 8], gamma=0.1
    )

    criterion = weighted_mse_loss

    # Train
    print("Training...")
    model, optimizer, _ = training_loop(
        train,
        validate,
        model,
        ct_list,
        args.seed,
        criterion,
        optimizer,
        scheduler,
        train_loader,
        val_loader,
        N_EPOCHS,
        args.patience,
        save_path,
        args.device,
        args.min_stripe_signal,
    )


if __name__ == "__main__":
    main()
