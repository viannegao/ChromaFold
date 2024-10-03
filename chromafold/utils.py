import pickle
import os

import numpy as np
import pandas as pd
import scipy
from scipy import sparse 
import torch


def get_effective_genome_size(genome):
    if "hg" in genome:
        effective_genome_size = 2913022398
    elif "mm" in genome:
        effective_genome_size = 2652783500
    else:
        raise ValueError(
            "Please compute effective_genome_size manuallly and modify this function accordingly."
        )
    return effective_genome_size


def get_metacell_profile(sc_atac_dict, nbrs):
    metacell_sc_atac_dict = {}
    metacell = nbrs
    for chrom in list(sc_atac_dict.keys()):
        metacell_sc_atac_dict[chrom] = (
            sparse.csr_matrix(metacell) * sc_atac_dict[chrom]
        )
        metacell_sc_atac_dict[chrom] = sparse.csr_matrix(
            metacell_sc_atac_dict[chrom].toarray().transpose()
        )

    return metacell_sc_atac_dict


def load_hic(data_path, ct):
    zscore_dict = pickle.load(
        open(os.path.join(data_path, "hic/{}_hic_zscore_dict.p".format(ct)), "rb")
    )
    # qvalue_dict = pickle.load(
    #     open(os.path.join(data_path, "hic/{}_hic_qvalue_dict.p".format(ct)), "rb")
    # )
    qvalue_dict = pickle.load(
        open(os.path.join(data_path, "hic/{}_hic_zscore_dict.p".format(ct)), "rb")
    )

    return zscore_dict, qvalue_dict


def load_atac(data_path, ct, pbulk_only=False):
    bulk_atac_dict = pickle.load(
        open(os.path.join(data_path, "atac/{}_tile_50bp_dict.p".format(ct)), "rb")
    )
    if pbulk_only:
        sc_atac_dict = None
    else:
        sc_atac_dict = pickle.load(
            open(os.path.join(data_path, "atac/{}_tile_500bp_dict.p".format(ct)), "rb")
        )
        metacell_mask = pd.read_csv(
            os.path.join(data_path, "atac/{}_metacell_mask.csv".format(ct)), index_col=0
        ).values
        sc_atac_dict = get_metacell_profile(sc_atac_dict, metacell_mask)

    return bulk_atac_dict, sc_atac_dict


def load_atac_eval(data_path, ct, pbulk_only=False):
    bulk_atac_dict = pickle.load(
        open(os.path.join(data_path, "atac/{}_tile_pbulk_50bp_dict.p".format(ct)), "rb")
    )
    
    for chrom in bulk_atac_dict.keys():
        bulk_atac_dict[chrom] = sparse.csr_matrix(bulk_atac_dict[chrom])
        
    if pbulk_only:
        sc_atac_dict = None
    else:
        sc_atac_dict = pickle.load(
            open(os.path.join(data_path, "atac/{}_tile_500bp_dict.p".format(ct)), "rb")
        )
        metacell_mask = pd.read_csv(
            os.path.join(data_path, "atac/{}_metacell_mask.csv".format(ct)), index_col=0
        ).values
        sc_atac_dict = get_metacell_profile(sc_atac_dict, metacell_mask)

    return bulk_atac_dict, sc_atac_dict


def load_ctcf_motif(data_path, genome, use_chip = False):
    if not use_chip:
        ctcf_dict = pickle.load(
            open(os.path.join(data_path, "dna/{}_ctcf_motif_score.p".format(genome)), "rb")
        )
    else:
        ctcf_dict = pickle.load(
            open(os.path.join(data_path, "dna/{}_CTCF_50bp_dict.p".format(genome)), "rb")
        )
    for chrom in ctcf_dict.keys():
        ctcf_dict[chrom] = ctcf_dict[chrom].T

    return ctcf_dict


def get_num_cells(atac_dict, dim=0):
    return atac_dict["chr1"].shape[dim]


def get_libsize(bulk_atac_dict):
    libsize_cell = np.zeros(get_num_cells(bulk_atac_dict))
    for chrom in bulk_atac_dict.keys():
        libsize_cell += np.array(bulk_atac_dict[chrom].sum(1).T)[0]

    return libsize_cell


def normalize(tensor, effective_genome_size=0, libsize=0, training=False):
    tensor = torch.tensor(tensor)
    if training:
        assert libsize != 0, "Libsize cannot be 0 during training."
        scale_factor = (libsize * 150) / effective_genome_size
        tensor.div_(torch.as_tensor(scale_factor))
    tensor = torch.log(tensor + 1)
    return tensor


def normalize_bulk_dict(pbulk_dict, effective_genome_size):
    normalized_pbulk_dict = {}
    libsize = 0
    for chrom in pbulk_dict.keys():
        libsize += pbulk_dict[chrom].toarray().sum()
    scale_factor = (libsize * 150) / (effective_genome_size)
    for chrom in pbulk_dict.keys():
        y = pbulk_dict[chrom].toarray()
        y_t = y / scale_factor
        normalized_pbulk_dict[chrom] = sparse.csr_matrix(y_t)
    return normalized_pbulk_dict


def get_y_vstripe(start, hic, input_size, resolution=1e04):
    hic_start_row = start
    hic_start_col = start + 200 * np.int64(resolution)
    ind_row = np.int64(hic_start_row // resolution)
    ind_col = np.int64(hic_start_col // resolution)
    vert = hic[ind_row : ind_row + 200, ind_col].toarray()

    hic_start_row = start + 200 * np.int64(resolution)
    hic_start_col = start + 201 * np.int64(resolution)
    ind_row = np.int64(hic_start_row // resolution)
    ind_col = np.int64(hic_start_col // resolution)

    hori = hic[ind_row, ind_col : ind_col + 200].toarray()

    y = np.append(vert, hori).astype(np.float32)
    y = y.clip(-16, 16)
    return y


def get_y_vstripe_eval(start, hic, input_size, resolution=1e04):
    hic_start_row = start
    hic_start_col = start + 200 * np.int64(resolution)
    ind_row = np.int64(hic_start_row // resolution)
    ind_col = np.int64(hic_start_col // resolution)
    vert = hic[max(0, ind_row) : ind_row + 200, ind_col].toarray()

    desired_length = 200
    num_zeros = desired_length - vert.shape[0]
    vert = np.pad(vert, ((num_zeros, 0), (0, 0)), "constant")

    hic_start_row = start + 200 * np.int64(resolution)
    hic_start_col = start + 201 * np.int64(resolution)
    ind_row = np.int64(hic_start_row // resolution)
    ind_col = np.int64(hic_start_col // resolution)

    hori = hic[ind_row, ind_col : ind_col + 200].toarray()
    num_zeros = desired_length - hori.shape[1]
    hori = np.pad(hori, ((0, 0), (0, num_zeros)), "constant")

    y = np.append(vert, hori).astype(np.float32)
    return y


def cpu_jaccard_vstripe(x):
    # calculate jaccard similarity of rows
    scatac_res = 500
    size = x.shape[1]
    eps = 1e-8
    i = np.int16(1000 / scatac_res)

    x = torch.where(x > 0.0, torch.tensor([1.0]), torch.tensor([0.0]))
    num = torch.mm(
        x[2000 * i : 2010 * i, :],
        x[np.r_[:, 0 : 2000 * i, 2010 * i : 4010 * i]].transpose(0, 1),
    )

    x = torch.where(x == 0.0, torch.tensor([1.0]), torch.tensor([0.0]))
    denom = torch.mm(
        x[2000 * i : 2010 * i, :],
        x[np.r_[:, 0 : 2000 * i, 2010 * i : 4010 * i]].transpose(0, 1),
    )
    denom = size - denom

    num = torch.div(num, torch.max(denom, eps * torch.ones_like(denom)))

    return num


def cpu_batch_corcoeff_vstripe(x):
    c = cpu_jaccard_vstripe(x.permute(1, 0))
    c[c != c] = 0
    return c


def get_starts(
    chrom_list,
    start_dict,
    end_dict,
    chrom_start_offset,
    chrom_end_offset,
    multi_ct=False,
    ct_list=[],
    step=5e4,
):
    """
    Function to get the start and chrom list for train/test/val
    """
    ctl = []
    startl = []
    chroml = []
    for chrom in chrom_list:
        chrom = "chr{}".format(chrom)
        cur_starts = list(
            np.arange(
                start_dict[chrom] + chrom_start_offset,
                end_dict[chrom] - chrom_end_offset,
                step,
            ).astype(int)
        )
        startl = startl + cur_starts
        chroml = chroml + list(np.repeat(chrom, len(cur_starts)))

    if multi_ct:
        for ct in ct_list:
            ctl += list(np.repeat(ct, len(chroml)))
        startl = startl * len(ct_list)
        chroml = chroml * len(ct_list)
        assert (len(startl) == len(chroml)) & (len(ctl) == len(chroml))

    return ctl, startl, chroml
