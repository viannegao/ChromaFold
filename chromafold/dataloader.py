import numpy as np
import torch
import torch.nn as nn

from utils import *


"""
Dataset with only aggregated accessibility
"""


class Dataset_bulk_only(torch.utils.data.Dataset):
    def __init__(
        self,
        input_size,
        effective_genome_size,
        hic_dict,
        pbulk_dict,
        ctcf_dict,
        n_cells,
        libsize_cell,
        list_cts,
        list_starts,
        list_chroms,
        transform,
        is_training=False,
    ):
        self.input_size = input_size
        self.effective_genome_size = effective_genome_size
        self.hic_dict = hic_dict
        self.pbulk_dict = pbulk_dict
        self.ctcf_dict = ctcf_dict
        self.n_cells = n_cells
        self.libsize_cell = libsize_cell
        self.list_cts = list_cts
        self.list_starts = list_starts
        self.list_chroms = list_chroms
        self.is_training = is_training
        self.transform = transform

    def __len__(self):
        return len(self.list_starts)

    def __getitem__(self, index):
        ct = self.list_cts[index]
        chrom = self.list_chroms[index]
        start = self.list_starts[index]

        shift_input_by = 0
        train_step = 5e04
        pbulk_res = 50

        if self.is_training:
            start = (
                np.int64(np.random.randint(0, np.int64(train_step / 1e4)) * 1e4) + start
            )
            # inject noise by adding shifts to the input
            shift_input_by = np.random.randint(-1, 2)

        input_window = np.int64(self.input_size / pbulk_res)
        atac_start_ind = np.int64(start / pbulk_res)

        hic = self.hic_dict[ct][chrom]
        X1 = self.pbulk_dict[ct][chrom]
        X1 = X1[
            :,
            atac_start_ind
            + shift_input_by : atac_start_ind
            + shift_input_by
            + input_window,
        ]
        X2 = self.ctcf_dict[chrom]
        X2 = X2[atac_start_ind : atac_start_ind + input_window]
        X2 = torch.tensor(X2.toarray().astype(np.float32).T)

        if self.is_training:
            # Choose the number of cells to sample
            size = np.random.randint(500, 5000)
            mask = np.random.choice(
                np.arange(self.n_cells[ct]), size=size, replace=True
            )
            X1 = X1[mask]
            cur_libsize = self.libsize_cell[ct][mask].sum()
        else:
            cur_libsize = self.libsize_cell[ct].sum()

        X1 = np.array(X1.sum(0)).astype(np.float32)
        X1 = self.transform(
            X1,
            effective_genome_size=self.effective_genome_size,
            libsize=cur_libsize,
            training=True,
        )

        y = get_y_vstripe(start, hic, self.input_size)
        y = torch.tensor(y)

        X1_rev = torch.empty_like(X1).copy_(X1)
        X1_rev = torch.flip(X1_rev, [1])
        X2_rev = torch.empty_like(X2).copy_(X2)
        X2_rev = torch.flip(X2_rev, [1])
        y_rev = torch.flip(y, [0])

        y = y[0:200]
        y_rev = y_rev[0:200]

        return X1, X1_rev, X2, X2_rev, y, y_rev


"""
Dataset with both aggregated accessibility and coaccessibility
"""


class Dataset(torch.utils.data.Dataset):
    def __init__(
        self,
        input_size,
        effective_genome_size,
        hic_dict,
        pbulk_dict,
        scatac_dict,
        ctcf_dict,
        n_cells,
        n_metacells,
        libsize_cell,
        list_cts,
        list_starts,
        list_chroms,
        transform,
        is_training=False,
    ):
        self.input_size = input_size
        self.effective_genome_size = effective_genome_size
        self.hic_dict = hic_dict
        self.pbulk_dict = pbulk_dict
        self.scatac_dict = scatac_dict
        self.ctcf_dict = ctcf_dict
        self.n_cells = n_cells
        self.n_metacells = n_metacells
        self.libsize_cell = libsize_cell
        self.list_cts = list_cts
        self.list_starts = list_starts
        self.list_chroms = list_chroms
        self.is_training = is_training
        self.transform = transform

    def __len__(self):
        return len(self.list_starts)

    def __getitem__(self, index):
        ct = self.list_cts[index]
        chrom = self.list_chroms[index]
        start = self.list_starts[index]

        shift_input_by = 0
        train_step = 5e04
        pbulk_res = 50
        scatac_res = 500

        if self.is_training:
            start = (
                np.int32(np.random.randint(0, np.int32(train_step / 1e4)) * 1e4) + start
            )
            # inject noise by adding shifts to the input
            shift_input_by = np.random.randint(-1, 2)

        hic = self.hic_dict[ct][chrom]
        input_window = np.int32(self.input_size / pbulk_res)
        atac_start_ind = np.int32(start / pbulk_res)
        sc_input_window = np.int32(self.input_size / scatac_res)
        sc_atac_start_ind = np.int32(start / scatac_res)

        X1 = self.pbulk_dict[ct][chrom]
        X1 = X1[
            :,
            atac_start_ind
            + shift_input_by : atac_start_ind
            + shift_input_by
            + input_window,
        ]
        Z1 = self.scatac_dict[ct][chrom]
        Z1 = Z1[
            sc_atac_start_ind
            + shift_input_by : sc_atac_start_ind
            + shift_input_by
            + sc_input_window
        ]
        X2 = self.ctcf_dict[chrom]
        X2 = X2[atac_start_ind : atac_start_ind + input_window]
        X2 = torch.tensor(X2.toarray().astype(np.float32).T)

        if self.is_training:
            # Choose the number of cells to sample
            size = np.random.randint(500, 5000)
            mask = np.random.choice(
                np.arange(self.n_cells[ct]), size=size, replace=True
            )
            X1 = X1[mask]
            cur_libsize = self.libsize_cell[ct][mask].sum()
            # Choose the number of metacells to sample
            size2 = np.random.randint(400, 1000)
            mask2 = np.random.choice(
                np.arange(self.n_metacells[ct]), size=size2, replace=True
            )
            Z1 = Z1.T[mask2]

        else:
            Z1 = Z1.T
            cur_libsize = self.libsize_cell[ct].sum()

        X1 = np.array(X1.sum(0)).astype(np.float32)
        X1 = self.transform(
            X1,
            effective_genome_size=self.effective_genome_size,
            libsize=cur_libsize,
            training=True,
        )

        Z1 = Z1.toarray().astype(np.float32)
        Z1 = torch.tensor(Z1)

        Z1 = cpu_batch_corcoeff_vstripe(Z1)

        y = get_y_vstripe(start, hic, self.input_size)
        y = torch.tensor(y)

        X1_rev = torch.empty_like(X1).copy_(X1)
        X1_rev = torch.flip(X1_rev, [1])

        Z1_rev = torch.empty_like(Z1).copy_(Z1)
        Z1_rev = torch.flip(Z1_rev, [1])

        X2_rev = torch.empty_like(X2).copy_(X2)
        X2_rev = torch.flip(X2_rev, [1])

        y_rev = torch.flip(y, [0])

        y = y[0:200]
        y_rev = y_rev[0:200]

        return Z1, Z1_rev, X1, X1_rev, X2, X2_rev, y, y_rev


"""
Dataset for testing with both aggregated accessibility and coaccessibility
"""


class test_Dataset(torch.utils.data.Dataset):
    def __init__(
        self,
        input_size,
        hic_dict,
        pbulk_dict,
        scatac_dict,
        ctcf_dict,
        list_cts,
        list_starts,
        list_chroms,
        transform,
        is_training=False,
    ):
        self.input_size = input_size
        self.hic_dict = hic_dict
        self.pbulk_dict = pbulk_dict
        self.scatac_dict = scatac_dict
        self.ctcf_dict = ctcf_dict
        self.list_cts = list_cts
        self.list_starts = list_starts
        self.list_chroms = list_chroms
        self.is_training = is_training
        self.transform = transform

    def __len__(self):
        return len(self.list_starts)

    def __getitem__(self, index):
        chrom = self.list_chroms[index]
        start = self.list_starts[index]

        shift_input_by = 0
        pbulk_res = 50
        scatac_res = 500

        hic = hic_dict[chrom]
        input_window = int(input_size / pbulk_res)
        atac_start_ind = int(start / pbulk_res)
        sc_input_window = int(input_size / scatac_res)
        sc_atac_start_ind = int(start / scatac_res)

        X1 = lib_normalizer[chrom]
        X1 = X1[
            max(0, atac_start_ind + shift_input_by) : atac_start_ind
            + shift_input_by
            + input_window
        ]
        X1 = self.transform(X1.toarray().astype(np.float32).T)

        Z1 = tile_dict[chrom]
        Z1 = Z1[
            max(0, sc_atac_start_ind + shift_input_by) : sc_atac_start_ind
            + shift_input_by
            + sc_input_window,
            :,
        ]

        Z1 = Z1.toarray().astype(np.float32).T
        Z1 = torch.tensor(Z1)

        X2 = ctcf_dict[chrom]
        X2 = X2[max(0, atac_start_ind) : atac_start_ind + input_window]
        X2 = torch.tensor(X2.toarray().astype(np.float32).T)

        if start < 0:
            desired_shape = (1, 80200)
            # Determine the number of zeros to add on the left side
            num_zeros = desired_shape[1] - X1.shape[1]
            # Pad the tensor with zeros on the left side
            X1 = torch.cat([torch.zeros((X1.shape[0], num_zeros)), X1], dim=1)
            X2 = torch.cat([torch.zeros((X2.shape[0], num_zeros)), X2], dim=1)

            desired_shape = (Z1.shape[0], 8020)
            num_zeros = desired_shape[1] - Z1.shape[1]
            Z1 = torch.cat([torch.zeros((Z1.shape[0], num_zeros)), Z1], dim=1)

        Z1 = cpu_batch_corcoeff_vstripe(Z1)

        y = get_y_vstripe_eval(start, hic, input_size)
        y = torch.tensor(y)

        X1_rev = torch.empty_like(X1).copy_(X1)
        X1_rev = torch.flip(X1_rev, [1])
        X2_rev = torch.empty_like(X2).copy_(X2)
        X2_rev = torch.flip(X2_rev, [1])
        Z1_rev = torch.empty_like(Z1).copy_(Z1)
        Z1_rev = torch.flip(Z1_rev, [1])

        y_rev = torch.flip(y, [0])

        y = y[0:200]
        y_rev = y_rev[0:200]

        return Z1, Z1_rev, X1, X1_rev, X2, X2_rev, y, y_rev
