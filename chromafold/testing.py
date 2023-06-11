import os
import numpy as np

import torch
import torch.nn as nn

from tqdm import tqdm
from datetime import datetime

from utils import *

"""
Dataset for testing with both aggregated accessibility and coaccessibility
"""


class test_Dataset(torch.utils.data.Dataset):
    def __init__(
        self,
        input_size,
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

        input_window = int(self.input_size / pbulk_res)
        atac_start_ind = int(start / pbulk_res)
        sc_input_window = int(self.input_size / scatac_res)
        sc_atac_start_ind = int(start / scatac_res)

        X1 = self.pbulk_dict[chrom]
        X1 = X1[
            max(0, atac_start_ind + shift_input_by) : atac_start_ind
            + shift_input_by
            + input_window
        ]
        X1 = self.transform(X1.toarray().astype(np.float32).T)

        Z1 = self.scatac_dict[chrom]
        Z1 = Z1[
            max(0, sc_atac_start_ind + shift_input_by) : sc_atac_start_ind
            + shift_input_by
            + sc_input_window,
            :,
        ]

        Z1 = Z1.toarray().astype(np.float32).T
        Z1 = torch.tensor(Z1)

        X2 = self.ctcf_dict[chrom]
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

        X1_rev = torch.empty_like(X1).copy_(X1)
        X1_rev = torch.flip(X1_rev, [1])
        X2_rev = torch.empty_like(X2).copy_(X2)
        X2_rev = torch.flip(X2_rev, [1])
        Z1_rev = torch.empty_like(Z1).copy_(Z1)
        Z1_rev = torch.flip(Z1_rev, [1])

        return Z1, Z1_rev, X1, X1_rev, X2, X2_rev


def test(test_loader, model, device):
    '''
    Function for the validation step of the training loop
    '''
    y_hat_list = []
    y_true_list = []
    model.eval()
    
    with tqdm(test_loader, unit="batch") as tepoch:
        for Z2, Z2_rev, X1, X1_rev, X2, X2_rev in tepoch:
            # left interactions
            X1 = X1.to(device)
            X2 = X2.to(device)
            X = torch.cat((X1, X2), axis = 1)
            Z2 = Z2.to(device)
            y_hat = model(Z2, X) 
            
            # right interactions
            X1_rev = X1_rev.to(device)
            X2_rev = X2_rev.to(device)
            X_rev = torch.cat((X1_rev, X2_rev), axis = 1)
            Z2_rev = Z2_rev.to(device)
            y_hat_rev = model(Z2_rev, X_rev) 

            y_hat = torch.cat((y_hat, torch.flip(y_hat_rev, [1])), 1)
            y_hat_list.append(y_hat.detach().cpu())
        
    return y_hat_list

def test_model(model, test_loader, device):
    with torch.no_grad():
        y_hat_list = test(test_loader, model, device)
    return y_hat_list
