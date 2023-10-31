#!/usr/bin/env python

'''Script for processing scATAC fragment files.

Usage example: 
screen
# bsub -n 2 -W 4:00 -R 'span[hosts=1] rusage[mem=512]' -Is /bin/bash
bsub -q gpuqueue -gpu - -W 4:00 -n 2 -R 'span[hosts=1] rusage[mem=128]' -Is /bin/bash
source /home/yangr2/dnabert_environment/bin/activate
python /lila/data/leslie/yangr2/chromafold/scripts/scATAC_preparation.py \
--cell_type_prefix daniel_acute_d8_Effector \
--fragment_file /data/leslie/tyler/data/daniel-2022/acute-d8-tetpos/outs/fragments.tsv.gz \
--barcode_file /data/leslie/yangr2/chromafold/data/test_input/archr_data/archr_filtered_barcode.csv \
--lsi_file /data/leslie/yangr2/chromafold/data/test_input/archr_data/archr_filtered_lsi.csv \
--genome_assembly "mm10" \
--save_path /data/leslie/yangr2/chromafold/data/test_input/atac_data 

'''

import random
import sys
import numpy as np
import pysam
import pandas as pd
import os, argparse
import itertools
import scipy
from scipy.sparse import csr_matrix
from scipy import sparse
import pickle
from scipy.sparse import coo_matrix, csr_matrix, find
from sklearn.neighbors import NearestNeighbors

parser = argparse.ArgumentParser(description="Set-up data preparations",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--cell_type_prefix", help="Cell type and its prefix for the file names.",
                    type=str, default="")
parser.add_argument("--fragment_file", help="Path and name of the input fragment file.",
                    type=str, default="")
parser.add_argument("--barcode_file", help="Path and name of the valid barcode file.",
                    type=str, default="")
parser.add_argument("--genome_assembly", help="Genome assembly of the scATAC files.",
                    type=str, default="")
parser.add_argument("--lsi_file", help="Path and name of the LSI file from ArchR running.",
                    type=str, default="")
parser.add_argument("--save_path", help="Path to save all the intermediate and final files.",
                    type=str, default="")

args = parser.parse_args()
config = vars(args)
print(config)

CELL_TYPE = config['cell_type_prefix']
FRAG_FILE = config['fragment_file']
BARCODE_FILE = config['barcode_file']
GENOME_ASSEMBLY = config['genome_assembly']
LSI_PATH = config['lsi_file']
SAVE_PATH = config['save_path']

chrom_size = pd.read_csv(f'/home/yangr2/assembly/{GENOME_ASSEMBLY}.chrom.sizes', 
                         sep = '\t', header = None, index_col = 0)
valid_chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5','chr6','chr7',
               'chr8', 'chr9','chr10','chr11','chr12','chr13','chr14',
               'chr15','chr16','chr17','chr18','chr19','chr20',
               'chr21','chr22','chrX']
valid_chrom = [i for i in valid_chrom if i in chrom_size.index.to_list()]


#########################
#     Load functions    #
#########################

def atac_processing(frag_path, nrows_skip, valid_barcode_path,
                    save_path, cell_type_prefix):
    #1. read fragment file and filter for valid cell barcode
    frag = pd.read_csv(frag_path, sep = '\t', header=None, skiprows = nrows_skip)
    dist = frag[2]-frag[1]
    valid_barcode = pd.read_csv(valid_barcode_path, index_col = 0)
    # valid_barcode = valid_barcode.index
    valid_barcode = valid_barcode.iloc[:,0]
    valid_barcode = [str.split(i,"#")[1] for i in valid_barcode.to_list()]
    frag = frag[frag[3].isin(valid_barcode)]
    print("Finished loading fragment file.")
    
    ########### Create aggregated tile file for 500bp tiles ###########
    
    #2a. aggregate for each tile
    frag_stacked = pd.concat([frag[[0,1,3,4]],frag[[0,2,3,4]].rename(columns={2:1})], ignore_index=True)
    frag_stacked[4] = 1 #binarize count
    frag_stacked[1] = frag_stacked[1]//500 # aggregate counts per 500bp tile
    barcodes, barcode_id = np.unique(frag_stacked[3], return_inverse=True)
    frag_stacked[5] = barcode_id
    #2b. create tile files
    tile_dict = {}
    for chrom in list(set(frag_stacked[0])):
        if chrom in valid_chrom:
            frag_grouped = frag_stacked[frag_stacked[0]==chrom][[1,5,4]].groupby([5,1]).sum().reset_index()
            frag_grouped = csr_matrix((frag_grouped[4].astype(np.float64), 
                                       (frag_grouped[5].values, frag_grouped[1].values)), 
                                       shape=(barcodes.shape[0], chrom_size.loc[chrom].values[0]//500+1))
            tile_dict[chrom] = frag_grouped
            print(chrom)
    #2c. save to file
    pickle.dump(tile_dict, open(f"{save_path}/atac/{cell_type_prefix}_tile_500bp_dict.p", "wb"))
    np.save(f"{save_path}/atac/{cell_type_prefix}_tile_500bp_barcode.npy", barcodes.astype(str))
    del frag_stacked, tile_dict
    print("500bp tile data processed.")
    
    ########### Create pseudo-bulk data for evaluation ###########
    
    #3a. aggregate for 50bp tiles
    frag_stacked = pd.concat([frag[[0,1,4]],frag[[0,2,4]].rename(columns={2:1})], ignore_index=True)
    frag_stacked[4] = 1 #binarize counte
    frag_stacked[1] = frag_stacked[1]//50 # aggregate counts per 50bp tile
    tile_dict = {}
    for chrom in list(set(frag_stacked[0])):
        if chrom in valid_chrom:
            frag_grouped = frag_stacked[frag_stacked[0]==chrom][[1,4]].groupby([1]).sum().reset_index()
            frag_grouped = pd.DataFrame(frag_grouped[4].values, 
                                        index = frag_grouped[1].values).reindex(np.arange(0,
                                                                                          chrom_size.loc[chrom].values[0]//50+1)).fillna(0).values
            tile_dict[chrom] = frag_grouped
    #3b. save
    pickle.dump(tile_dict, open(f"{save_path}/atac/{cell_type_prefix}_tile_pbulk_50bp_dict.p", "wb"))
    del frag_stacked, tile_dict
    print("pseudo-bulk data processed.")
    
    ########### Create aggregated tile file for 50bp tiles ###########
    
    #4a. aggregate for each tile
    frag_stacked = pd.concat([frag[[0,1,3,4]],frag[[0,2,3,4]].rename(columns={2:1})], ignore_index=True)
    frag_stacked[4] = 1 #binarize counte
    frag_stacked[1] = frag_stacked[1]//50 # aggregate counts per 50bp tile
    barcodes, barcode_id = np.unique(frag_stacked[3], return_inverse=True)
    frag_stacked[5] = barcode_id
    #4b. create tile files
    tile_dict = {}
    for chrom in valid_chrom:
        frag_grouped = frag_stacked[frag_stacked[0]==chrom][[1,5,4]].groupby([5,1]).sum().reset_index()
        frag_grouped = csr_matrix((frag_grouped[4].astype(np.float64), 
                                   (frag_grouped[5], frag_grouped[1])), 
                                   shape=(barcodes.shape[0], chrom_size.loc[chrom].values[0]//50+1))
        tile_dict[chrom] = frag_grouped
    #4c. save
    pickle.dump(tile_dict, open(f"{save_path}/atac/{cell_type_prefix}_tile_50bp_dict.p", "wb"))  # save it into a file named save.p
    np.save(f"{save_path}/atac/{cell_type_prefix}_tile_50bp_barcode.npy",barcodes.astype(str))
    del frag_stacked, tile_dict
    print("50bp tile file processed.")
    
    return None

def get_max_overlap(y,x):
    return np.max(x.dot(y))

def generate_cicero_metacell(nbrs, max_overlap, sampled_id = [0]):
    order = np.arange(nbrs.shape[0])
    np.random.seed(10)
    np.random.shuffle(order)

    selected_idx = [False]*nbrs.shape[0]
    for i in sampled_id:
        selected_idx[i] = True
    
    for idx in order:
        selected_cells = nbrs[selected_idx]
        candidate = nbrs[idx]
        overlap = get_max_overlap(candidate,selected_cells)
        if overlap < max_overlap:
            selected_idx[idx] = True
    return selected_idx    

def metacell_computing(lsi_path, n_neighbors = 100, max_overlap = 33,
                      save_path = SAVE_PATH, cell_type_prefix = CELL_TYPE):
    lsi = pd.read_csv(lsi_path, index_col = 0)
    lsi.index = [x.split('#')[1] for x in lsi.index]
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, metric='euclidean').fit(lsi.values)
    nbrs = nbrs.kneighbors_graph(lsi.values).toarray()
    selected_idx = generate_cicero_metacell(nbrs, max_overlap = max_overlap)
    metacell_assignment = nbrs[selected_idx,:]
    pd.DataFrame(metacell_assignment).to_csv(f"{save_path}/atac/{cell_type_prefix}_metacell_mask.csv")
    return None

#########################
#     Script running    #
#########################

def main():
    args = parser.parse_args()
    config = vars(args)
    print(config)
    CELL_TYPE = config['cell_type_prefix']
    FRAG_FILE = config['fragment_file']
    BARCODE_FILE = config['barcode_file']
    GENOME_ASSEMBLY = config['genome_assembly']
    LSI_PATH = config['lsi_file']
    SAVE_PATH = config['save_path']

    chrom_size = pd.read_csv(f'/home/yangr2/assembly/{GENOME_ASSEMBLY}.chrom.sizes', 
                         sep = '\t', header = None, index_col = 0)
    
    _ = atac_processing(frag_path = FRAG_FILE, nrows_skip = 0, 
                        valid_barcode_path = BARCODE_FILE,
                        save_path = SAVE_PATH, cell_type_prefix = CELL_TYPE)
    _ = metacell_computing(lsi_path = LSI_PATH, save_path = SAVE_PATH, cell_type_prefix = CELL_TYPE)
    
if __name__ == '__main__':
    main()