#!/usr/bin/env python

'''Script for processing scATAC fragment files.

Usage example: 
screen

python /chromafold/scripts/fragment_celltype_merge.py \
--cell_type selected_cell_type \
--fragment_list /data/fragments_1.tsv.gz /data/fragments_2.tsv.gz \
--data_prefix_list "data_prefix" \
--save_name /data/merged_fragments.tsv.gz \
--cell_type_file /data/cell_type_file.csv \
--genome_assembly mm10

'''

import pandas as pd
import numpy as np 
import argparse

parser = argparse.ArgumentParser(description="Set-up data preparations",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--cell_type", help="Which cell type to use in the fragment files.",
                    type=str, default="")
parser.add_argument("--fragment_list", help="Path and name of the input fragment files.",
                    nargs='+', type=str, default="")
parser.add_argument("--data_prefix_list", help="Prefix for the fragment files (length must match fragment_list).",
                    nargs='+', type=str, default="")
parser.add_argument("--save_name", help="Path and name of the filtered fragment file for specific cell type.",
                    type=str, default="")
parser.add_argument("--cell_type_file", help="General cell type list for the fragment files.",
                    type=str, default="")
parser.add_argument("--genome_assembly", help="Genome assembly of the scATAC files.",
                    type=str, default="")

args = parser.parse_args()
config = vars(args)
print(config)

CELL_TYPE = config['cell_type']
FRAG_LIST = config['fragment_list']
DATA_PREFIX_LIST = config['data_prefix_list']
SAVE_NAME = config['save_name']
CELL_TYPE_FILE = config['cell_type_file']
GENOME_ASSEMBLY = config['genome_assembly']

#########################
#     Load functions    #
#########################

def filter_celltype(FRAG_LIST, DATA_PREFIX_LIST, cell_type, CELL_TYPE_FILE, SAVE_NAME, nrows_skip = 52):
    cell_type_path = CELL_TYPE_FILE
    cell_types = pd.read_csv(cell_type_path)
    cell_types.columns = ["barcode_id", "cell_type_id"]
    cell_type_subdf = cell_types[cell_types["cell_type_id"] == cell_type]
    cell_barcode_list = cell_type_subdf.loc[:,"barcode_id"].to_list()
    i = 0
    for INPUT_FRAG in FRAG_LIST:
        DATA_PREFIX = DATA_PREFIX_LIST[i]
        print(INPUT_FRAG, DATA_PREFIX)
        frag = pd.read_csv(INPUT_FRAG, sep = '\t', header=None, skiprows = nrows_skip)
        frag.columns = ["chrom","start","end","cell_id","counts"]
        frag.loc[:,"cell_barcode"] = [f"{DATA_PREFIX}#{i}" for i in frag.loc[:,"cell_id"].to_list()]
        print(frag.shape)
        frag_sub = frag[frag["cell_barcode"].isin(cell_barcode_list)]
        print(frag_sub.shape)
        frag_sub.to_csv(SAVE_NAME, mode="a", sep="\t", index=None, header=False)
        del frag, frag_sub
        i += 1
    return None


#########################
#     Script running    #
#########################

def main():
    args = parser.parse_args()
    config = vars(args)
    print(config)
    CELL_TYPE = config['cell_type']
    FRAG_LIST = config['fragment_list']
    DATA_PREFIX_LIST = config['data_prefix_list']
    SAVE_NAME = config['save_name']
    CELL_TYPE_FILE = config['cell_type_file']
    GENOME_ASSEMBLY = config['genome_assembly']
    
    _ = filter_celltype(FRAG_LIST, DATA_PREFIX_LIST,
                        cell_type=CELL_TYPE, CELL_TYPE_FILE = CELL_TYPE_FILE, 
                        SAVE_NAME = SAVE_NAME) #f"{GENOME_ASSEMBLY}_{SAVE_NAME}"
    
if __name__ == '__main__':
    main()