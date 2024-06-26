#!/usr/bin/env python

'''Script for processing scATAC fragment files.

Usage example: 
screen
bsub -q gpuqueue -gpu - -W 4:00 -n 2 -R 'span[hosts=1] rusage[mem=128]' -Is /bin/bash
python ./hicdcplus/step2_python_save.py \
--assembly 'hg38'

'''
import pandas as pd
import argparse
import pickle
import numpy as np
from scipy.sparse import csr_matrix

parser = argparse.ArgumentParser(description="Set-up data preparations",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--assembly", help="Assembly of the Hi-C library.",
                    type=str, default="hg38")
args = parser.parse_args()
config = vars(args)
print(config)
ASSEMBLY = config['assembly']

def main():

    args = parser.parse_args()
    config = vars(args)
    print(config)
    ASSEMBLY = config['assembly']

    PTH = './intermediate'
    SAVE_PTH = './zvalue_intermediate'

    hic_dict = {}
    hic_resolution = 10000
    if ASSEMBLY in ['hg19', 'hg38']:
        chrom_length = 22
    elif ASSEMBLY in ['mm9', 'mm10']:
        chrom_length = 19
    for i in list(range(1,int(chrom_length+1))) + ['X']:
        chrom = f'chr{i}'
        print(chrom)
        a = pd.read_csv(f'{PTH}/normalized_{chrom}.txt', sep='\t')
        a.loc[:,'zvalue'] = [(i-j)/k for i,j,k in zip(a.loc[:,'counts'].to_list(),
                                                    a.loc[:,'mu'].to_list(),
                                                    a.loc[:,'sdev'].to_list())]
        a1 = a.loc[:,['start1','start2','zvalue']]
        a1.to_csv(f'{SAVE_PTH}/normalized_{chrom}_simplified.txt', sep='\t',header=None,index=None)
        del a1
        
        # save into dictionary
        hic = pd.read_table(f'{SAVE_PTH}/normalized_{chrom}_simplified.txt', header = None, sep='\t')
        start = hic[[0,1]].values.min()
        end = hic[[0,1]].values.max()+hic_resolution
        hic = hic.pivot(0,1,2).fillna(0)
        hic = hic.reindex(columns=np.arange(0,end, hic_resolution).astype(int), 
                        index=np.arange(0,end, hic_resolution).astype(int),
                        fill_value = 0)
        hic_dict[chrom] = csr_matrix(hic)

    SAVE_PTH2 = './zvalue'
    pickle.dump(hic_dict, open(f"{SAVE_PTH2}/normalized_hic_zscore_dict.p", "wb"), protocol=4)


if __name__ == '__main__':
    main()
