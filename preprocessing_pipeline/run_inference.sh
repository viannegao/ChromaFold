#/bin/bash

screen
bsub -q gpuqueue -gpu - -W 4:00 -n 2 -R 'span[hosts=1] rusage[mem=256]' -Is /bin/bash
source /home/yangr2/dnabert_environment/bin/activate

# Need to run: 
#              
#              
#              
#              
# Already run: 'Effector', 'EfflikeI', 'EfflikeII',
#              'Mem', 'MemCTL', 
#              'ExhPre', 'ExhTrans', 'ExhIntermediate', 'ExhKLR', 'ExhProgenitor', 'ExhTerminal', 
#              'Naive', 'CTL', 'MP'
#              'TransMem', 'TransCTLI', 'TransCTLII'

python /data/leslie/yangr2/chromafold/ChromaFold/chromafold/inference.py \
--data-path /data/leslie/yangr2/chromafold/data/giles_paper \
--save-path /data/leslie/yangr2/chromafold/data/giles_paper/predictions \
-ct mm10_giles_ExhTrans_fragments \
--model-path /data/leslie/yangr2/chromafold/ChromaFold/checkpoints/chromafold_CTCFmotif.pth.tar \
-chrom 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 \
-offset -2000000 \
--genome mm10
