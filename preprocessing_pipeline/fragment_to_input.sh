#/bin/bash

# Copy from chromafold_data_preparation.sh

# Usage:
# 
# screen
# bsub -n 2 -W 10:00 -R 'span[hosts=1] rusage[mem=256]' -Is /bin/bash
# bash /data/leslie/yangr2/chromafold/scripts/pipeline/fragment_to_input.sh

#######################################################
#       Step 0. Filter and merge fragment files       #
#######################################################

bsub -n 2 -W 10:00 -R 'span[hosts=1] rusage[mem=256]' -Is /bin/bash
source /home/yangr2/dnabert_environment/bin/activate

# Already run: 'Effector', 'Eff-like.I', 'Eff-like.II', 
#              'Mem', 'Mem.CTL', 
#              'Exh.Pre', 'Exh.Trans', 'Exh.Intermediate', 'Exh.KLR', 'Exh.Progenitor', 'Exh.Terminal', 
#              'Naive', 'CTL', 'MP'
#              'Trans.Mem', 'Trans.CTL.I', 'Trans.CTL.II'

python /data/leslie/yangr2/chromafold/scripts/pipeline/fragment_celltype_merge.py \
--cell_type Exh.Trans \
--fragment_list /data/leslie/tyler/data/giles-2022/SRR18505566/outs/fragments.tsv.gz \
/data/leslie/tyler/data/giles-2022/SRR18505565/outs/fragments.tsv.gz \
/data/leslie/tyler/data/giles-2022/SRR18505564/outs/fragments.tsv.gz \
/data/leslie/tyler/data/giles-2022/SRR18505563/outs/fragments.tsv.gz \
/data/leslie/tyler/data/giles-2022/SRR18505562/outs/fragments.tsv.gz \
/data/leslie/tyler/data/giles-2022/SRR18505561/outs/fragments.tsv.gz \
/data/leslie/tyler/data/giles-2022/SRR18505560/outs/fragments.tsv.gz \
--data_prefix_list "giles_naive_d0" "giles_acute_d8" "giles_acute_d15" "giles_acute_d30" \
"giles_chronic_d8" "giles_chronic_d15" "giles_chronic_d30" \
--save_name /data/leslie/yangr2/chromafold/data/giles_paper/fragments_by_celltype/mm10_giles_ExhTrans_fragments.tsv.gz \
--cell_type_file /data/leslie/tyler/t-cell-exhaustion/cellspace-experiments3/arrow-giles/cell-md-from-paper/cell-md-fig1.csv \
--genome_assembly mm10

#######################################
#       Step 1. Calculate ArchR       #
#######################################

# Needs to run: 
#              
#              
#              
#              
# Already run: 'Effector', 'EfflikeI', 'EfflikeII',
#              'Mem', 'MemCTL',
#              'ExhPre', 'ExhTrans', 'ExhIntermediate', 'ExhKLR', 'ExhProgenitor', 'ExhTerminal', 
#              'Naive', 'CTL', 'MP'
#              'TransMem', 'TransCTLI', 'TransCTLII'

#--------------- Define arguments ---------------#

SAVE_LOC="/data/leslie/yangr2/chromafold/data/giles_paper"
DATA_PREFIX="mm10_giles_TransCTLII_fragments"
FRAG_LOC="/data/leslie/yangr2/chromafold/data/giles_paper/fragments_by_celltype"
FRAG_FILE_PREFIX="${DATA_PREFIX}"
GENOME_ASSEMBLY="mm10"

#----------------- Run scripts -----------------#

mkdir -p /data/leslie/yangr2/chromafold/data/test_input/archr_data/"${DATA_PREFIX}"
ARCHR_LOC=/data/leslie/yangr2/chromafold/data/test_input/archr_data/"${DATA_PREFIX}"

export PATH=/data/leslie/yangr2/chromafold/packages/samtools/bin:$PATH
export PATH=/data/leslie/yangr2/chromafold/packages/samtools/samtools/bin:$PATH
export PATH=/data/leslie/yangr2/setd2/packages/bedtools2/bin:$PATH
gunzip -c "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}.tsv.gz" | bgzip  > "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz.tsv.gz" # convert gzip to bgzip
sortBed -i "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz.tsv.gz" > "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv" # sort merged bgzip file
htsfile "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv" # check whether file is in bed format
bgzip "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv" # bgzip file for tabix

# count number of unique cells
# zcat mm10_daniel_acute_exhterminal_fragments.tsv.gz |cut -f 4 | sort | uniq | wc -l

rm "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv.gz.tbi" # remove previously calculated .tbi file
# rm "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}.tsv.gz" # remove intermediate files
rm "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz.tsv.gz" # remove intermediate files

cd "${SAVE_LOC}"
mkdir -p atac
mkdir -p dna
mkdir -p predictions

# Copy CTCF motif score to the designated folder
cp /data/leslie/shared/forRui/fromYuchen/dna/* "${SAVE_LOC}"/dna/

# Run R to create LSI file using ArchR 
source /home/yangr2/miniconda3/etc/profile.d/conda.sh # load conda
conda activate /data/leslie/yangr2/seqmodel/environment/chromafold_environment # activate conda environment for R
cd /data/leslie/yangr2/chromafold/scripts 

Rscript /data/leslie/yangr2/chromafold/scripts/pipeline/ArchR_preparation.R \
"${DATA_PREFIX}" \
"${ARCHR_LOC}" \
"${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv.gz" \
"${GENOME_ASSEMBLY}"
# exit

############################################
#       Step 2. Calculate tile files       #
############################################

# Run Python to create tile files for scATAC data 
source /home/yangr2/dnabert_environment/bin/activate # activate python environment
python /data/leslie/yangr2/chromafold/scripts/pipeline/scATAC_preparation.py \
--cell_type_prefix "${DATA_PREFIX}" \
--fragment_file  "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv.gz" \
--barcode_file "${ARCHR_LOC}"/archr_filtered_barcode.csv \
--lsi_file "${ARCHR_LOC}"/archr_filtered_lsi.csv \
--genome_assembly "${GENOME_ASSEMBLY}" \
--save_path "${SAVE_LOC}"

