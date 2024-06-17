#/bin/bash

# Copy from chromafold_data_preparation.sh

# Usage:
# 
# screen
# bash /chromafold/scripts/pipeline/fragment_to_input.sh

#######################################################
#       Step 0. Filter and merge fragment files       #
#######################################################

'''
This section is for when you have multiple fragment files covering multiple cell types. 

Input: fragments_1.tsv.gz, fragments_2.tsv.gz, cell_type_file.csv

Task:
- specify selected cell type name 
- select from the corresponding fragments files
- merge all fragments for the specific cell type
'''

bsub -n 2 -W 10:00 -R 'span[hosts=1] rusage[mem=256]' -Is /bin/bash

python /chromafold/scripts/pipeline/fragment_celltype_merge.py \
--cell_type selected_cell_type \
--fragment_list /data/fragments_1.tsv.gz \
/data/fragments_2.tsv.gz \
--data_prefix_list "selected_cell_type" \
--save_name /data/selected_fragments.tsv.gz \
--cell_type_file /data/cell_type_file.csv \
--genome_assembly hg38

#######################################
#       Step 1. Calculate ArchR       #
#######################################

#--------------- Define arguments ---------------#

SAVE_LOC="/chromafold/data/giles_paper"
DATA_PREFIX="mm10_giles_TransCTLII_fragments"
FRAG_LOC="/chromafold/data/giles_paper/fragments_by_celltype"
FRAG_FILE_PREFIX="${DATA_PREFIX}"
GENOME_ASSEMBLY="mm10"

#----------------- Run scripts -----------------#

mkdir -p /chromafold/data/test_input/archr_data/"${DATA_PREFIX}"
ARCHR_LOC=/chromafold/data/test_input/archr_data/"${DATA_PREFIX}"

export PATH=/chromafold/packages/samtools/bin:$PATH
export PATH=/chromafold/packages/samtools/samtools/bin:$PATH
export PATH=/setd2/packages/bedtools2/bin:$PATH
gunzip -c "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}.tsv.gz" | bgzip  > "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz.tsv.gz" # convert gzip to bgzip
sortBed -i "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz.tsv.gz" > "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv" # sort merged bgzip file
htsfile "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv" # check whether file is in bed format
bgzip "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv" # bgzip file for tabix

rm "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv.gz.tbi" # remove previously calculated .tbi file
# rm "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}.tsv.gz" # remove intermediate files
rm "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz.tsv.gz" # remove intermediate files

cd "${SAVE_LOC}"
mkdir -p atac
mkdir -p dna
mkdir -p predictions

# Copy CTCF motif score to the designated folder
# CTCF motif can be downloaded from google drive: https://drive.google.com/drive/folders/1TlfXGix2U-K1UIrOYv8gWGIiSx10dP3M?usp=sharing
cp /dna/* "${SAVE_LOC}"/dna/

# Run R to create LSI file using ArchR 
source /miniconda3/etc/profile.d/conda.sh # load conda
conda activate /environment/chromafold_env # activate conda environment for R
cd /chromafold/scripts 

Rscript /chromafold/scripts/pipeline/ArchR_preparation.R \
"${DATA_PREFIX}" \
"${ARCHR_LOC}" \
"${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv.gz" \
"${GENOME_ASSEMBLY}"
# exit

############################################
#       Step 2. Calculate tile files       #
############################################

# Run Python to create tile files for scATAC data 
python /chromafold/scripts/pipeline/scATAC_preparation.py \
--cell_type_prefix "${DATA_PREFIX}" \
--fragment_file  "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv.gz" \
--barcode_file "${ARCHR_LOC}"/archr_filtered_barcode.csv \
--lsi_file "${ARCHR_LOC}"/archr_filtered_lsi.csv \
--genome_assembly "${GENOME_ASSEMBLY}" \
--save_path "${SAVE_LOC}"

