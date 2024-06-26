#!/bin/bash

# Usage
# bash ChromaFold/process_input/hicdc_normalization/hicdcplus_normalization.sh \
# hic_file='imr90.hic'
# resolution=10000
# assembly='hg38'

SCRIPT_DIR="/chromafold/process_input/hic_normalization/hicdcplus"

# Ensure required variables are set
if [ -z "$hic_file" ] || [ -z "$resolution" ] || [ -z "$assembly" ]; then
  echo "Error: Missing required environment variables."
  echo "Usage: hic_file='your_file.hic' resolution=10000 assembly='hg38' ./your_script.sh"
  exit 1
fi

mkdir -p ./hicdc_normalization
mkdir -p ./hicdc_normalization/zvalue_intermediate
mkdir -p ./hicdc_normalization/zvalue
mkdir -p ./hicdc_normalization/features
mkdir -p ./hicdc_normalization/intermediate

cd ./hicdc_normalization

#1. Step 1. Use HiC-DC+ to calculate normalization

# Usage
# screen
source /miniconda3/etc/profile.d/conda.sh
conda activate chromafold_env
cd /chromafold/scripts

Rscript "${SCRIPT_DIR}"/step1_hicdcplus_normalization_run.R \
"$hic_file" \
"$resolution" \
"$assembly"

# exit from R env
conda deactivate

#2. Convert calculated z-values to input Hi-C format 

# load python env
conda activate cuda111_torch
python "${SCRIPT_DIR}"/step2_python_save.py \
--assembly "$assembly"

#3. Normalized .p file should be saaved into ./hicdc_normalization/zvalue