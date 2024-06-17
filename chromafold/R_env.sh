#/bin/bash

mkdir -p /chromafold/packages
cd /chromafold/packages
# /miniconda3/condabin/conda install -c r r-base

#1. Create conda environment 
cd /chromafold/packages
/miniconda3/bin/conda create -n chromafold_env

#2. activate conda env and install essential packages 
source /miniconda3/etc/profile.d/conda.sh
conda init bash
conda activate chromafold_env
conda install -c r r-essentials

#3. going into R, and install packages
# https://github.com/GreenleafLab/ArchR

R

install.packages("Rcpp", dependencies = TRUE)
install.packages('git2r')

q()

# exit R, in bash
conda install conda-forge::r-rcpp
conda install -c conda-forge binutils
conda install -c conda-forge libgit2
conda install -c conda-forge r-gert
conda install -c conda-forge r-usethis
conda install -c conda-forge r-devtools
conda install -c conda-forge r-biocmanager

conda install conda-forge::r-rcpp
conda install -c conda-forge r-devtools
conda install bioconda::bioconductor-chromvar
conda install bioconda::bioconductor-dirichlemialtinom
conda install bioconda::bioconductor-motifmatchrols
conda install -c bioconda bioconductor-tfbstools

#4. Install ArchR (in R, the R installed from conda r, not independent R)

R

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

Sys.setenv(CONDA_BUILD_SYSROOT="/")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

library(ArchR)
ArchR::installExtraPackages()
