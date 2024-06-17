#/bin/bash

cd /chromafold/packages
/miniconda3/condabin/conda install -c r r-base

# https://unix.stackexchange.com/questions/149451/install-r-in-my-own-directory 
# wget http://cran.rstudio.com/src/base/R-3/R-3.6.3.tar.gz
# tar xvf R-3.6.3.tar.gz
# cd R-3.6.3
# ./configure --prefix=/chromafold/packages/R-4.2.3
# make && make install
# export PATH=/chromafold/packages/R-3.6.3/bin:$PATH

#1. Create conda environment 
cd /chromafold/packages
/miniconda3/bin/conda create -n chromafold_env

#2. activate conda env and install essential packages 
# related links: https://github.com/conda/conda/issues/7980, 
# https://stackoverflow.com/questions/69187142/solving-environment-failed-with-current-repodata-json-will-retry-with-next-rep
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

# exit R, in bash: https://github.com/r-lib/gert/issues/140 
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

# start R
# export LD_LIBRARY_PATH=/home/libgit2-1.1.0/build:${LD_LIBRARY_PATH}

#4. Install ArchR (in R, the R installed from conda r, not independent R)

R

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

Sys.setenv(CONDA_BUILD_SYSROOT="/")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

library(ArchR)
ArchR::installExtraPackages()
