#!/usr/bin/env Rscript

# Run ArchR to generate metacell information for running ChromaFold
#
# Installation of ArchR
# install.packages("Rcpp", dependencies = TRUE) # nolint
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools") # nolint
# install.packages("devtools") # nolint
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") # nolint
# devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories()) # nolint
# library(ArchR) # nolint
# ArchR::installExtraPackages() # nolint
#
# Usage
# screen
# bsub -n 2 -W 10:00 -R 'span[hosts=1] rusage[mem=64]' -Is /bin/bash
# bsub -q gpuqueue -gpu - -W 4:00 -n 2 -R 'span[hosts=1] rusage[mem=128]' -Is /bin/bash #or using gpu # nolint
# source /home/yangr2/miniconda3/etc/profile.d/conda.sh
# conda activate chromafold_environment
# cd /data/leslie/yangr2/chromafold/scripts
#
# Rscript /data/leslie/yangr2/chromafold/scripts/ArchR_preparation.R \
# "daniel_acute_d8_tetpos" \
# "/data/leslie/yangr2/chromafold/data/test_input/archr_data" \
# "/data/leslie/tyler/data/daniel-2022/acute-d8-tetpos/outs/fragments.tsv.gz" \
# "mm10"


library(ArchR)
library(tidyr)

args <- commandArgs(trailing = TRUE)

#####################################
#      Step 0. Initial set-ups      #
#####################################

#1. Initial set-up
sample_name <- args[1]
archr_path <- args[2]
input_path <- args[3]
genome_assembly <- args[4]

# sample_name <- "daniel_acute_d8_tetpos"  # nolint
# archr_path <- "/data/leslie/yangr2/chromafold/data/test_input/archr_data"  # nolint
# input_path <- "/data/leslie/tyler/data/daniel-2022/acute-d8-tetpos/outs/fragments.tsv.gz" # nolint
setwd(archr_path)

#2. Hyper-parameters
# addArchRGenome("mm10")  # nolint
addArchRGenome(genome_assembly)
addArchRThreads(threads = 1)
tile_mat_params <- list()
tile_mat_params$tileSize <- 500

#########################################
#      Step 1. Load fragmens files      #
#########################################

#1. Load fragments files

ArrowFiles <- createArrowFiles( # nolint
    inputFiles = input_path,
    sampleNames = sample_name,
    TileMatParams = tile_mat_params,
    filterTSS = 4, #Dont set this too high because you can always increase later
    filterFrags = 1000,
    addTileMat = TRUE,
    addGeneScoreMat = FALSE
)

print("Finished creating ArrowFiles")

#2. Create archr project
proj <- ArchRProject(
    ArrowFiles = c(paste0(archr_path, "/", sample_name, ".arrow")),
    # outputDirectory = paste0(archr_path, "/archr_out/"), #nolint
    outputDirectory = paste0(archr_path, "/"),
    copyArrows = FALSE
)

print("Finished creating Archr project")

#3. Run lsi
proj <- addIterativeLSI(
  ArchRProj = proj, useMatrix = "TileMatrix",
  name = "IterativeLSI", iterations = 2,
  clusterParams = list(
    resolution = c(0.2), sampleCells = 10000, n.start = 10
  ),
  varFeatures = 25000, dimsToUse = 1:30
)

print("Finished running lsi")

#3.1. Clustering
proj <- addClusters(
    input = proj, reducedDims = "IterativeLSI",
    method = "Seurat", name = "Clusters",
    resolution = 2, force = TRUE
)

print("Finished clustering")

#3.2. tsne
# proj <- addTSNE(
#     ArchRProj = proj, reducedDims = "IterativeLSI",
#     name = "TSNE", perplexity = 30
# )

####################################
#      Step 2. Save all files      #
####################################

#1. Save
final_bc <- rownames(proj@cellColData)
write.csv(final_bc, file = paste0(archr_path,  "/archr_filtered_barcode.csv"))
lsi <- getReducedDims(proj)
write.csv(lsi, file = paste0(archr_path, "/archr_filtered_lsi.csv"))
saveArchRProject(ArchRProj = proj, outputDirectory = archr_path, load = FALSE)
