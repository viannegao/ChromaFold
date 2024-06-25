#hicdc+ for normalization
#
# Usage
# screen
# bsub -n 2 -W 10:00 -R 'span[hosts=1] rusage[mem=128]' -Is /bin/bash
# source /miniconda3/etc/profile.d/conda.sh
# conda activate chromafold_env
# cd ./hic_normalization
#
# Rscript ./hicdcplus/step1_hicdcplus_normalization_run.R \
# imr90.hic \
# 10000 \
# 'hg38'

library(HiCDCPlus)

args <- commandArgs(trailing = TRUE)

#####################################
#      Step 0. Initial set-ups      #
#####################################

#1. Initial set-up
hicfile_path <- args[1] #location of the raw .hic file to be normalized
resolution <- args[2] #resolution of Hi-C
assembly <- args[3]

#2. Variable set-ups

if (assembly %in% c("mm9", "mm10")) {
  species <- "Mmusculus"
  chrom_length <- 19
} else if (assembly %in% c("hg19", "hg38")) {
  species <- "Hsapiens"
  chrom_length <- 22
} else {
  species <- "Unknown"
  print("Species not match")
}

outdir <- "./features/"
outpth <- "./intermediate/"

########################################
#      Step 1. Construct features      #
########################################

#Step 1. construct features (resolution, chrom specific)
chromosomes <- paste("chr", seq(1, chrom_length, 1), sep = "")
chromosomes[[chrom_length+1]] <- "chrX"
construct_features(output_path = paste0(outdir, assembly, "_", as.integer(resolution/1000),"kb_GATC_GANTC", sep = ""), # nolint
gen = species, gen_ver = assembly, sig = c("GATC", "GANTC"),
bin_type = "Bins-uniform", binsize = resolution, chrs = chromosomes)

#Step 2. generate gi_list instance

set.seed(1010)

#generate gi_list instance
gi_list <- generate_bintolen_gi_list(
  bintolen_path = paste0(outdir, "/", assembly,
  "_", as.integer(resolution / 1000), "kb_GATC_GANTC_bintolen.txt.gz"))
#add .hic counts
gi_list <- add_hic_counts(gi_list, hic_path = hicfile_path, chrs = names(gi_list)) #nolint

#############################################################
#      Step 2. Modeling and generate normalized values      #
#############################################################

#expand features for modeling
gi_list <- expand_1D_features(gi_list)
#run HiC-DC+ on 2 cores
gi_list <- HiCDCPlus_parallel(gi_list, ncore = 2)
head(gi_list$chr1)

##################################
#      Step 3. Save results      #
##################################

#write normalized counts (observed/expected) to a .hic file
hicdc2hic(gi_list,
hicfile = paste0(outpth, "normalized.hic"),
          mode = "zvalue", gen_ver = assembly)
#write results to a text file
gi_list_write(gi_list,
fname = paste0(outpth, "normalized.txt.gz"))

for (chrom in names(gi_list)) {
  print(chrom)
  # Construct the file name
  file_name <- paste0(outpth, "normalized_", chrom, ".txt")
  # Extract the data for the current chromosome
  data <- gi_list[[chrom]]
  # Save the data to a text file
  write.table(data, file = file_name, sep = "\t",
  row.names = FALSE, col.names = TRUE, quote = FALSE)
}
