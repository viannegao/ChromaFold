suppressMessages(library(AnnotationHub))
ah <- AnnotationHub()
#> snapshotDate(): 2022-10-31
query_data <- subset(ah, preparerclass == "CTCF")
# Explore the AnnotationHub object
subset(query_data, species == "Homo sapiens" & 
         genome == "hg38")

CTCF_hg38_all <- query_data[["AH104727"]]
CTCF_hg38 <- query_data[["AH104729"]]
suppressMessages(library(plyranges))
CTCF_hg38_all <- CTCF_hg38_all %>% keepStandardChromosomes() %>% sort()
CTCF_hg38 <- CTCF_hg38 %>% keepStandardChromosomes() %>% sort()

# Check length before filtering
print(paste("Number of CTCF motifs at the default 1e-4 threshold:", length(CTCF_hg38)))
#> [1] "Number of CTCF motifs at the default 1e-4 threshold: 887980"
# Filter and check length after filtering
CTCF_hg38_filtered <- CTCF_hg38 %>% plyranges::filter(pvalue < 2)
print(paste("Number of CTCF motifs at the 1e-6 threshold:", length(CTCF_hg38_filtered)))
#> [1] "Number of CTCF motifs at the 1e-6 threshold: 21671"
# Similarly, filter
CTCF_hg38_all_filtered <- CTCF_hg38_all %>% plyranges::filter(pvalue < 2)

CTCF_hg38_all_filtered_forward <- CTCF_hg38_all_filtered[CTCF_hg38_all_filtered@strand == '+',]
CTCF_hg38_all_filtered_reverse <- CTCF_hg38_all_filtered[CTCF_hg38_all_filtered@strand == '-',]

# write R code to overlap a GRanges file with genome-wide 50bp bins.
# Load required packages
library(GenomicRanges)

# Load the hg38 genome assembly
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Extract the chromosome names and lengths from the hg38 genome assembly
chromosomes <- unique(seqlevels(txdb))[0:25]
chrom_lengths <- seqlengths(txdb)[chromosomes][0:25]


bin_size <- 50
for (chrom in chromosomes){
  bins <- GRanges(seqnames = Rle(rep(chrom, floor(chrom_lengths[chrom] / bin_size)+1)),
                  ranges = IRanges(start = seq(0, floor(chrom_lengths[chrom] / bin_size)) * bin_size,
                                   end = seq(0, floor(chrom_lengths[chrom] / bin_size)) * bin_size + bin_size))
  
  # Find the overlapping bins for each element in the GRanges file
  overlapping_bins <- findOverlaps(subject = bins, query = CTCF_hg38_all_filtered_forward, minoverlap = 10, ignore.strand=T)
  scores <- CTCF_hg38_all_filtered_forward[queryHits(overlapping_bins)]$score
  
  # Create a binary vector of whether there is an overlap
  overlap_vector <- integer(length(bins))
  overlap_vector[subjectHits(overlapping_bins)] <- scores
  
  # Save the overlap vector to a binary file
  write.csv(overlap_vector,file=paste0('~/data/dna/hg38_ctcf_motif_forward_', chrom,'.csv'),row.names=F)
  
}

for (chrom in chromosomes){
  bins <- GRanges(seqnames = Rle(rep(chrom, floor(chrom_lengths[chrom] / bin_size)+1)),
                  ranges = IRanges(start = seq(0, floor(chrom_lengths[chrom] / bin_size)) * bin_size,
                                   end = seq(0, floor(chrom_lengths[chrom] / bin_size)) * bin_size + bin_size))
  
  # Find the overlapping bins for each element in the GRanges file
  overlapping_bins <- findOverlaps(subject = bins, query = CTCF_hg38_all_filtered_reverse, minoverlap = 10, ignore.strand=T)
  scores <- CTCF_hg38_all_filtered_reverse[queryHits(overlapping_bins)]$score
  
  # Create a binary vector of whether there is an overlap
  overlap_vector <- integer(length(bins))
  overlap_vector[subjectHits(overlapping_bins)] <- scores
  
  # Save the overlap vector to a binary file
  write.csv(overlap_vector,file=paste0('~/data/dna/hg38_ctcf_motif_reverse_', chrom,'.csv'),row.names=F)
  
}


