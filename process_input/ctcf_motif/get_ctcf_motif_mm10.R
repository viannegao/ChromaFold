suppressMessages(library(AnnotationHub))
ah <- AnnotationHub()
#> snapshotDate(): 2022-10-31
query_data <- subset(ah, preparerclass == "CTCF")
# Explore the AnnotationHub object
subset(query_data, 
         genome == "mm10")
CTCF_mm10_all <- query_data[["AH104753"]]
suppressMessages(library(plyranges))
CTCF_mm10_all <- CTCF_mm10_all %>% keepStandardChromosomes() %>% sort()

#> [1] "Number of CTCF motifs at the 1e-6 threshold: 21671"
# Similarly, filter
CTCF_mm10_all_filtered <- CTCF_mm10_all %>% plyranges::filter(pvalue < 2)
print(paste("Number of CTCF motifs at the 1e-6 threshold:", length(CTCF_mm10_all)))

CTCF_mm10_all_filtered_forward <- CTCF_mm10_all_filtered[CTCF_mm10_all_filtered@strand == '+',]
CTCF_mm10_all_filtered_reverse <- CTCF_mm10_all_filtered[CTCF_mm10_all_filtered@strand == '-',]

# write R code to overlap a GRanges file with genome-wide 50bp bins.
# Load required packages
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Mmusculus.UCSC.mm10)

mm10_seqinfo <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)

# Extract the chromosome names and lengths from the hg38 genome assembly
chromosomes <- c(mm10_seqinfo@seqnames[0:21])
chrom_lengths <- c(mm10_seqinfo@seqlengths[0:21])
names(chrom_lengths) <- chromosomes
bin_size <- 50
for (chrom in chromosomes){
  bins <- GRanges(seqnames = Rle(rep(chrom, floor(chrom_lengths[chrom] / bin_size)+1)),
                  ranges = IRanges(start = seq(0, floor(chrom_lengths[chrom] / bin_size)) * bin_size,
                                   end = seq(0, floor(chrom_lengths[chrom] / bin_size)) * bin_size + bin_size))
  
  # Find the overlapping bins for each element in the GRanges file
  overlapping_bins <- findOverlaps(subject = bins, query = CTCF_mm10_all_filtered_forward, minoverlap = 10, ignore.strand=T)
  scores <- CTCF_mm10_all_filtered_forward[queryHits(overlapping_bins)]$score
  
  # Create a binary vector of whether there is an overlap
  overlap_vector <- integer(length(bins))
  overlap_vector[subjectHits(overlapping_bins)] <- scores
  
  # Save the overlap vector to a binary file
  write.csv(overlap_vector,file=paste0('~/data/dna/mm10_ctcf_motif_forward_', chrom,'.csv'),row.names=F)
  
}

for (chrom in chromosomes){
  bins <- GRanges(seqnames = Rle(rep(chrom, floor(chrom_lengths[chrom] / bin_size)+1)),
                  ranges = IRanges(start = seq(0, floor(chrom_lengths[chrom] / bin_size)) * bin_size,
                                   end = seq(0, floor(chrom_lengths[chrom] / bin_size)) * bin_size + bin_size))
  
  # Find the overlapping bins for each element in the GRanges file
  overlapping_bins <- findOverlaps(subject = bins, query = CTCF_mm10_all_filtered_reverse, minoverlap = 10, ignore.strand=T)
  scores <- CTCF_mm10_all_filtered_reverse[queryHits(overlapping_bins)]$score
  
  # Create a binary vector of whether there is an overlap
  overlap_vector <- integer(length(bins))
  overlap_vector[subjectHits(overlapping_bins)] <- scores
  
  # Save the overlap vector to a binary file
  write.csv(overlap_vector,file=paste0('~/data/dna/mm10_ctcf_motif_reverse_', chrom,'.csv'),row.names=F)
  
}

