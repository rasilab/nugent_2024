options(warn = -1)
suppressPackageStartupMessages({
  library(Biostrings)
  library(plyranges)
  library(tidyverse)
  library(GenomicAlignments)
})

args <- commandArgs(trailingOnly = TRUE)

# intron_file <- "../data/mane/MANE.GRCh38.v1.3.ensembl_genomic.tsv"
# alignment_file <- "../data/alignments//genome_mane_183p1/Aligned.sortedByCoord.out.bam"
intron_file <- args[1]
alignment_file <- args[2]
output_file <- args[3]

sj <- read_tsv(intron_file, show_col_types = F) %>% 
  GRanges() %>% 
  setNames(seq(1, length(.))) %>% 
  filter(junction_type == "annotated") 

aln <- readGAlignmentPairs(alignment_file) %>% 
  invertStrand()

sj %>% 
  findOverlaps(aln, ., minoverlap = 10L, ignore.strand = F) %>% 
  subjectHits() %>% 
  sj[.] %>%
  names() %>%
  table() %>%
  as.data.frame() %>%
  setNames(c("intron", "count")) %>%
  write_csv(output_file)
