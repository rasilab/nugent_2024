#!/usr/bin/env Rscript

# Create input files for MaGeCK

# Load libraries ----
suppressPackageStartupMessages(library(tidyverse))

# Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)
treatment_counts_file <- args[1]
control_counts_file <- args[2]
insert_annotation_file <- args[3]
umi_cutoff <- as.numeric(args[4])
output_file <- args[5]

# Load insert annotations ----
insert_annotations <- read_csv(insert_annotation_file, show_col_types = FALSE) %>%
  select(insert_num, gene, sgrna)

# Process counts data ----
counts <- list(
  "treatment" = read_csv(treatment_counts_file, show_col_types = FALSE),
  "control" = read_csv(control_counts_file, show_col_types = FALSE)
) %>%
  bind_rows(.id = "sample") %>%
  select(insert_num, sample, umi_count) %>%
  pivot_wider(names_from = sample, values_from = umi_count, values_fill = umi_cutoff) %>%
  left_join(insert_annotations, by = "insert_num") %>%
  select(sgrna, gene, treatment, control) %>%
  write_tsv(output_file)