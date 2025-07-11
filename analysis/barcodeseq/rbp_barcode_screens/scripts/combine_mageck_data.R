#!/usr/bin/env Rscript

# Combine MaGeCK tables for supplementary table generation

# Load libraries ----
library(tidyverse)

# Get a list of all tsv files in the subfolders ----
gene_files <- list.files(path = "../data/mageck/", pattern = "mageck.gene_summary.tsv$", recursive = TRUE, full.names = TRUE)
sgrna_files <- list.files(path = "../data/mageck/", pattern = "mageck.sgrna_summary.tsv$", recursive = TRUE, full.names = TRUE)

# Read and concatenate the tsv files ----
gene_data <- read_tsv(gene_files, col_types = cols(.default = "c"), id = "sample_name") %>%
  mutate(sample_name = str_extract(sample_name, "[^/]+(?=/mageck.gene_summary)"))

sgrna_data <- read_tsv(sgrna_files, col_types = cols(.default = "c"), id = "sample_name" ) %>%
  mutate(sample_name = str_extract(sample_name, "[^/]+(?=/mageck.sgrna_summary)"))

# Write output files ----
write_csv(gene_data, "../data/mageck/gene_summary_table.csv.gz")
write_csv(sgrna_data, "../data/mageck/sgrna_summary_table.csv.gz")