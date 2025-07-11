#!/usr/bin/env Rscript

# Extract library statistics from barcode-UMI counts

# Load libraries ----
options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))

# Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

raw_count_file <- args[1]
linked_count_file <- args[2]
umi_cutoff <- as.numeric(args[3])
output_file <- args[4]

# Load sample annotations ----
sample_annotations <- read_csv("../annotations/sample_annotations.csv", show_col_types = F)

# Load data ----
raw_counts <- read_csv(raw_count_file, show_col_types = FALSE)
linked_counts <- read_csv(linked_count_file, show_col_types = FALSE)

# Calculate raw statistics ----
raw_stats <- raw_counts %>%
  filter(umi_count >= umi_cutoff) %>%
  summarize(
    total_barcodes = dplyr::n(), 
    total_umi_count = sum(umi_count), 
    total_read_count = sum(read_count)) %>%
    pivot_longer(everything())

# Load insert annotations ----
insert_annotations <- read_csv("../../rbp_sgrna_barcode_linkage/annotations/insert_annotations.csv", show_col_types = F) %>%
  select(insert_num, gene, sgrna)

# Calculate linked statistics ----
linked_stats <- linked_counts %>%
  filter(umi_count >= umi_cutoff) %>%
  left_join(insert_annotations, by = "insert_num") %>%
  group_by(insert_num) %>%
  summarize(n_barcodes = dplyr::n(), umi_count = sum(umi_count), read_count = sum(read_count), gene = first(gene), .groups = "drop") %>%
  summarize(
    total_genes = dplyr::n_distinct(gene),
    total_inserts = dplyr::n(),
    total_linked_read_count = sum(read_count), 
    total_linked_umi_count = sum(umi_count), 
    total_linked_barcodes = sum(n_barcodes),
    median_barcodes_per_insert = median(n_barcodes), 
    median_reads_per_insert = median(read_count), 
    median_umis_per_insert = median(umi_count)
  ) %>%
  pivot_longer(everything())

# Combine and write output ----
bind_rows(raw_stats, linked_stats) %>%
  write_csv(output_file)