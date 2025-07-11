#!/usr/bin/env Rscript

# Extract only reads with right subpool barcodes (upto Hamming distance of 2) and add up UMIs and reads per barcode

# Load libraries ----
options(warn = -1)
suppressPackageStartupMessages(library(tidyverse))

# Parse command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

subpool_barcode <- args[1]
sample_id <- args[2]
linked_count_file <- args[3]
output_file <- args[4]

# Helper functions ----
generate_hamming_variants <- function(barcode, characters) {
  n <- nchar(barcode)
  variants <- tibble(pos = 1:n, orig_char = strsplit(barcode, "")[[1]])

  variants %>%
    expand(nesting(pos, orig_char), new_char = characters) %>%
    filter(orig_char != new_char) %>%
    mutate(new_barcode = map2_chr(pos, new_char, ~{
      str_c(str_sub(barcode, 1, .x - 1), .y, str_sub(barcode, .x + 1, n))
    })) %>%
    pull(new_barcode)
}

# Generate barcode variants ----
characters <- c("A", "T", "G", "C")

subpool_barcodes <- subpool_barcode %>%
  as_tibble_col("subpool_barcode") %>%
  mutate(subpool_barcode1 = map(subpool_barcode, ~(generate_hamming_variants(.x,characters)))) %>%
  unnest(subpool_barcode1) %>%
  mutate(subpool_barcode2 = map(subpool_barcode1, ~(generate_hamming_variants(.x,characters)))) %>%
  unnest(subpool_barcode2) %>%
  distinct(subpool_barcode2, .keep_all = T) %>%
  mutate(subpool_barcode = subpool_barcode2) %>%
  select(subpool_barcode)

# Process counts and write output ----
counts <- read_csv(linked_count_file, show_col_types = F) %>%
  right_join(subpool_barcodes, by = "subpool_barcode") %>%
  group_by(insert_num, barcode_num) %>%
  summarise(umi_count = sum(umi_count), read_count = sum(read_count), .groups = "drop") %>%
  drop_na() %>%
  arrange(insert_num, desc(umi_count)) %>%
  mutate(barcode_group = if_else(barcode_num %% 2 == 1, "A", "B")) %>%
  write_csv(output_file)