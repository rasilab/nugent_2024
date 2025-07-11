#!/usr/bin/env Rscript

# Compare barcode partitions
# Converted from compare_barcode_partitions.ipynb

# Load libraries
options(warn = -1)

suppressPackageStartupMessages({
library(tidyverse)
library(rasilabRtemplates)
})

theme_set(theme_rasilab() + 
 theme(
  axis.line = element_line(color = "grey"), 
 axis.title.y = element_text(margin = margin(r=10)),
 axis.title.x = element_text(margin = margin(t=10))
))

fdr_cutoff <- 0.05
sgrna_cutoff <- 3

# Load MaGeCK gene splicing data
mageck_gene <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc"), str_detect(file, "retained|skipped")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names() %>%
  select(-neg_lfc) %>%
  rename(lfc = pos_lfc, gene = id) %>%
  filter(!str_detect(gene, "MCHERRY|PURO")) %>%
  mutate(pos_fdr = if_else(is.na(pos_fdr), 1, pos_fdr)) %>%
  arrange(pos_rank) %>% 
  filter(str_detect(treatment, "day7")) %>%
  mutate(day = str_replace(str_extract(treatment, "day\\d+"), "day", "day ")) %>%
  mutate(isoform = str_extract(treatment, "(e2|i1|i2)_.+"))

gene_hits <- mageck_gene %>% 
  filter(pos_fdr <= fdr_cutoff & pos_goodsgrna >= sgrna_cutoff) %>% 
  distinct(gene) %>% 
  print()

# Load MaGeCK splicing partitioned gene data
mageck_gene_random_partition <- list.files("../data/mageck_random_partition/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc"), str_detect(file, "retained|skipped")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck_random_partition//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names() %>%
  select(-neg_lfc) %>%
  rename(lfc = pos_lfc, gene = id) %>%
  filter(!str_detect(gene, "MCHERRY|PURO")) %>%
  arrange(pos_rank) %>% 
  filter(str_detect(treatment, "day7")) %>%
  mutate(day = str_replace(str_extract(treatment, "day\\d+"), "day", "day ")) %>%
  mutate(isoform = str_extract(treatment, "(e2|i1|i2)_.+"))

# Compare splicing barcode partitions
plot_data <- mageck_gene_random_partition %>% 
  separate(gene, c("gene", "barcode_group"), sep = "_") %>%
  mutate(isoform = str_replace(isoform, "_retained|_skipped", "")) %>% 
  inner_join(gene_hits, by = "gene") %>%
  select(isoform, lfc, gene, barcode_group) %>% 
  pivot_wider(names_from = barcode_group, values_from = lfc) %>%
  write_csv("../../../../source_data/figure_s3b.csv")

correlation  <- plot_data %>% 
  group_by(isoform) %>% 
  nest() %>% 
  mutate(cor = map(data, function(df) cor.test(~ A + B, data=df, method = "pearson", exact = F))) %>% 
  mutate(cor = map(cor, broom::tidy)) %>%
  select(-data) %>% 
  unnest(cor) 

plot_data %>% 
  ggplot(aes(x = A, y = B)) +
  facet_wrap(~fct_rev(isoform)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point() +
  geom_text(data = correlation, aes(label = paste0("r = ", round(estimate, 2))), x = -2, y = 5, hjust = 0, vjust = 0, size = 4) +
  labs(x = "Barcode Set A", y = "Barcode Set B")

ggsave("../figures/mageck_splicing_scatter.pdf", width = 4, height = 2)

# Heat map of splicing screen correlation coefficients
corr_data <- mageck_gene_random_partition %>% 
  separate(gene, c("gene", "barcode_group"), sep = "_") %>%
  inner_join(gene_hits, by = "gene") %>%
  mutate(isoform = str_replace(isoform, "_retained|_skipped", "")) %>% 
  select(isoform, lfc, gene, barcode_group) %>% 
  pivot_wider(names_from = barcode_group, values_from = lfc) %>% 
  rename(sample_name = isoform)

plot_data <- corr_data %>% 
  distinct(sample_name) %>% 
  rename(sample1 = sample_name) %>% 
  mutate(sample2 = sample1) %>% 
  expand(sample1, sample2) %>% 
  mutate(r = map2_dbl(sample1, sample2, function(x,y) {
    xvec <- corr_data %>% filter(sample_name == x) %>% select(gene, A)
    yvec <- corr_data %>% filter(sample_name == y) %>% select(gene, B)
    joint <- inner_join(xvec, yvec, by = "gene")
    cor.test(joint$A, joint$B, method = "pearson")$estimate
  })) %>%
  write_csv("../../../../source_data/figure_s3c.csv")

plot_data %>% 
  ggplot(aes(x = sample1, y = sample2, fill = r)) +
  geom_tile() +
  geom_text(aes(label = round(r, 2)), size = 3) +
  scale_fill_gradient2(low = "#132B43", high = "#56B1F7") +
  theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  legend.text = element_text(size = 8),
  ) +
  labs(x = "Barcode Set A", y = "Barcode Set B")

ggsave("../figures/mageck_splicing_correlation.pdf", width = 2.8, height = 2, units = "in")

# Load MaGeCK fitness data
mageck_gene <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc.+bc1.+(total|grna).+ntc.+bc1.+(total|grna)")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names() %>%
  select(-pos_lfc) %>%
  rename(lfc = neg_lfc, gene = id) %>%
  arrange(neg_rank)

# Load MaGeCK fitness partitioned gene data
mageck_gene_random_partition <- list.files("../data/mageck_random_partition/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc.+bc1.+(total|grna).+ntc.+bc1.+(total|grna)"), !str_detect(file, "total.+grna")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck_random_partition//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names() %>%
  select(-pos_lfc) %>%
  rename(lfc = neg_lfc, gene = id) %>%
  filter(!str_detect(gene, "MCHERRY|PURO")) %>%
  arrange(neg_rank) %>% 
  filter(str_detect(treatment, "day5|day13|day21")) %>%
  mutate(day = str_replace(str_extract(treatment, "day\\d+"), "day", "day "))

# Compare fitness data for barcode partitions
plot_data <- mageck_gene_random_partition %>% 
  separate(gene, c("gene", "barcode_group"), sep = "_") %>%
  mutate(type = str_extract(treatment, "grna|total")) %>% 
  select(type, day, lfc, gene, barcode_group) %>% 
  pivot_wider(names_from = barcode_group, values_from = lfc) %>% 
  mutate(day = fct_relevel(day, "day 5", "day 13", "day 21")) %>% 
  mutate(type = fct_relevel(type, "total", "grna")) %>%
  write_csv("../../../../source_data/figure_s1f.csv")

correlation  <- plot_data %>% 
  group_by(type, day) %>% 
  nest() %>% 
  mutate(cor = map(data, function(df) cor.test(~ A + B, data=df, method = "pearson", exact = F))) %>% 
  mutate(cor = map(cor, broom::tidy)) %>%
  select(-data) %>% 
  unnest(cor) 

plot_data %>% 
  ggplot(aes(x = A, y = B)) +
  facet_wrap(~type + day, ncol = 3, scales = "fixed") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_text(data = correlation, aes(label = paste0("r = ", round(estimate, 2))), x = -1, y = -6, hjust = 0, vjust = 0, size = 4) +
  labs(x = "log2 fold change w.r.t day 1\nBarcode Set A", y = "log2 fold change w.r.t day 1\nBarcode Set B") +
  scale_x_continuous(limits = c(-6, 4)) +
  scale_y_continuous(limits = c(-6, 4))

ggsave("../figures/mageck_fitness_compare_barcode_scatter.pdf", width = 4.5, height = 4.5)

print("Analysis complete!")