#!/usr/bin/env Rscript

# Analyze polysome ReLiC data
# Converted from plot_polysome_relic_data.ipynb

# Load libraries
options(warn = -1)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(tidyverse)
  library(rasilabRtemplates)
  library(ggpubr)
})

cbPalette_12 <- c(
  "#DDCC78", "#CC6677", "#6699CC", "#661100", "#117733", "#999933",
  "#332288", "#AA4499", "#44AA99", "#882255", "#88CCEE", "#999999"
)

theme_set(theme_rasilab() +
  theme(
    axis.line = element_line(color = "grey"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  ))

fdr_cutoff <- 0.05

set.seed(111)

# Load ribosome biogenesis annotations
orgdb <- org.Hs.eg.db

ribi <- AnnotationDbi::select(orgdb, keytype = "GOALL", keys = c("GO:0042273", "GO:0042274"), columns = c("SYMBOL")) %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  distinct(goall, symbol) %>%
  dplyr::rename(gene = symbol, go_id = goall) %>%
  mutate(annotation = if_else(go_id == "GO:0042273", "60S biogenesis", "40S biogenesis")) %>%
  filter(!str_detect(gene, "^RPL|^RPS|GNB2L1|RACK1"))

# Read essential gene list
essential_genes <- read_csv("../annotations/DepMapCRISPRInferredCommonEssentials_23Q2.csv", show_col_types = F) %>% 
  separate(Essentials, c("gene", "num")) %>% 
  select(gene) %>%
  mutate(essential = TRUE)

# Load MaGeCK gene-level results
mageck_gene <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc.+polysome")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(rep = str_extract(sample_name, "rep1|rep2")) %>% 
  mutate(treatment = str_extract(treatment, "heavy|light|mono|supernatant")) %>%
  mutate(control = str_extract(control, "mono|supernatant")) %>%
  filter(!is.na(treatment), !is.na(control), !is.na(rep)) %>% 
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names() %>%
  select(-neg_lfc) %>%
  rename(lfc = pos_lfc, gene = id) %>%
  mutate(p_value = if_else(lfc > 0, pos_p_value, neg_p_value)) %>%
  mutate(fdr = if_else(lfc > 0, pos_fdr, neg_fdr)) %>%
  mutate(goodsgrna = if_else(lfc > 0, pos_goodsgrna, neg_goodsgrna)) %>%
  mutate(fdr = if_else(is.na(fdr), 1, fdr))

# Load MaGeCK sgRNA-level results
mageck_sgrna <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "sgrna_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc.+polysome")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(rep = str_extract(sample_name, "rep1|rep2")) %>% 
  mutate(treatment = str_extract(treatment, "heavy|light|mono|supernatant")) %>%
  mutate(control = str_extract(control, "mono|supernatant")) %>%
  filter(!is.na(treatment), !is.na(control), !is.na(rep)) %>% 
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names()

# How many hits overlap between the heavy-mono and heavy-superntant comparisons?
plot_data <- mageck_gene %>% 
  filter(treatment == "heavy", fdr < fdr_cutoff) %>% 
  mutate(rep = as.integer(str_extract(rep, ".$"))) %>%
  group_by(gene, control) %>% 
  summarize(total_rep  = sum(rep), .groups = "drop") %>% 
  filter(total_rep == 3)

# Look at genes that selectively affect heavy polysome / supernatant ratios
interesting_supernatant_genes <- c("XRN1", "DDX6", "SNTB1", "EIF3A", "EIF2S1", "EIF4E", "EIF5B", "EIF5A")

mageck_gene  %>% 
  filter(gene %in% interesting_supernatant_genes, treatment == "heavy") %>% 
  select(-sample_name, -matches("score")) %>% 
  arrange(gene) %>% 
  select(control, rep, gene, lfc) %>% 
  mutate(rep = str_replace(rep, "rep", "")) %>% 
  pivot_wider(names_from = control, values_from = lfc) %>%
  write_csv("../../../../source_data/figure_s2d.csv") %>%
  ggplot(aes(x = supernatant, y = mono, color = gene, shape = rep)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 3)  +
  scale_x_continuous(limits = c(-1.65, 0.4)) +
  scale_y_continuous(limits = c(-1.25, 1)) +
  labs(x = "log2 heavy polysome\n/ supernatant", y = "log2 heavy polysome\n/ monosome", rep = "Replicate") 

ggsave("../figures/interesting_supernatant_genes.pdf", width = 4.2, height = 2.8)

# Calculate correlation between replicates
corr_data <- mageck_sgrna %>% 
  select(sgrna, treatment, control, rep, lfc) %>% 
  pivot_wider(names_from = rep, values_from = lfc) %>% 
  filter(!is.na(rep1), !is.na(rep2)) %>% 
  unite(sample_name, treatment, control, sep = "_vs_")

# Scatter plot of correlation between replicates
corr_data %>% 
  drop_na() %>% 
  filter(!str_detect(sample_name, "supernatant")) %>%
  mutate(sample_name = str_extract(sample_name, "^[:alpha:]+")) %>%
  write_csv("../../../../source_data/figure_s2a_scatter.csv") %>%
  ggplot(aes(x = rep1, y = rep2)) +
  facet_wrap(~ sample_name, ncol = 1, scales = "free") +
  geom_point(color = "black", alpha = 0.1, stroke = NA) +
  labs(x = "Replicate 1", y = "Replicate 2") +
  scale_x_continuous(limits = c(-6, 3), breaks = c(-6, -3, 0, 3)) +
  scale_y_continuous(limits = c(-6, 3), breaks = c(-6, -3, 0, 3)) 

ggsave("../figures/mageck_polysome_correlation.pdf", width = 1.5, height = 3)

# Density plot of polysome to monosome ratios
corr_data %>% 
  filter(str_detect(sample_name, "mono$")) %>% 
  drop_na() %>% 
  mutate(sample_name = str_extract(sample_name, "^[:alpha:]+")) %>%
  write_csv("../../../../source_data/figure_s2a_histo.csv") %>%
  ggplot(aes(x = rep1)) +
  facet_wrap(~ sample_name, ncol = 1) +
  geom_area(aes(y = ..count../1000), stat = "bin", bins = 50, alpha = 0.5, color = "black") +
  scale_y_continuous(breaks = c(0,3)) +
  coord_flip()

ggsave("../figures/mageck_polysome_histogram.pdf", width = 1, height = 2.9)

# Volcano plot of screen results
plot_data <- mageck_gene %>%
  filter(control == "mono", rep == "rep1") %>% 
  mutate(annotation = case_when(
    str_detect(gene, "^RPL") ~ "1 60S ribosome",
    str_detect(gene, "^RPS|RACK1") ~ "2 40S ribosome",
    gene %in% (ribi %>% filter(annotation == "60S biogenesis") %>% pull(gene)) ~ "3 60S biogenesis",
    gene %in% (ribi %>% filter(annotation == "40S biogenesis") %>% pull(gene)) ~ "4 40S biogenesis",
    TRUE ~ NA_character_
  )) %>%
  mutate(annotation = if_else(fdr > fdr_cutoff, NA_character_, annotation)) %>%
  mutate(annotation = if_else(is.na(annotation), "7 Other", as.character(annotation))) %>%
  write_csv("../../../../source_data/figure_2c.csv")

plot_data %>%
  ggplot(aes(x = lfc, y = -log10(p_value))) +
  facet_wrap(~treatment, ncol = 3, scales = "free_x") +
  geom_vline(aes(xintercept = 0), linetype = 2, linewidth = 0.5, color = "grey") +
  geom_point(aes(color = annotation, fill = annotation),
    size = 0.5, alpha = 0.5, shape = 'circle',
    data = plot_data  %>% filter(annotation == "7 Other"),
  ) +
  geom_point(aes(color = annotation, fill = annotation),
    size = 1, alpha = 1,
    data = plot_data  %>% filter(annotation != "7 Other"),
  ) +
  scale_color_manual(values = cbPalette_12[c(10,5,2,6,12)]) +
  scale_fill_manual(values = cbPalette_12[c(10,2,5,6,12)]) +
  scale_shape_manual(values = rep(c("cross", "plus", "square open", "diamond open", "triangle down filled", "triangle filled", "circle"), 5)) +
  labs(x = "log2 polysome / monosome", y = "-log10 P-value") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = c(0,2,4,6), limits = c(0,6)) +
  theme(legend.title = element_blank())

ggsave("../figures/polysome_volcano.pdf", width = 5.2, height = 2.2, units = "in")

# Volcano plot of supernatant comparisons
plot_data <- mageck_gene %>%
  filter(control == "supernatant", rep == "rep1") %>% 
  mutate(annotation = case_when(
    str_detect(gene, "^RPL") ~ "1 60S ribosome",
    str_detect(gene, "^RPS|RACK1") ~ "2 40S ribosome",
    gene %in% (ribi %>% filter(annotation == "60S biogenesis") %>% pull(gene)) ~ "3 60S biogenesis",
    gene %in% (ribi %>% filter(annotation == "40S biogenesis") %>% pull(gene)) ~ "4 40S biogenesis",
    TRUE ~ NA_character_
  )) %>%
  mutate(annotation = if_else(fdr > fdr_cutoff, NA_character_, annotation)) %>%
  mutate(annotation = if_else(is.na(annotation), "7 Other", as.character(annotation))) %>%
  filter(treatment == "mono") %>%
  write_csv("../../../../source_data/figure_2e_volcano.csv")

plot_data %>%
  ggplot(aes(x = lfc, y = -log10(p_value))) +
  geom_vline(aes(xintercept = 0), linetype = 2, linewidth = 0.5, color = "grey") +
  geom_point(aes(color = annotation, fill = annotation),
    size = 0.5, alpha = 0.5, shape = 'circle',
    data = plot_data  %>% filter(annotation == "7 Other"),
  ) +
  geom_point(aes(color = annotation, fill = annotation),
    size = 1, alpha = 1,
    data = plot_data  %>% filter(annotation != "7 Other"),
  ) +
  scale_color_manual(values = cbPalette_12[c(10,5,2,6,12)]) +
  scale_fill_manual(values = cbPalette_12[c(10,2,5,6,12)]) +
  scale_shape_manual(values = rep(c("cross", "plus", "square open", "diamond open", "triangle down filled", "triangle filled", "circle"), 5)) +
  labs(x = "log2 (monosome/supernatant)", y = "-log10 P-value") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = c(0,2,4,6), limits = c(0,6)) +
  theme(legend.title = element_blank())

ggsave("../figures/polysome_volcano_supernatant.pdf", width = 4.5, height = 2.2, units = "in")

# Plot ribosomal effects for monosome comparisons
plot_data <- mageck_gene %>%
  filter(control == "mono", rep == "rep1") %>% 
  mutate(annotation = case_when(
    str_detect(gene, "^RPL") ~ "4 60S ribosome",
    str_detect(gene, "^RPS|RACK1") ~ "3 40S ribosome",
    gene %in% (ribi %>% filter(annotation == "60S biogenesis") %>% pull(gene)) ~ "2 60S biogenesis",
    gene %in% (ribi %>% filter(annotation == "40S biogenesis") %>% pull(gene)) ~ "1 40S biogenesis",
    TRUE ~ NA_character_
  )) %>%
  mutate(annotation = if_else(is.na(annotation), "7 Other", as.character(annotation))) %>%
  mutate(sig = if_else(fdr < fdr_cutoff, "sig", "notsig"))

mylabel <- function (x) {
  str_replace(x, ". ", "")
}

plot_data %>%
  filter(!str_detect(annotation, "Other")) %>% 
  filter(treatment == "heavy") %>%
  write_csv("../../../../source_data/figure_2d.csv") %>%
  ggplot(aes(x = annotation, y = lfc, fill = annotation, color = annotation)) +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.5, color = "grey") +
  stat_summary(fun = median, width = 0.5, linewidth = 1, color = "grey", geom = "tile", fill = "grey") +
  geom_point(aes(shape = sig), size = 1, position = position_jitter(height = 0, width = 0.3, seed = 11)) +
  scale_x_discrete(labels = mylabel) +
  scale_color_manual(values = cbPalette_12[c(6,2,5,10,12)]) +
  scale_fill_manual(values = cbPalette_12[c(6,2,5,10,12)]) +
  scale_size_manual(values = c(1, 2)) +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_shape_manual(values = c("circle open", "circle")) +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  labs(y = "log2 (heavy \n polysome / monosome)") +
  coord_flip()

ggsave("../figures/polysome_ribosome_groups.pdf", width = 2, height = 2, units = "in")

# Plot ribosomal effects for supernatant comparisons
plot_data <- mageck_gene %>%
  filter(control == "supernatant", rep == "rep1") %>% 
  mutate(annotation = case_when(
    str_detect(gene, "^RPL") ~ "4 60S ribosome",
    str_detect(gene, "^RPS|RACK1") ~ "3 40S ribosome",
    gene %in% (ribi %>% filter(annotation == "60S biogenesis") %>% pull(gene)) ~ "2 60S biogenesis",
    gene %in% (ribi %>% filter(annotation == "40S biogenesis") %>% pull(gene)) ~ "1 40S biogenesis",
    TRUE ~ NA_character_
  )) %>%
  mutate(annotation = if_else(is.na(annotation), "7 Other", as.character(annotation))) %>%
  mutate(sig = if_else(fdr < fdr_cutoff, "sig", "notsig"))

plot_data %>%
  filter(!str_detect(annotation, "Other")) %>%
  write_csv("../../../../source_data/figure_2e_s2c_ribosome_groups.csv") %>%
  ggplot(aes(x = annotation, y = lfc, fill = annotation, color = annotation)) +
  facet_wrap(~treatment, ncol = 3, scales = "free") +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.5, color = "grey") +
  stat_summary(fun = median, width = 0.5, linewidth = 1, color = "grey", geom = "tile", fill = "grey") +
  geom_point(aes(shape = sig), size = 1, position = position_jitter(height = 0, width = 0.3, seed = 11)) +
  scale_x_discrete(labels = mylabel) +
  scale_color_manual(values = cbPalette_12[c(6,2,5,10,12)]) +
  scale_fill_manual(values = cbPalette_12[c(6,2,5,10,12)]) +
  scale_size_manual(values = c(1, 2)) +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_shape_manual(values = c("circle open", "circle")) +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  labs(y = "log2 ribosome fraction / supernatant") +
  coord_flip()

ggsave("../figures/supernatant_ribosome_groups.pdf", width = 6, height = 2, units = "in")

# Plot other translation groups
plot_data <- mageck_gene %>%
  filter(control == "mono", rep == "rep1", treatment == "heavy") %>% 
  mutate(annotation = case_when(
    str_detect(gene, "^CNOT") ~ "4 CCR4-NOT deadenylase",
    str_detect(gene, "^PSM") ~ "2 Proteasome",
    str_detect(gene, "^CCT|TCP1") ~ "3 TRiC chaperonin",
    str_detect(gene, "^POLR2") ~ "1 RNA polymerase II",
    str_detect(gene, "^SF3A|^SF3B|^PHF5") ~ "0 Splicing Factor 3a/b",
    str_detect(gene, "^EIF2S") ~ "9 EIF2",
    str_detect(gene, "^EIF2B") ~ "8 EIF2B",
    str_detect(gene, "^EIF3") ~ "7 EIF3",
    str_detect(gene, "EIF4A1$|EIF4G1|EIF4E$") ~ "6 EIF4F",
    str_detect(gene, ".ARS.") ~ "5 tRNA synthetases",
    TRUE ~ NA_character_
  )) %>%
  mutate(annotation = if_else(is.na(annotation), "7 Other", as.character(annotation))) %>%
  mutate(sig = if_else(fdr < fdr_cutoff, "sig", "notsig"))

plot_data %>%
  filter(!str_detect(annotation, "Other")) %>%
  write_csv("../../../../source_data/figure_2f.csv") %>%
  ggplot(aes(x = annotation, y = lfc, fill = annotation, color = annotation, shape = sig)) +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.5, color = "grey") +
  geom_point(position = position_jitter(height = 0, width = 0.2, seed = 11)) +
  scale_color_manual(values = cbPalette_12[c(2,3,4,5,6,9,8,7,10,11)]) +
  scale_x_discrete(labels = mylabel) +
  scale_fill_manual(values = cbPalette_12[c(2,3,4,5,6,9,8,7,10,11)]) +
  scale_size_manual(values = c(1, 2)) +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_shape_manual(values = c("circle open", "circle")) +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  labs(y = "log2 (heavy \n polysome/monosome)") +
  coord_flip()

ggsave("../figures/polysome_translation_groups_2.pdf", width = 3, height = 3.5, units = "in")

# Load MaGeCK fitness results for mRNA
mageck_fitness <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc.+day(7|13).+total.+total")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names() %>%
  select(-neg_lfc) %>%
  rename(lfc = pos_lfc, gene = id) %>%
  mutate(p_value = if_else(lfc > 0, pos_p_value, neg_p_value)) %>%
  mutate(fdr = if_else(lfc > 0, pos_fdr, neg_fdr)) %>%
  mutate(goodsgrna = if_else(lfc > 0, pos_goodsgrna, neg_goodsgrna)) %>%
  mutate(day = str_extract(treatment, "day\\d+"))

# Load MaGeCK genomic DNA fitness results
mageck_grna_fitness <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc.+day(5|13).+grna.+grna")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names() %>%
  select(-neg_lfc) %>%
  rename(lfc = pos_lfc, gene = id) %>%
  mutate(p_value = if_else(lfc > 0, pos_p_value, neg_p_value)) %>%
  mutate(fdr = if_else(lfc > 0, pos_fdr, neg_fdr)) %>%
  mutate(goodsgrna = if_else(lfc > 0, pos_goodsgrna, neg_goodsgrna)) %>%
  mutate(day = str_extract(treatment, "day\\d+"))

# Plot polysome vs mRNA fitness
plot_data <- mageck_gene %>% 
  filter(control == "mono", treatment == "heavy", rep == "rep1") %>% 
  select(treatment, control, rep, gene, lfc, fdr, p_value) %>% 
  left_join(mageck_fitness  %>% select(day, gene, lfc, fdr), by = "gene", suffix = c ("", ".fitness")) %>%
  write_csv("../../../../source_data/figure_s2f_mrna.csv")

plot_data %>% 
  ggplot(aes(x = lfc, y = lfc.fitness)) +
  facet_wrap(~fct_rev(day), ncol = 2, scales = "fixed") + 
  geom_point(size = 0.5, alpha = 0.2) +
  labs(x = "log2 heavy polysome / monosome", y = "log2 mRNA depletion\nw.r.t. day 1") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

ggsave("../figures/polysome_vs_mrna_fitness_all.pdf", width = 4.2, height = 2.5)

# Plot polysome vs gDNA fitness
plot_data <- mageck_gene %>% 
  filter(control == "mono", treatment == "heavy", rep == "rep1") %>% 
  select(treatment, control, rep, gene, lfc, fdr, p_value) %>% 
  left_join(mageck_grna_fitness  %>% select(day, gene, lfc, fdr), by = "gene", suffix = c ("", ".fitness")) %>%
  write_csv("../../../../source_data/figure_s2f_gdna.csv")

plot_data %>% 
  ggplot(aes(x = lfc, y = lfc.fitness)) +
  facet_wrap(~fct_rev(day), ncol = 2, scales = "fixed") + 
  geom_point(size = 0.5, alpha = 0.2) +
  labs(x = "log2 heavy polysome / monosome", y = "log2 genomic DNA depletion\nw.r.t. day 1") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(limits = c(NA, 2.5))

ggsave("../figures/polysome_vs_grna_fitness_all.pdf", width = 4.2, height = 2.5)

# Plot mRNA fitness vs polysome monosome ratio for specific gene groups
plot_data <- mageck_gene %>% 
  filter(control == "mono", treatment == "heavy", rep == "rep1") %>% 
  select(treatment, control, rep, gene, lfc, fdr, p_value) %>% 
  left_join(mageck_fitness  %>% select(day, gene, lfc, fdr), by = "gene", suffix = c ("", ".fitness")) %>% 
  mutate(annotation = case_when(
    str_detect(gene, "^RPL") ~ "1 60S ribosome \n& biogenesis",
    str_detect(gene, "^RPS") ~ "2 40S ribosome \n& biogenesis",
    gene %in% (ribi %>% filter(annotation == "60S biogenesis") %>% pull(gene)) ~ "1 60S ribosome \n& biogenesis",
    gene %in% (ribi %>% filter(annotation == "40S biogenesis") %>% pull(gene)) ~ "2 40S ribosome \n& biogenesis",
    str_detect(gene, "^EIF3") ~ "5 EIF3",
    str_detect(gene, "^PSM") ~ "6 Proteasome",
    str_detect(gene, "^POLR2") ~ "7 RNA polymerase II",
    TRUE ~ NA_character_
  )) %>%
  mutate(annotation = if_else(is.na(annotation), "0 Other", as.character(annotation))) %>%
  filter(annotation != "0 Other") %>%
  write_csv("../../../../source_data/figure_2g_s2g_mrna.csv")

fit_data <- plot_data  %>% 
  group_by(annotation, day) %>% 
  nest() %>% 
  mutate(fit = map(data, ~lm(lfc ~ lfc.fitness, data = .x))) %>%
  mutate(fit = map(fit, .  %>% predict(interval = "confidence") %>% as_tibble())) %>% 
  unnest(data, fit)

plot_data %>% 
  ggplot(aes(x = lfc, y = lfc.fitness, color = annotation, fill = annotation)) +
  facet_wrap(~fct_rev(day), ncol = 2, scales = "free") + 
  scale_color_manual(values = cbPalette_12[c(2,6,7,4,3)]) +
  scale_fill_manual(values = cbPalette_12[c(2,6,7,4,3)]) +
  scale_shape_manual(values = c("circle", "cross", "plus", "square", "diamond filled", "triangle filled", "asterisk", "plus")) +
  geom_point(size = 0.5, alpha = 0.2, data = plot_data  %>%  filter(annotation == "0 Other")) +
  geom_ribbon(aes(xmin = lwr, xmax = upr, y = lfc.fitness), data = fit_data, alpha = 0.3, linewidth = 0, show.legend = F) +
  geom_point(size = 1, data = plot_data  %>%  filter(annotation != "0 Other")) +
  labs(x = "log2 heavy polysome / monosome", y = "log2 mRNA depletion\nw.r.t. day 1") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))

ggsave("../figures/polysome_vs_mrna_fitness.pdf", width = 5.2, height = 4.5)

# Plot polysome ReLiC against gRNA fitness
plot_data <- mageck_gene %>% 
  filter(control == "mono", treatment == "heavy", rep == "rep1") %>% 
  select(treatment, control, rep, gene, lfc, fdr, p_value) %>% 
  left_join(mageck_grna_fitness  %>% select(day, gene, lfc, fdr), by = "gene", suffix = c ("", ".fitness")) %>% 
  mutate(annotation = case_when(
    str_detect(gene, "^RPL") ~ "1 60S ribosome \n& biogenesis",
    str_detect(gene, "^RPS") ~ "2 40S ribosome \n& biogenesis",
    gene %in% (ribi %>% filter(annotation == "60S biogenesis") %>% pull(gene)) ~ "1 60S ribosome \n& biogenesis",
    gene %in% (ribi %>% filter(annotation == "40S biogenesis") %>% pull(gene)) ~ "2 40S ribosome \n& biogenesis",
    str_detect(gene, "^EIF3") ~ "5 EIF3",
    str_detect(gene, "^PSM") ~ "6 Proteasome",
    str_detect(gene, "^POLR2") ~ "7 RNA polymerase II",
    TRUE ~ NA_character_
  )) %>%
  mutate(annotation = if_else(is.na(annotation), "0 Other", as.character(annotation))) %>%
  filter(annotation != "0 Other") %>%
  write_csv("../../../../source_data/figure_s2g_gdna.csv")

fit_data <- plot_data  %>% 
  group_by(annotation, day) %>% 
  nest() %>% 
  mutate(fit = map(data, ~lm(lfc ~ lfc.fitness, data = .x))) %>%
  mutate(fit = map(fit, .  %>% predict(interval = "confidence") %>% as_tibble())) %>% 
  unnest(data, fit)

plot_data %>% 
  ggplot(aes(x = lfc, y = lfc.fitness, color = annotation, fill = annotation)) +
  facet_wrap(~fct_rev(day), ncol = 2, scales = "free") + 
  scale_color_manual(values = cbPalette_12[c(2,6,7,4,3)]) +
  scale_fill_manual(values = cbPalette_12[c(2,6,7,4,3)]) +
  scale_shape_manual(values = c("circle", "cross", "plus", "square", "diamond filled", "triangle filled", "asterisk", "plus")) +
  geom_point(size = 0.5, alpha = 0.2, data = plot_data  %>%  filter(annotation == "0 Other")) +
  geom_ribbon(aes(xmin = lwr, xmax = upr, y = lfc.fitness), data = fit_data, alpha = 0.3, linewidth = 0, show.legend = F) +
  geom_point(size = 1, data = plot_data  %>%  filter(annotation != "0 Other")) +
  labs(x = "log2 heavy polysome / monosome", y = "log2 genomic DNA depletion\nw.r.t. day 1") +
  theme(legend.position = "bottom", legend.direction = "vertical") +
  guides(fill = guide_legend(ncol = 2), color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))

ggsave("../figures/polysome_vs_grna_fitness.pdf", width = 5.2, height = 4.5)

# Plot sgRNA effects for elongation hits
plot_data <- mageck_sgrna %>% 
  filter(rep == "rep1", treatment == "heavy", control == "mono") %>%
  filter(gene %in% c("ASCC3", "ZNF598", "EIF5A", "METAP2", "EEF2", "EIF4A1", "EIF2S1", "EIF3A", "EIF5B", "EEF1A2", "EEF1A1")) %>%
  mutate(annotation = case_when(
    gene %in% c("EEF2", "EIF5A", "ASCC3", "METAP2", "ZNF598", "EEF1A1", "EEF1A2") ~ "Translation elongation",
    gene %in% c("EIF4A1", "EIF2S1", "EIF3A", "EIF5B") ~ "Translation initiation"
  )) %>%
  mutate(gene = fct_reorder(gene, -lfc)) %>%
  write_csv("../../../../source_data/figure_2h.csv")

plot_data  %>% 
  ggplot(aes(x = gene, y = lfc, color = annotation, fill = annotation, shape = annotation)) +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.5, color = "grey") +
  stat_summary(fun = median, width = 0.5, color = "grey", geom = "tile", linewidth = 1) +
  geom_point(position = position_jitter(height = 0, width = 0.1, seed = 11), size = 1, alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.33)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_color_manual(values = cbPalette_12[c(4,7)]) +
  scale_fill_manual(values = cbPalette_12[c(4,7)]) +
  scale_shape_manual(values = rep(c("triangle filled", "triangle down filled"), 5)) +
  labs(y = "log2 (heavy \n polysome/monosome", x = "sgRNA") +
  theme(legend.position = "none")

ggsave("../figures/elongation_vs_initiation_sgrna.pdf", width = 2.5, height = 1.8)

# Read GORILLA GO terms
gorilla <- read_csv("../data/gorilla/ntc_polysome_bc1_heavy_rep1_vs_ntc_polysome_bc1_mono_rep1.csv", show_col_types = F) %>% 
  janitor::clean_names() %>% 
  mutate(enrichment = str_extract(enrichment_n_b_n_b, "\\d+\\.\\d+")) %>%
  select(-genes, -enrichment_n_b_n_b, -p_value) %>% 
  write_csv("../../../../source_data/figure_s2b.csv")

# Show selected GO terms enriched in downregulated genes
subset_go_terms <- c(
  "GO:0006364" = "rRNA processing",
  "GO:0006413" = "translational initiation",
  "GO:0022625" = "cytosolic large ribosomal subunit",
  "GO:0022627" = "cytosolic small ribosomal subunit",
  "GO:0000502" = "proteasome complex",
  "GO:0005832" = "chaperonin-containing T-complex"
)

selected_go_down <- gorilla  %>% 
  filter(go_term %in% names(subset_go_terms)) %>%
  select(go_term, description, fdr_q_value, enrichment)

print("Selected GO terms (downregulated):")
print(selected_go_down)

# Show selected GO terms enriched in upregulated genes
subset_go_terms <- c(
  "GO:0008380" = "RNA splicing",
  "GO:0031124" = "mRNA 3'-end processing",
  "GO:0006406" = "mRNA export from nucleus"
)

selected_go_up <- gorilla  %>% 
  filter(go_term %in% names(subset_go_terms)) %>%
  select(go_term, description, fdr_q_value, enrichment)

print("Selected GO terms (upregulated):")
print(selected_go_up)

print("Analysis complete!")