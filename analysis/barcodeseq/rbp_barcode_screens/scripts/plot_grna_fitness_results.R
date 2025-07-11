# Plot results of barcode depletion screens

# Load libraries
options(warn = -1)
suppressPackageStartupMessages({
  library(tidyverse)
  library(rasilabRtemplates)
})

cbPalette_12 <- c(
  "#88CCEE", "#CC6677", "#117733", "#999933", "#332288", "#AA4499",
  "#661100", "#44AA99", "#882255", "#6699CC", "#DDCC77", "#888888"
)

theme_set(theme_rasilab() + 
  theme(
    axis.line = element_line(color = "grey"), 
    axis.title.y = element_text(margin = margin(r=10)),
    axis.title.x = element_text(margin = margin(t=10))
  )
)

umi_cutoff <- 20
umi_grna_cutoff <- 1

# Read sample annotations
sample_annotations <- read_csv("../annotations/sample_annotations.csv", show_col_types = F) %>%
  print()

# Read barcode counts for day 1 gDNA and mRNA
barcode_counts <- c("../data/subpool_barcode_counts/219p138.csv","../data/subpool_barcode_counts/219p27.csv") %>%
  as_tibble_col("file") %>%
  mutate(sample_id = str_extract(file, "219p\\d+")) %>%
  inner_join(sample_annotations %>% select(sample_name, sample_id), by = "sample_id") %>%
  mutate(data = map(file, ~read_csv(.x, show_col_types = F))) %>%
  select(-file) %>%
  unnest() %>%
  print()

# CDF of barcode counts on day 1 for gDNA and cDNA
plot_data <- barcode_counts %>%
  mutate(cutoff = if_else(str_detect(sample_name, "grna"), read_count >= umi_grna_cutoff, umi_count >= umi_cutoff)) %>%
  filter(cutoff) %>%
  group_by(sample_name, insert_num) %>%
  summarize(n_barcodes = n(), .groups = "drop") %>%
  group_by(sample_name) %>%
  arrange(desc(n_barcodes)) %>%
  mutate(x = 1:n()) %>%
  ungroup() %>%
  filter(n_barcodes < 100) %>%
  mutate(type = if_else(str_detect(sample_name, "total"), "mRNA", "gDNA"))

median_data <- plot_data %>%
  group_by(sample_name) %>%
  summarize(n_barcodes = median(n_barcodes), .groups = "drop") %>%
  write_csv("../../../../source_data/figure_s1e.csv")

plot_data %>%
  ggplot(aes(x = x, y = n_barcodes, color = type, group = type)) +
  geom_hline(aes(yintercept = n_barcodes), linetype = "dashed", data = median_data) + 
  geom_line(linewidth = 1) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0, 8800)) +
  scale_color_manual(values = cbPalette[c(3,7)]) +
  theme_rasilab() +
  labs(x = "sgRNA ID", y = "Number of barcodes", color = "")

ggsave("../figures/cdf_sgrna_n_barcodes.pdf", width = 4, height = 2)

# PDF of read counts for sgRNA
plot_data <- barcode_counts %>%
  mutate(umi_count = if_else(str_detect(sample_name, "grna"), read_count, umi_count)) %>%
  group_by(sample_name, insert_num) %>%
  summarize(umi_count = sum(umi_count), .groups = "drop") %>%
  mutate(type = if_else(str_detect(sample_name, "total"), "mRNA", "gDNA")) %>%
  write_csv("../../../../source_data/figure_s1d.csv")

plot_data %>% 
  ggplot(aes(x = log10(umi_count), color = type)) +
  geom_freqpoly(binwidth=0.1, linewidth = 1) +
  scale_color_manual(values = rasilabRtemplates::cbPalette[c(3,7)]) +
  theme_rasilab() +
  labs(x = "log10 read counts", y = "Number of sgRNAs", color = "") +
  scale_x_continuous(limits = c(0, 5))

ggsave("../figures/pdf_sgrna_umi_counts.pdf", width = 3.6, height = 2)

# Load MaGeCK gene hit data
mageck_gene <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc.+bc1.+(total|grna).+ntc.+bc1.+(total|grna)"), !str_detect(file, "total.+grna")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names() %>%
  select(-pos_lfc) %>%
  rename(lfc = neg_lfc, gene = id) %>%
  mutate(lfc = round(lfc, 2)) %>%
  arrange(neg_rank)

# Read essential gene list
essential_genes <- read_csv("../annotations/DepMapCRISPRInferredCommonEssentials_23Q2.csv", show_col_types = F) %>% 
  separate(Essentials, c("gene", "num")) %>% 
  select(gene) %>%
  mutate(essential = TRUE) %>%
  print()

# Plot mRNA fold-change vs. gDNA fold-change after Cas9 induction
plot_data = mageck_gene %>% 
  mutate(day = str_extract(treatment, "day\\d+"), type = str_extract(treatment, "total|grna")) %>%
  filter(day %in% c("day5", "day13", 'day21'), str_detect(sample_name, "total.+total|grna.+grna")) %>%
  mutate(type = if_else(type == "total", "mRNA", "gDNA")) %>%
  select(gene, day, type, lfc) %>% 
  pivot_wider(names_from = type, values_from = lfc) %>%
  mutate(daynum = as.integer(str_replace(day, "day", ""))) %>%
  mutate(day = fct_reorder(paste0("day ", daynum), daynum)) %>%
  select(-daynum) %>%
  arrange(gene, day) %>%
  write_csv("../../../../source_data/figure_1d.csv")

correlation <- plot_data %>% 
  group_by(day) %>% 
  nest() %>%
  mutate(cor = map(data, ~cor.test(~ mRNA + gDNA, data=.x, method = "spearman", exact = F))) %>%
  mutate(cor = map(cor, broom::tidy)) %>%
  select(-data) %>%
  unnest(cor)

plot_data %>%
  ggplot(aes(x = gDNA, y = mRNA)) +
  facet_wrap(~day) + 
  geom_point(size = 0.5, alpha = 0.2) +
  geom_text(data = correlation, aes(label = paste0("r = ", round(estimate, 2))), x = -1, y = -6, hjust = 0, vjust = 0, size = 4) +
  scale_x_continuous(breaks = c(-6, -3, 0, 3), limits = c(-5, 3)) +
  scale_y_continuous(breaks = c(-6, -3, 0, 3), limits = c(-6, 3)) +
  labs(
    x = "gDNA fold-change (log2, median-centered)",
    y = "mRNA fold-change\n(log2, median-centered)")

ggsave("../figures/compare_mrna_gdna_foldchange_with_time.pdf", width = 4, height = 2)

# NTC mRNA fold change with time
plot_data = mageck_gene %>% 
  filter(str_detect(treatment, "total"), control == "ntc_day1_bc1_total") %>%
  mutate(day = str_extract(treatment, "(?<=day)\\d+")) %>%
  filter(day %in% c("5", "13", "21")) %>%
  select(gene, day, lfc) %>%
  left_join(essential_genes, by = "gene") %>%
  group_by(day, essential) %>% 
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(essential = if_else(is.na(essential), "non-essential", "essential")) %>%
  mutate(day = fct_reorder(day, -as.integer(day))) %>%
  arrange(gene, day, essential) %>%
  write_csv("../../../../source_data/figure_1e_mrna.csv")

plot_data %>%
  ggplot(aes(x = lfc, y = day, fill = fct_rev(essential))) +
  ggridges::geom_density_ridges(alpha = 0.8) +
  labs(x = "mRNA fold-change\n(log2, arb.units)", y = "days post Cas9", fill = "") +
  theme(axis.line = element_line(color = "grey"),
        axis.ticks = element_line(color = "grey"), 
        legend.position = "top", legend.direction = "vertical", legend.key.height = unit(0.2, "cm")) +
  scale_y_discrete(expand = c(0.1, 0.1)) +
  scale_x_continuous(limits = c(-8,4), breaks = c(-8, -4, 0, 4))

ggsave("../figures/ntc_total_mrna_foldchange_with_time.pdf", width = 2.4, height = 3.2, units = "in")

# NTC genomic DNA fold change with time
plot_data = mageck_gene %>% 
  filter(str_detect(treatment, "grna"), control == "ntc_day1_bc1_grna") %>%
  mutate(day = str_extract(treatment, "(?<=day)\\d+")) %>%
  filter(day %in% c("5", "13", "21")) %>%
  select(gene, day, lfc) %>%
  left_join(essential_genes, by = "gene") %>%
  group_by(day, essential) %>% 
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(essential = if_else(is.na(essential), "non-essential", "essential")) %>%
  mutate(day = fct_reorder(day, -as.integer(day))) %>%
  arrange(gene, day, essential) %>%
  write_csv("../../../../source_data/figure_1e_gdna.csv")

plot_data %>%
  ggplot(aes(x = lfc, y = day, fill = fct_rev(essential))) +
  ggridges::geom_density_ridges(alpha = 0.8) +
  labs(x = "gDNA fold-change\n(log2, arb.units)", y = "days post Cas9", fill = "") +
  theme(axis.line = element_line(color = "grey"),
        axis.ticks = element_line(color = "grey"), 
        legend.position = "top", legend.direction = "vertical", legend.key.height = unit(0.2, "cm")) +
  scale_y_discrete(expand = c(0.1, 0.1)) +
  scale_x_continuous(limits = c(-6,4), breaks = c(-6, -3, 0, 3))

ggsave("../figures/gdna_foldchange_with_time.pdf", width = 2.4, height = 3.2)

# Plot gRNA depletion of SF3b RNA-seq candidates
subset_data <- mageck_gene %>%
  filter(str_detect(treatment, "grna")) %>%
  mutate(gene = if_else(gene == "PHF5A", "SF3B7", gene)) %>%
  filter(gene %in% c("NLUC", "FLUC", "EGFP", "SF3B5", "SF3B6") | str_detect(gene, "SF3B")) %>%
  mutate(day = as.numeric(str_extract(treatment, "(?<=day)\\d+"))) %>%
  select(day, gene, neg_p_value, lfc)
  
plot_data <- bind_rows(subset_data, subset_data %>% distinct(gene) %>% mutate(day = 1, neg_p_value = 1, lfc = 0)) %>% 
  mutate(sig = case_when(
            neg_p_value >= 0.05 ~ 1,
            neg_p_value >= 0.01 & neg_p_value < 0.05 ~ 2,
            neg_p_value >= 0.001 & neg_p_value < 0.01 ~ 3,
            neg_p_value < 0.001 ~ 4
          )) %>%
  group_by(gene) %>%
  mutate(ordering = lfc[day == 21]) %>%
  ungroup() %>%
  mutate(gene = fct_reorder(gene, -ordering)) %>%
  arrange(gene, day) %>%
  write_csv("../../../../source_data/figure_s3d.csv")

plot_data %>%
  ggplot(aes(x = day, y = lfc, color = gene, group = gene, size = as.factor(sig))) +
  geom_point() +
  scale_color_manual(values = cbPalette_12[c(1,8,3,4,5,6,9,7,2,10,11)]) +
  scale_size_manual(values = c(0.5, 1, 1.5, 2)) +
  geom_line(size = 1, alpha = 0.5) +
  labs(x = "Days post Cas9", y = "log2 genomic DNA\ndepletion w.r.t day 1", color = "sgRNA") +
  theme(legend.box = "horizontal") +
  guides(color = guide_legend(ncol = 2), size = "none")

ggsave("../figures/sf3b_fitness.pdf", width = 4.5, height = 2.2, units = "in")

# Plot gRNA depletion of eIF2,3,4 complex subunits
set.seed(1010)

subset_data <- mageck_gene %>%
  filter(str_detect(treatment, "grna")) %>%
  mutate(gene_group = case_when(
    str_detect(gene, "FLUC|NLUC") ~ 'CONTROL',
    str_detect(gene, "EIF4(E$|A1|G1)") ~ 'eIF4F',
    str_detect(gene, "EIF2S") ~ 'eIF2',
    str_detect(gene, "EIF2B") ~ 'eIF2B',
    str_detect(gene, "EIF3") ~ 'eIF3'
  )) %>%
  filter(!is.na(gene_group)) %>%
  mutate(day = as.numeric(str_extract(treatment, "(?<=day)\\d+"))) %>%
  select(day, gene, neg_p_value, lfc, gene_group)

plot_data <- subset_data %>% 
  bind_rows(subset_data %>% distinct(gene, gene_group) %>% mutate(day = 1, neg_p_value = 1, lfc = 0)) %>%
  group_by(gene) %>%
  mutate(ordering = lfc[day == 21]) %>%
  ungroup() %>%
  mutate(gene = fct_reorder(gene, -ordering)) %>%
  arrange(gene, day) %>%
  write_csv("../../../../source_data/figure_s4e.csv")

control_data <- plot_data %>% 
  filter(gene_group == "CONTROL") %>%
  select(day, gene, lfc)

plot_data %>%
  filter(gene_group != "CONTROL") %>%
  ggplot(aes(x = day, y = lfc, group = gene, color = gene)) +
  facet_wrap(~gene_group, scales = "free", ncol = 2) +
  geom_point() +
  geom_text(aes(label = gene), data = . %>% filter(day == 21), size = 3, hjust = -0.1, alpha = 0.5, position = position_jitter(width = 0, height = 0.5)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  scale_size_manual(values = c(1, 2, 3, 4)) +
  scale_color_viridis_d(end = 0.8) +
  geom_line(size = 1, alpha = 0.5) +
  geom_line(data = control_data, size = 1, alpha = 0.5, color = "#AAAAAA") +
  geom_point(data = control_data, color = "#888888") +
  geom_text(aes(label = gene), data = control_data %>% filter(day == 21), size = 3, hjust = -0.1, color = "#888888") +
  labs(x = "Days post Cas9", y = "log2 genomic DNA\ndepletion (w.r.t day 1)") +
  theme(legend.position = "none", legend.box = "horizontal", panel.spacing = unit(0.5, "in")) +
  coord_cartesian(clip = "off", xlim = c(0, 21), ylim = c(-4.5, 2.5))

ggsave("../figures/eif_fitness_mrna.pdf", width = 5, height = 5, units = "in")