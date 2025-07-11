#!/usr/bin/env Rscript

# Analyze NMD screen
# Converted from plot_nmd_results.ipynb

# Load libraries and define analysis parameters
options(warn = -1)

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(tidyverse)
  library(rasilabRtemplates)
})

cbPalette_12 <- c(
  "#DDCC77", "#CC6677", "#6699CC", "#661100", "#117733", "#999933",
  "#332288", "#AA4499", "#44AA99", "#882255", "#88CCEE", "#999999"
)

theme_set(theme_rasilab() +
  theme(
    axis.line = element_line(color = "grey"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  ))

fdr_cutoff <- 0.05
p_value_cutoff <- 0.05

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

# Load gene annotations for NMD screen
nmd_hit_annotations <- read_csv("../annotations/nmd_hit_annotations.csv", show_col_types = FALSE)

annotation_order <- c(
"Core NMD",
"Translation initiation",
"tRNA metabolism",
"Mitochondria",
"RNA exosome",
"Proteasome",
"60S ribosome",
"40S ribosome",
"Ribosome biogenesis"
) %>%
as_tibble_col("annotation") %>%
mutate(priority = row_number())

# Load MaGeCK gene-level results for NMD screens
mageck_gene <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc|ptc"), !str_detect(file, "retained|skipped"), !str_detect(file, "ptc_bc1_total_vs"),
  !str_detect(file, "ntc.+day"), !str_detect(file, ".+polysome.+vs.+(total|light)")) %>%
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
  mutate(fdr = if_else(is.na(fdr), 1, fdr))

# Load MaGeCK sgRNA-level results for NMD screens
mageck_sgrna <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "sgrna_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc|ptc"), !str_detect(file, "retained|skipped"), !str_detect(file, "ptc_bc1_total_vs"),
  !str_detect(file, "ntc.+day"), !str_detect(file, ".+polysome.+vs.+(total|light)")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names()

# Plot alphabetical list of NMD hits
plot_data <- mageck_gene %>%
  filter(str_detect(sample_name, "bc1.+bc2")) %>%
  mutate(reporter = str_extract(sample_name, "ptc|ntc")) %>%
  left_join(nmd_hit_annotations, by = "gene") %>%
  left_join(annotation_order, by = "annotation") %>%
  mutate(priority = if_else(is.na(annotation), as.integer(10), priority)) %>%
  mutate(annotation = if_else(is.na(annotation), "Other", as.character(annotation))) %>%
  mutate(annotation = fct_reorder(annotation, priority)) %>%
  group_by(reporter) %>%
  arrange(gene) %>%
  mutate(x = row_number()) %>%
  ungroup() %>%
  mutate(reporter = toupper(reporter)) %>%
  write_csv("../../../../source_data/figure_4b_s4c.csv")

plot_data %>%
  ggplot(aes(x = x, y = -log10(pos_p_value), color = annotation, fill = annotation)) +
  facet_wrap(~ fct_rev(reporter), ncol = 1) +
  geom_point(size = 1, alpha = 1, data = plot_data %>% filter(fdr > fdr_cutoff | lfc < 0), color = "grey", fill = "grey") +
  geom_point(aes(shape = annotation),
    size = 2, alpha = 1,
    data = plot_data %>% filter(fdr < fdr_cutoff & lfc > 0),
    position = position_jitter(height = 0, width = 50, seed = 11)
  ) +
  scale_color_manual(values = cbPalette_12[c(2,5,4, 3,7, 7, 6,6,6,12)]) +
  scale_fill_manual(values = cbPalette_12[c(2,5,4, 3, 7, 7, 6,6,6,12)]) +
  scale_shape_manual(values = rep(c("circle filled", "square filled", "diamond filled", "cross", "circle filled", "plus", "triangle filled", "triangle down filled", "diamond filled"), 5)) +
  labs(x = "Gene (alphabetical)", y = "-log10 P-value") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank())

ggsave("../figures/nmd_alphabetic.pdf", width = 6, height = 4, units = "in")

# Scatter plot to show effect size
plot_data %>%
  ggplot(aes(x = lfc, y = -log10(p_value), color = annotation, fill = annotation, shape = annotation)) +
  facet_wrap(~ fct_rev(reporter), ncol = 1) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color = "grey") +
  scale_color_manual(values = cbPalette_12[c(2,5,4, 3,7, 7, 6,6,6,12, 12)]) +
  scale_fill_manual(values = cbPalette_12[c(2,5,4, 3, 7, 7, 6,6,6,12, 12)]) +
  scale_shape_manual(values = rep(c("circle filled", "square filled", "diamond filled", "cross", "circle filled", "plus", "triangle filled", "triangle down filled", "diamond filled"), 5)) +
  labs(x = "-log2 mRNA level\n(median-centered)", y = "-log10 P-value")

ggsave("../figures/nmd_volcano.pdf", width = 4.5, height = 4.5, units = "in")

# Highlight data for EIF2, EIF3, EIF4F
mageck_sgrna  %>% 
  filter(str_detect(gene, "^EIF2S|^EIF2B|^EIF2B|^EIF3|EIF4A1|EIF4E$|EIF4G1")) %>% 
  filter(str_detect(sample_name, "ptc.+bc1.+bc2")) %>% 
  group_by(gene) %>% 
  summarize(mean_lfc = mean(lfc),  sd_lfc = sd(lfc), n = n()) %>%
  left_join(mageck_gene %>% filter(str_detect(sample_name, "ptc.+bc1.+bc2"))  %>% select(gene, pos_p_value), by = "gene") %>%
  mutate(sig = case_when(
    pos_p_value < 0.001 ~ "***",
    pos_p_value < 0.01 ~ "**",
    pos_p_value < 0.05 ~ "*",
    TRUE ~ "")
  ) %>% 
  write_csv("../../../../source_data/figure_s4d.csv") %>%
  ggplot(aes(x = gene, y = mean_lfc)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  geom_point() +
  geom_errorbar(aes(ymin = mean_lfc - sd_lfc, ymax = mean_lfc + sd_lfc), width = 0.2) +
  geom_text(aes(label = sig, y = mean_lfc + sd_lfc + 0.5), angle = 90) +
  labs(x = "sgRNA", y = "log2 PTC mRNA level\n(w.r.t control mRNA)") +
  coord_flip() +
  theme(axis.text.y = element_text(size = 8))

ggsave("../figures/ptc_mrna_eif2_eif3_eif4.pdf", width = 2.5, height = 4.5, units = "in")

# Highlight fold changes for NMD hit gene groups
plot_data %>%
  filter(reporter == "PTC") %>%
  ggplot(aes(x = annotation, y = lfc, color = annotation)) +
  geom_point(data = plot_data  %>% filter(annotation != "Elongation factor"), show.legend = F, shape = "circle", alpha = 0) +
  geom_boxplot(data = plot_data  %>% filter(annotation == "Other"), 
  show.legend = F, outlier.shape = NA, width = 0.4, alpha = 0.2) +
  geom_point(aes(shape = annotation, fill = annotation),
    size = 1,
    data = plot_data %>% filter(annotation != "Other", pos_fdr < fdr_cutoff),
    position = position_jitter(height = 0, width = 0.2, seed = 111),
    show.legend = F
  ) +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.5, color = "grey") +
  scale_color_manual(values = cbPalette_12[c(2,5,4, 3,7, 7, 6,6,6,12)]) +
  scale_fill_manual(values = cbPalette_12[c(2,5,4, 3, 7, 7, 6,6,6,12)]) +
  scale_shape_manual(values = rep(c("circle filled", "square filled", "diamond filled", "cross", "circle filled", "plus", "triangle filled", "triangle down filled", "diamond filled"), 5)) +
  scale_y_continuous(limits = c(-1.5, 4)) +
  labs(y = "log2 mRNA level", x = "Gene group") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("../figures/nmd_hits_jitter.pdf", width = 2.5, height = 2.5, units = "in")

# NMD + ISRIB / GCN2 hits as volcano plot
fdr_cutoff <- 0.01
plot_data <- mageck_gene %>%
  filter(str_detect(sample_name, "gcn2_|isrib_1")) %>%
  mutate(reporter = str_extract(sample_name, "ptc|ntc")) %>%
  mutate(treatment = str_extract(sample_name, "isrib_1|sggcn2")) %>%
  left_join(nmd_hit_annotations, by = "gene") %>%
  left_join(annotation_order, by = "annotation") %>%
  mutate(priority = if_else(is.na(annotation), as.integer(10), priority)) %>%
  mutate(annotation = if_else(is.na(annotation), "Other", as.character(annotation))) %>%
  mutate(annotation = fct_reorder(annotation, priority)) %>%
  mutate(reporter = toupper(reporter)) %>%
  write_csv("../../../../source_data/figure_4c_4d.csv")

plot_data %>%
  ggplot(aes(x = lfc, y = -log10(neg_p_value))) +
  facet_wrap(~treatment, ncol = 3, scales = "free") +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color = "grey") +
  geom_point(size = 1, alpha = 1, data = plot_data %>% filter(fdr > fdr_cutoff), color = "grey") +
  geom_point(aes(shape = annotation, color = annotation, fill = annotation),
    size = 2, alpha = 1,
    data = plot_data %>% filter(fdr < fdr_cutoff),
    position = position_jitter(height = 0.3, width = 0.3, seed = 11)
  ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_color_manual(values = cbPalette_12[c(2,4, 3,7, 6,6,5,12)]) +
  scale_fill_manual(values = cbPalette_12[c(2,4, 3, 7, 6,6,5,12)]) +
  scale_shape_manual(values = rep(c("circle filled", "diamond filled", "cross", "plus", "triangle filled", "triangle down filled", "diamond filled"), 5)) +
  labs(x = "log2 mRNA ratio\n(treatment / control)", y = "-log10 P-value") +
  theme(legend.title = element_blank())

ggsave("../figures/nmd_gcn1_isrib_volcano.pdf", width = 6, height = 2.5, units = "in")

# Read GORILLA GO terms
gorilla <- read_csv("../data/gorilla/ptc_day7_dmso_bc1_total_vs_ptc_day7_dmso_bc2_total.csv", show_col_types = F) %>% 
  janitor::clean_names() %>% 
  mutate(enrichment = str_extract(enrichment_n_b_n_b, "\\d+\\.\\d+")) %>%
  select(-genes, -enrichment_n_b_n_b, -p_value) %>% 
  write_csv("../../../../source_data/figure_s4c.csv")

# Show selected GO terms enriched in downregulated genes
subset_go_terms <- c(
'GO:0006413' = 'translational initiation',
'GO:0000184' = 'nuclear-transcribed mRNA catabolic process, nonsense-mediated decay',
'GO:0034198' = 'cellular response to amino acid starvation',
'GO:0036499' = 'PERK-mediated unfolded protein response',
'GO:0042254' = 'ribosome biogenesis',
'GO:0044391' = 'ribosomal subunit',
'GO:0005851' = 'eukaryotic translation initiation factor 2B complex',
'GO:0005850' = 'eukaryotic translation initiation factor 2 complex'
)

selected_go_terms <- gorilla  %>% 
  filter(go_term %in% names(subset_go_terms)) %>%
  select(go_term, description, fdr_q_value, enrichment)

print(selected_go_terms)
print("Analysis complete!")