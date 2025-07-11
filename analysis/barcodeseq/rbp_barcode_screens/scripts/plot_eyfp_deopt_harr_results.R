#!/usr/bin/env Rscript

# Analyze EYFP harringtonine screen
# Converted from plot_eyfp_deopt_harr_results.ipynb

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

# Load MaGeCK gene-level results for EYFP screens
mageck_gene <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "eyfp")) %>%
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

# Load MaGeCK sgRNA-level results for EYFP screens
mageck_sgrna <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "sgrna_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "eyfp")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names()

# Plot volcano plot
plot_data <- mageck_gene %>%
  filter(str_detect(sample_name, "harr")) %>%
  filter(!str_detect(gene, "MCHERRY|PURO")) %>%
  write_csv("../../../../source_data/figure_5b.csv")

plot_data %>%
  ggplot(aes(x = lfc, y = -log10(p_value))) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, color = "grey") +
  geom_point(size = 1, alpha = 1, data = plot_data %>% filter(fdr > fdr_cutoff | lfc < 0), color = "grey", fill = "grey") +
  geom_point(
    size = 2, alpha = 1,
    data = plot_data %>% filter(fdr < fdr_cutoff & lfc > 0), show.legend = F
  ) +
  geom_text(aes(label = gene), data = plot_data %>% filter(fdr < fdr_cutoff & lfc > 0), nudge_x = 0, nudge_y = -0.5, size = 3) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_color_manual(values = cbPalette_12[c(3,12)]) +
  scale_fill_manual(values = cbPalette_12[c(3,12)]) +
  scale_shape_manual(values = rep(c("circle filled", "square filled", "diamond filled", "circle", "plus", "triangle filled", "triangle down filled", "diamond filled"), 5)) +
  labs(x = "log2 mRNA ratio\n(HHT / DMSO)", y = "-log10 P-value") +
  theme(legend.title = element_blank())

ggsave("../figures/eyfp_harr_volcano.pdf", width = 2.4, height = 2.5, units = "in")

# Plot sgRNA effects for specific genes
plot_data <- mageck_sgrna %>%
  filter(str_detect(sample_name, "harr")) %>%
  filter(str_detect(gene, "GCN1|EIF2AK4|ZNF598|GIGYF2|DDX6$|NLUC$")) %>%
  mutate(gene = if_else(gene == "EIF2AK4", "GCN2", gene)) %>%
  mutate(gene = fct_relevel(gene, "GCN1", "GCN2", "ZNF598", "GIGYF2", "DDX6", "NLUC")) %>%
  mutate(sgrna = case_when(
    str_detect(sgrna, "1_2") ~ 1,
    str_detect(sgrna, "2_3") ~ 2,
    str_detect(sgrna, "3_4") ~ 3,
    str_detect(sgrna, "4_1") ~ 4,
  ))

mean_data <- plot_data %>% 
    group_by(gene) %>% 
    summarize(lfc_se = sd(lfc) / sqrt(n()), lfc = mean(lfc)) %>%
    write_csv("../../../../source_data/figure_5c.csv")

plot_data %>%
  ggplot(aes(x = gene, y = lfc, color = gene)) +
  geom_jitter(height = 0, width = 0.2, size = 1, show.legend = F) +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.5, color = "grey") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(-0.5, 1)) +
  scale_color_manual(values = cbPalette_12[c(3,4,5,6,7,12)]) +
  labs(y = "log2 mRNA ratio\n(HHT / DMSO)", x = "sgRNA")

ggsave("../figures/gcn1_sgrna_hht.pdf", width = 2.5, height = 1.8, units = "in")

# Calculate t-test statistics
test_results <- plot_data %>%
  select(sgrna, gene, lfc) %>%
  group_by(gene) %>%
  nest() %>%
  mutate(t = map_dbl(data, ~ t.test(.$lfc, mu = 0)$p.value))

print(test_results)
print("Analysis complete!")