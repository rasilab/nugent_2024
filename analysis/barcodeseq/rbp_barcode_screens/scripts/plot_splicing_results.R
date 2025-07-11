#!/usr/bin/env Rscript

# Analyze splicing ReLiC screen
# Converted from plot_splicing_results.ipynb

# Define libraries and analysis-specific parameters
options(warn = -1, repr.matrix.max.rows = 10)

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rasilabRtemplates))

cbPalette_12 <- c(
  "#88CCEE", "#CC6677", "#117733", "#999933", "#332288", "#AA4499",
  "#661100", "#44AA99", "#882255", "#6699CC", "#DDCC77", "#888888"
)

theme_set(theme_rasilab() + 
 theme(
  axis.line = element_line(color = "grey"), 
 axis.title.y = element_text(margin = margin(r=10)),
 axis.title.x = element_text(margin = margin(t=10))
))

fdr_cutoff <- 0.05

# Load MaGeCK gene hit data
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
  arrange(pos_rank) 

print(mageck_gene)

# Load MaGeCK sgrna hit data
mageck_sgrna <- list.files("../data/mageck/", full.names = T, recursive = T, pattern = "sgrna_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc"), str_detect(file, "retained|skipped")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = F) %>%
  mutate(data = map(file, . %>% read_tsv(show_col_types = F))) %>%
  select(-file) %>%
  unnest(data) %>%
  janitor::clean_names()

print(mageck_sgrna)

# Plot number of hits on different days (Figure S3a)
plot_data <- mageck_gene  %>% 
  group_by(treatment, control) %>%
  summarize(n_fdr = sum(pos_fdr < fdr_cutoff), .groups = "drop") %>% 
  filter(str_detect(treatment, "day")) %>%
  mutate(day = str_extract(treatment, "(?<=day)\\d+")) %>%
  mutate(isoform = str_extract(treatment, "e2|i1|i2(?=.+)")) %>% 
  mutate(isoform = case_when(
    isoform == "e2" ~ "e13",
    isoform == "i1" ~ "i12",
    isoform == "i2" ~ "i23",
    TRUE ~ "other"
  )) %>% 
  mutate(isoform = fct_relevel(isoform, "i23", "i12", "e13", "other")) %>%
  write_csv("../../../../source_data/figure_s3a.csv")  

plot_data %>% 
  ggplot(aes(x = day, y = n_fdr, fill = isoform, group = isoform, shape = isoform)) +
  geom_point(size = 2) +
  geom_line() +
  labs(x = "Days post Cas9", y = "Number of genes \n with FDR < 0.25") +
  theme(legend.position = "top", legend.direction = "horizontal", legend.title = element_blank())

ggsave("../figures/splicing_n_hits.pdf", width = 2.5, height = 2.7, units = "in")

# Show the number of hits as a table
options(repr.matrix.max.rows = 12)
print(plot_data)

# Load Splicing hit annotations
splicing_hit_annotations <- read_csv("../annotations/splicing_hit_annotations.csv", show_col_types = FALSE) %>%
  print()

# Volcano plot of P value vs fold change for different isoforms (Figure 3e)
plot_data <- mageck_gene %>%
  filter(str_detect(treatment, "day")) %>%
  mutate(day = str_replace(str_extract(treatment, "day\\d+"), "day", "day ")) %>%
  mutate(isoform = str_extract(treatment, "(e2|i1|i2)_.+")) %>%
  group_by(isoform) %>% 
  mutate(fdr_cutoff = fdr_cutoff) %>%
  ungroup() %>% 
  left_join(splicing_hit_annotations, by = "gene") %>%
  mutate(annotation = fct_relevel(annotation, "Spliceosome/splicing-associated", "SF3a/b", "Translation", "RNA exosome")) %>%
  arrange(annotation) %>%
  write_csv("../../../../source_data/figure_3e.csv")

plot_data %>%
  ggplot(aes(x = lfc, y = -log10(pos_p_value), color = annotation, shape = isoform)) +
  facet_grid(fct_rev(isoform) ~ day) +
  geom_point(size = 1, alpha = 0.5, data = plot_data %>% filter(is.na(annotation) | pos_fdr > fdr_cutoff), color = "grey", shape = "circle") +
  geom_point(alpha = 1, data = plot_data %>% filter(!is.na(annotation), pos_fdr < fdr_cutoff)) +
  scale_x_continuous(breaks= c(-3, 0, 3)) +
  scale_y_continuous(breaks = c(0, 2, 4)) +
  scale_color_manual(values = cbPalette_12[c(11,7,2,3,5)]) +
  scale_shape_manual(values = c(15,17,16,18,3)) +
  labs(x = "log2 isoform / total", y = "-log10 P-value")

ggsave("../figures/splicing_volcano.pdf", width = 8, height = 4, units = "in")

# Check EIF3 genes
options(repr.matrix.max.rows = 19)
eif3_results <- plot_data %>% 
  filter(str_detect(gene, "EIF3"))  %>% 
  select(gene, day, isoform, lfc, pos_fdr) %>% 
  arrange(desc(pos_fdr)) %>% 
  filter(pos_fdr < 0.05) %>% 
  arrange(desc(gene))

print(eif3_results)

# LFC for SF3b complex (Figure 3f)
plot_data <- mageck_gene %>%
  filter(str_detect(treatment, "day")) %>%
  mutate(day = str_extract(treatment, "(?<=day)\\d+")) %>%
  mutate(isoform = str_extract(treatment, "(e2|i1|i2)_.+")) %>% 
  mutate(gene = if_else(gene == "PHF5A", "SF3B7", gene)) %>% 
  filter(str_detect(gene, "^SF3B|AQR")) %>% 
  mutate(sig = case_when(
                pos_fdr >= 0.05 ~ 1,
                pos_fdr < 0.05 ~ 2,
                )) %>% 
  mutate(gene = fct_relevel(gene, "SF3B1", "SF3B2", "SF3B3", "SF3B4", "SF3B5", "SF3B6", "SF3B7", "AQR", "XAB2", "ISY1", "PPIE")) %>%
  write_csv("../../../../source_data/figure_3f.csv")

plot_data %>%
  ggplot(aes(x = day, y = lfc, color = gene, group = gene, size = as.factor(sig), alpha = as.factor(sig), shape = fct_rev(isoform))) +
  facet_wrap(~fct_rev(isoform), ncol = 3, scales = "fixed") +
  geom_point(position = position_jitter(width = 0.2, height = 0, seed = 111)) +
  scale_color_manual(values = cbPalette_12[c(1,8,3,4,5,6,9,7,2,10,11)]) +
  scale_size_manual(values = c(0.5, 1.5)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(x = "Days post Cas9", y = "isoform / total") +
  theme(legend.box = "horizontal", legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 11))

ggsave("../figures/sf3b_lfc.pdf", width = 7, height = 1.8, units = "in")

# Load gene ontology analysis from GORILLA
subset_go_terms <- c(
  "GO:0005681" = 1,
  "GO:0046540" = 1.5,
  "GO:0005686" = 2,
  "GO:0022625" = 3,
  "GO:0006413" = 4,
  "GO:0051170" = 5,
  NULL
)

gorilla <- list.files("../data/gorilla/", full.names = T, pattern = "retained|skipped") %>% 
  read_csv(show_col_types = F, id = "filename") %>% 
  janitor::clean_names() %>% 
  mutate(isoform= str_extract(filename, "i2|i1|e2")) %>% 
  mutate(enrichment = as.numeric(str_extract(enrichment_n_b_n_b, "\\d+\\.\\d+"))) %>%
  select(isoform, go_term, description, fdr_q_value, enrichment, ontology) %>% 
  filter(go_term %in% names(subset_go_terms)) %>%
  print()

# Plot enrichment of GO terms (Figure 3d)
plot_data <- gorilla %>%
  complete(isoform, nesting(go_term, description), fill = list(fdr_q_value = NA, enrichment = NA)) %>% 
  mutate(y = str_c(description, " ", go_term)) %>% 
  mutate(y = fct_reorder(y, -subset_go_terms[go_term])) %>%
  write_csv("../../../../source_data/figure_3d.csv")

plot_data %>%
  ggplot(aes(x = fct_rev(isoform), y = y, size = enrichment, color = y)) +
  geom_point(shape = 16) +
  scale_color_manual(values = cbPalette_12[c(3,7,2,5,4,11)]) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0),
  panel.grid.minor = element_line(size = 1, color = "grey"),
  legend.text = element_text(size = 10),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "top",
  axis.text = element_text(size = 11, color = "black"),
  axis.title.x = element_blank(), axis.title.y = element_blank()) +
  guides(color = F) +
  scale_x_discrete(position = "top")

ggsave("../figures/splicing_go_enrichment.pdf", width = 4, height = 2.4, units = "in")

cat("Analysis complete. Generated figures:\n")
cat("- ../figures/splicing_n_hits.pdf (Figure S3a)\n")
cat("- ../figures/splicing_volcano.pdf (Figure 3e)\n") 
cat("- ../figures/sf3b_lfc.pdf (Figure 3f)\n")
cat("- ../figures/splicing_go_enrichment.pdf (Figure 3d)\n")