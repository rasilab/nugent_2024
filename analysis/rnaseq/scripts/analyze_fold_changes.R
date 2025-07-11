## Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(plyranges)
  library(rasilabRtemplates)
  library(DESeq2)
})

count_cutoff <- 1000

gene_subset <- c("FOS", "JUN", "ATF3", "EGR1", "DUSP1", "IER2", "MYC", "GPR50", "ADAMST1", "RHOB", "TIMP3", "JUND")

theme_set(theme_rasilab() +
  theme(
    axis.line = element_line(color = "grey"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  ))

## Read sample annotations
sample_annotations <- read_csv("annotations/sample_annotations.csv", show_col_types = F) %>%
  filter(str_detect(sample_id, "^220"))

## Read gene annotations
gene_annotations <- read_gff2("data/ensembl/Homo_sapiens.GRCh38.108.gtf") %>%
  as_tibble() %>%
  filter(type == "gene", gene_biotype == "protein_coding")

## Read in gene count data
count_data <- "data/alignments/" %>% 
  list.files(full.names = T, pattern = "ReadsPerGene.out.tab", recursive = T) %>% 
  as_tibble_col("file") %>% 
  filter(str_detect(file, "genome_220")) %>%
  mutate(sample_id = str_extract(file, "220p[:digit:]+")) %>%
  mutate(data = map(file, ~read_tsv(.x, show_col_types = F, col_names = c("gene_id", "total", "fwd", "rev")))) %>% 
  select(-file) %>%
  unnest(data) %>%
  select(-total, -fwd) %>%
  filter(str_detect(gene_id, "^ENSG"))

## Subset to transcripts that get summed counts across conditions greater than count cutoff
count_data_above_cutoff <- count_data %>% 
  group_by(gene_id) %>% 
  mutate(total_study_counts = sum(rev)) %>% 
  ungroup() %>% 
  filter(total_study_counts >= count_cutoff) %>% 
  select(-total_study_counts)

## Count data for GCN1-Harringtonine experiment
count_data <- count_data_above_cutoff %>% 
  left_join(sample_annotations %>% select(sample_id, sample_name), by = "sample_id") %>%
  select(-sample_id) %>%
  pivot_wider(values_from = "rev", names_from = "sample_name") %>% 
  column_to_rownames("gene_id")

## Col data for GCN1-Harringtonine experiment
col_data <- colnames(count_data) %>% 
  as_tibble_col("sample_name") %>% 
  separate(sample_name, c("ko", "sgrna1", "sgrna2", "rep", "time"), remove=F) %>%
  unite(sample, ko, sgrna1, sgrna2, time, sep = "_") %>%
  column_to_rownames("sample_name")

## Run DESeq2 for GCN1-Harringtonine experiment
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ sample)
dds <- DESeq(dds)

## Get DESeq2 results for specific comparisons
results_fluc <- results(dds, 
    contrast = c("sample", "fluc_1_2_6h", "fluc_1_2_0h")) %>%
  as.data.frame() %>% 
  dplyr::rename(lfc = log2FoldChange) %>%
  mutate(across(baseMean:pvalue, . %>% round(3))) %>% 
  arrange(-lfc) %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  inner_join(gene_annotations %>% select(gene_id, gene_name), by = "gene_id")

results_gcn1_1_2 <- results(dds, 
    contrast = c("sample", "gcn1_1_2_6h", "gcn1_1_2_0h")) %>%
  as.data.frame() %>% 
  dplyr::rename(lfc = log2FoldChange) %>%
  mutate(across(baseMean:pvalue, . %>% round(3))) %>% 
  arrange(-lfc) %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  inner_join(gene_annotations %>% select(gene_id, gene_name), by = "gene_id")

results_gcn1_3_4 <- results(dds, 
    contrast = c("sample", "gcn1_3_4_6h", "gcn1_3_4_0h")) %>%
  as.data.frame() %>% 
  dplyr::rename(lfc = log2FoldChange) %>%
  mutate(across(baseMean:pvalue, . %>% round(3))) %>% 
  arrange(-lfc) %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  inner_join(gene_annotations %>% select(gene_id, gene_name), by = "gene_id")

## Combine results
lfc <- bind_rows(
  results_fluc %>% mutate(compare = "sgfluc_1_2.6hvs0h"),
  results_gcn1_1_2 %>% mutate(compare = "sggcn1_1_2.6hvs0h"),
  results_gcn1_3_4 %>% mutate(compare = "sggcn1_3_4.6hvs0h")
)

## Scatter plot of harringtonine fold changes between sgGCN1 and sgFLUC
plot_data <- lfc %>% 
  filter(compare %in% c("sgfluc_1_2.6hvs0h", "sggcn1_1_2.6hvs0h", "sggcn1_3_4.6hvs0h")) %>%
  mutate(compare = str_extract(compare, "^[^\\.]+")) %>%
  mutate(compare = str_replace(compare, "sgfluc_1_2", "FLUC sgRNA 1,2")) %>%
  mutate(compare = str_replace(compare, "sggcn1_1_2", "GCN1 sgRNA 1,2")) %>%
  mutate(compare = str_replace(compare, "sggcn1_3_4", "GCN1 sgRNA 3,4")) %>%
  mutate(padj = if_else(padj < 1e-50, 1e-50, padj)) %>% 
  select(compare, padj, lfc, gene_id, gene_name, baseMean, lfcSE) %>% 
  filter(!is.na(gene_name))

plot_data_2 <- plot_data %>% 
  mutate(lfc = if_else(lfc < -5, -5, lfc)) %>%
  group_by(gene_id) %>% 
  filter(dplyr::n() == 3) %>%
  mutate(lfc_ctrl = lfc[compare == "FLUC sgRNA 1,2"]) %>% 
  ungroup() %>% 
  filter(compare != "FLUC sgRNA 1,2")

subset_data_2 <- plot_data_2 %>% 
  filter(gene_name %in% gene_subset)

plot_data_2 %>%
  ggplot(aes(x = lfc_ctrl, y = lfc)) +
  facet_wrap(~compare, ncol = 2) +
  geom_point(alpha = 0.2, color = "grey") +
  geom_point(data = subset_data_2, stroke = NA, size = 2, color = "black") +
  ggrepel::geom_text_repel(data = subset_data_2, aes(label = gene_name), color = "black", size = 3, max.overlaps = Inf) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "log2 mRNA fold-change (HHT/CTRL)\nsgFLUC", y = "log2 mRNA fold-change (HHT/CTRL)\nsgGCN1")

ggsave("figures/scatter_plot_hht_treatment_gcn1_vs_fluc.pdf", width = 4.5, height = 2.8)

## IEG alone plot
plot_data <- lfc %>% 
  filter(compare %in% c("sgfluc_1_2.6hvs0h", "sggcn1_1_2.6hvs0h", "sggcn1_3_4.6hvs0h")) %>%
  mutate(compare = str_extract(compare, "^[^\\.]+")) %>%
  mutate(compare = str_replace(compare, "sgfluc_1_2", "FLUC")) %>%
  mutate(compare = str_replace(compare, "sggcn1_1_2", "GCN1_1_2")) %>%
  mutate(compare = str_replace(compare, "sggcn1_3_4", "GCN1_3_4")) %>%
  filter(compare != "GCN1_3_4") %>% 
  mutate(sgrna = str_replace(compare, "GCN1_1_2", "GCN1")) %>%
  select(sgrna, padj, lfc, gene_id, gene_name, baseMean, lfcSE) %>% 
  filter(gene_name %in% gene_subset) %>% 
  group_by(gene_name) %>% 
  mutate(order = sum(lfc)) %>% 
  ungroup() %>% 
  mutate(gene_name = fct_reorder(gene_name, -order))

plot_data %>% 
  ggplot(aes(x = gene_name, y = lfc, ymin = lfc - lfcSE, ymax = lfc + lfcSE, color = sgrna)) +
  geom_point() +
  geom_errorbar(width = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5), legend.position = "top") +
  scale_y_continuous(breaks = c(0,4,8)) +
  labs(y = "log2 mRNA fold-change\nHHT /  CTRL", x = "Gene")

ggsave("figures/ieg_alone_plot_hht_treatment_gcn1_vs_fluc.pdf", width = 2.7, height = 2)