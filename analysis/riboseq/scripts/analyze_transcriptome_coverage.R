options(warn = -1, repr.matrix.max.rows = 15)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicAlignments)
  library(plyranges)
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

# Load annotations and data
sample_annotations <- read_csv("../annotations/sample_annotations.csv", show_col_types = F) %>% 
  filter(str_detect(sample_name, "mono")) %>% 
  extract(sample_name, c("sample_name", "replicate"), "(.*)_(.*)") %>% 
  select(sample_id, sample_name, replicate)

tx_annotations <- read_gff2("ftp://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz") 
seqlevelsStyle(tx_annotations) <- "UCSC"

tx_gene_names <- tx_annotations %>% 
  filter(type == "transcript") %>% 
  as_tibble() %>% 
  select(transcript_id, gene_name, transcript_name)

tx <- tx_annotations %>% 
  filter(type == "exon") %>% 
  filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
  split(.$transcript_id)

aln <- list.files("../data/alignments/", pattern = ".bam$", full.names = T) %>% 
  as_tibble_col("file") %>% 
  mutate(sample_id = str_extract(file, "225p\\d+")) %>% 
  inner_join(sample_annotations, by = "sample_id") %>% 
  mutate(aln = map(file, . %>% readGAlignments %>% as_data_frame)) %>% 
  select(-file) %>% 
  unnest(aln)

start_codons_tx <- tx_annotations %>% 
  filter(type == "start_codon") %>% 
  filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))) %>%
  GenomicFeatures::mapToTranscripts(tx) %>%
  as_tibble() %>% 
  select(seqnames, start, end)

utr_cds_annotations <- tx_annotations %>% 
  filter(type %in% c("UTR", "CDS")) %>%
  filter(seqnames %in% paste0("chr", c(1:22, "X", "Y", "M"))) 

utr_cds_length <- utr_cds_annotations %>%
  GenomicFeatures::mapToTranscripts(tx) %>%
  as_tibble() %>% 
  mutate(type = utr_cds_annotations$type[xHits]) %>%
  left_join(start_codons_tx, by = "seqnames", suffix = c("", "_start")) %>%
  left_join(tx_gene_names, by = c("seqnames" = "transcript_id")) %>%
  mutate(type = case_when(
    type == "UTR" & start < start_start ~ "UTR5",
    type == "UTR" & start > start_start ~ "UTR3",
    TRUE ~ "CDS"
  )) %>%
  group_by(gene_name, type) %>%
  summarize(start = min(start), .groups = "drop")

# Prepare plot data for Figure 5f (metadensity)
plot_data <- aln %>% 
  left_join(start_codons_tx, by = "seqnames", suffix = c("", "_start")) %>% 
  filter(!is.na(start_start)) %>% 
  mutate(relative_start = start - start_start) %>% 
  filter(relative_start < 500, relative_start > -50) %>% 
  filter(width > 20, width < 40) %>% 
  group_by(sample_name, width, relative_start) %>% 
  count()

# Generate Figure 5f (metadensity plot)
options(repr.plot.width = 4.5, repr.plot.height = 1.8)

fig5f_data <- plot_data  %>% 
  filter(width >= 27, width <= 33, relative_start < 200) %>%
  group_by(sample_name , relative_start) %>%
  summarize(n = sum(n), .groups = "drop") %>%
  group_by(sample_name) %>% 
  mutate(n = n / max(n)) %>%
  ungroup() %>% 
  mutate(sample_name = str_replace(sample_name, "_mono", ""))

fig5f_plot <- fig5f_data %>% 
  ggplot(aes(x = relative_start + 13, y = n, color = sample_name)) +
  scale_color_manual(values = cbPalette_12[c(12,3,6,1)]) +
  geom_line() +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "P-site position w.r.t start codon", y = "Normalized read count")

ggsave("../figures/riboseq_metadensity.pdf", plot = fig5f_plot, width = 4.5, height = 1.8)

# Save Figure 5f source data
write_csv(fig5f_data, "../../../source_data/figure_5f.csv")

# Prepare alignment counts for gene-specific plots
aln_counts <- aln %>% 
  GRanges() %>% 
  filter(width >= 27, width <= 33) %>% 
  narrow(start = 15, width = 1) %>%
  as_tibble() %>% 
  group_by(seqnames, start, end, sample_name) %>%
  summarize(score = sum(width), .groups = "drop") %>%
  ungroup() %>% 
  mutate(psite = start) %>%
  GRanges()

genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
tx_seq <- GenomicFeatures::extractTranscriptSeqs(genome, tx)

# Generate Figure 5g - JUN panel
options(repr.plot.width = 4, repr.plot.height = 2)

gene <- "JUN"
plot_data_jun <- aln_counts  %>%
  as_tibble() %>%
  left_join(tx_gene_names, by = c("seqnames" = "transcript_id")) %>%
  filter(gene_name == gene)

tx_width_jun <- tx_gene_names %>% 
  filter(gene_name == gene) %>% 
  pull(transcript_id) %>% 
  tx_seq[.] %>% 
  width()

plot_data_jun <- tx_width_jun %>% 
  seq(1, .) %>% 
  as_tibble_col(column_name = "start") %>% 
  left_join(plot_data_jun, by = "start") %>% 
  complete(sample_name, start) %>% 
  mutate(score = ifelse(is.na(score), 0, score)) %>% 
  filter(!is.na(sample_name)) %>% 
  mutate(drug = str_extract(sample_name, "dmso|hht")) %>%
  mutate(sgrna = str_extract(sample_name, "gcn1|fluc"))

gene_plot_jun <- plot_data_jun  %>% 
  inner_join(utr_cds_length, by = c("psite" = "start", "gene_name"))

jun_plot <- plot_data_jun %>%
  ggplot(aes(x = psite, y = score, fill = sample_name, color = sample_name)) +
  facet_wrap(~ drug + sgrna, scales = "free_y", ncol = 1) +
  geom_col() +
  geom_point(aes(x = psite, y = 0), size = 1, shape = "plus", color = "red", data = gene_plot_jun, vjust = 1, hjust = 0.5) +
  labs(title = gene) +
  theme(strip.text.x = element_blank(),
    axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank()
  ) +
  scale_fill_manual(values = cbPalette_12[c(12,3,6,7)]) +
  scale_color_manual(values = cbPalette_12[c(12,3,6,7)]) +
  scale_x_continuous(limits=c(1, tx_width_jun))

ggsave("../figures/jun_riboseq.pdf", plot = jun_plot, width = 4, height = 2)

# Save JUN data for Figure 5g
write_csv(plot_data_jun, "../../../source_data/figure_5g_jun.csv")

# Generate Figure 5g - MYC panel
gene <- "MYC"
plot_data_myc <- aln_counts  %>%
  as_tibble() %>%
  left_join(tx_gene_names, by = c("seqnames" = "transcript_id")) %>%
  filter(gene_name == gene)

tx_width_myc <- tx_gene_names %>% 
  filter(gene_name == gene) %>% 
  pull(transcript_id) %>% 
  tx_seq[.] %>% 
  width()

plot_data_myc <- tx_width_myc %>% 
  seq(1, .) %>% 
  as_tibble_col(column_name = "start") %>% 
  left_join(plot_data_myc, by = "start") %>% 
  complete(sample_name, start, gene_name) %>% 
  mutate(score = ifelse(is.na(score), 0, score)) %>% 
  filter(!is.na(sample_name), !is.na(gene_name)) %>% 
  mutate(drug = str_extract(sample_name, "dmso|hht")) %>%
  mutate(sgrna = str_extract(sample_name, "gcn1|fluc"))

gene_plot_myc <- plot_data_myc  %>% 
  inner_join(utr_cds_length, by = c("start", "gene_name"))

myc_plot <- plot_data_myc %>%
  ggplot(aes(x = start, y = score, fill = sample_name, color = sample_name)) +
  facet_wrap(~ drug + sgrna, scales = "free_y", ncol = 1) +
  geom_col() +
  geom_point(aes(x = start, y = 0), size = 1, shape = "plus", color = "red", data = gene_plot_myc, vjust = 1, hjust = 0.5) +
  labs(title = gene) +
  theme(strip.text.x = element_blank(),
    axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  axis.line.y = element_blank()
  ) +
  scale_fill_manual(values = cbPalette_12[c(12,3,6,7)]) +
  scale_color_manual(values = cbPalette_12[c(12,3,6,7)])

ggsave("../figures/myc_riboseq.pdf", plot = myc_plot, width = 4, height = 2)

# Save MYC data for Figure 5g
write_csv(plot_data_myc, "../../../source_data/figure_5g_myc.csv")