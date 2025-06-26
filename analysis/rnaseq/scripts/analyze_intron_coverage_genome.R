options(warn = -1)
suppressPackageStartupMessages({
  library(Biostrings)
  library(plyranges)
  library(tidyverse)
  library(GenomicAlignments)
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
))

read_length <- 120
overlap_min <- 20

subset_samples <- c("SF3B5", "SF3B6","AQR", "FLUC")
subset_genes <- c("RPL41", "RPL24")

sample_annotations <- read_csv("../annotations/sample_annotations.csv", show_col_types = F) %>% 
  filter(str_detect(sample_id, "226")) %>%
  select(sample_id, sample_name) %>% 
  mutate(sample_name = toupper(sample_name)) %>% 
  print()

sj <- read_tsv("../data/ensembl/Homo_sapiens.GRCh38.108.all.ss.tsv", show_col_types = F) %>% 
  mutate(sj_num = row_number()) %>%
  filter(junction_type == "annotated") %>%
  # select(sj_num, gene_id, gene_name, transcript_id) %>%
  print()

junction_counts <- list.files("../data/alignments/",
            recursive = T, full.names = T, pattern = "SJ.out.tab") %>%
  as_tibble_col("file") %>% 
  filter(str_detect(file, "genome.+226p")) %>%
  mutate(data = map(file, readSTARJunctions)) %>%
  mutate(data = map(data, as_tibble)) %>% 
  mutate(sample_id = str_extract(file, "226p\\d+")) %>% 
  select(-file) %>% 
  unnest(data) %>%
  left_join(sj, by = c("seqnames", "start", "end", "strand")) %>%
  right_join(sample_annotations, by = "sample_id") %>% 
  print()

intron_counts <- list.files("../data/intron_counts/",
            recursive = T, full.names = T, pattern = "genome_.*csv") %>%
  as_tibble_col("file") %>% 
  mutate(sample_id = str_extract(file, "226p\\d+")) %>% 
  filter(!is.na(sample_id)) %>%
  mutate(data = map(file, ~read_csv(.x, show_col_types = F))) %>%
  select(-file) %>% 
  unnest(data) %>%
  rename(sj_num = intron, intron_count = count) %>%
  left_join(sj, by = "sj_num") %>%
  mutate(gene_id = str_extract(gene_id, "ENSG\\d+")) %>%
  print()

gene_counts <- list.files("../data/alignments/",
            recursive = T, full.names = T, pattern = "ReadsPerGene.out.tab") %>%
  as_tibble_col("file") %>% 
  filter(str_detect(file, "genome")) %>%
  mutate(sample_id = str_extract(file, "226p\\d+")) %>% 
  filter(!is.na(sample_id)) %>%
  mutate(data = map(file, ~read_tsv(.x, col_names = c("gene_id", "total", "fwd", "rev"), show_col_types = F))) %>%
  select(-file) %>% 
  unnest(data) %>%
  select(-total,-fwd) %>%
  mutate(gene_id = str_extract(gene_id, "ENSG\\d+")) %>%
  filter(!is.na(gene_id)) %>%
  print()

psi <- intron_counts %>% 
  left_join(junction_counts %>% select(um_reads, sample_id, sj_num), by = c("sample_id", "sj_num")) %>%
  rename(sj_count = um_reads) %>%
  filter(sj_count >= 100) %>%
  mutate(width = end - start + 1) %>%
  mutate(intron_count = intron_count / (width + read_length) * read_length) %>%
  mutate(psi = 100 * intron_count / (sj_count + intron_count)) %>%
  left_join(sample_annotations, by = "sample_id") %>%
  print()

options(repr.plot.width = 3, repr.plot.height = 1.6)

plot_data <- psi %>%
  filter(gene_name %in% subset_genes) %>%
  filter(sample_name %in% subset_samples) %>%
  arrange(gene_name, sample_name) %>%
  filter(sj_num %in% c(636757, 1986993)) %>%
  select(sample_name, gene_name, sj_num, psi) %>%
  mutate(sample_name = fct_relevel(sample_name, subset_samples))

plot_data %>%
  ggplot(aes(x = sample_name, y = psi, shape = gene_name, group = gene_name)) +
  geom_point(size = 2) +
  geom_line(color = "#44aa00", size = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) +
  labs(x = "sgRNA", y = "Intron retained\nisoform (%)")

# Export source data
plot_data %>%
  mutate(across(where(is.numeric), ~signif(.x, 3))) %>%
  write_csv("../../../source_data/figure_3i.csv")

ggsave("../figures/rpl24_rpl41_intron_retained_fraction.png", width = 3, height = 1.6, units = "in")



plot_data

delta_psi <- psi %>%
  filter(str_detect(sample_name, "AQR|SF3B5|SF3B6|FLUC")) %>%
  select(sample_name, sj_num, psi, intron_count, sj_count) %>%
  complete(sample_name, sj_num, fill = list(psi = NA, intron_count = 0, sj_count = NA)) %>%
  group_by(sj_num) %>%
  mutate(psi_fluc = psi[sample_name == "FLUC"], 
  intron_count_fluc = intron_count[sample_name == "FLUC"],
  sj_count_fluc = sj_count[sample_name == "FLUC"],
  ) %>% 
  ungroup() %>%
  filter(!is.na(sj_count), !is.na(sj_count_fluc)) %>%
  filter(intron_count >= 20, intron_count_fluc >= 20) %>%
  mutate(delta_psi = psi - psi_fluc) %>%
  filter(sample_name %in% c("SF3B5", "SF3B6", "AQR")) %>%
  group_by(sample_name) %>%
  arrange(desc(delta_psi)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>% 
  left_join(sj, by = "sj_num")

plot_data <- delta_psi 

plot_data %>%
  ggplot(aes(x = sample_name, y = delta_psi)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(y = "Intron-retained isoform\nchange (%)", x = "sgRNA")

delta_psi %>%
  group_by(sample_name) %>%
  nest() %>%
  ungroup() %>%
  mutate(ttest = map(data, ~t.test(.x$psi, .x$psi_fluc, paired = T, alternative = "greater"))) %>%
  mutate(ttest = map(ttest, broom::tidy)) %>%
  mutate(wilcoxtest = map(data, ~wilcox.test(.x$psi, .x$psi_fluc, paired = T, alternative = "greater"))) %>%
  mutate(wilcoxtest = map(wilcoxtest, broom::tidy)) %>%
  select(-data) %>%
  unnest(ttest, wilcoxtest) 

options(repr.plot.width = 3.5, repr.plot.height = 2)

plot_data %>%
  ggplot(aes(x = rank, y = delta_psi, color = sample_name)) +
  geom_line(size = 2) +
  scale_color_manual(values = cbPalette_12) +
  # scale_x_continuous(limits = c(1,100)) +
  scale_y_continuous(limits = c(-20, 60)) +
  labs(y = "Intron-retained isoform (%)", x = "Intron") +
  theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(0.5, "cm"))

options(repr.plot.width = 1.75, repr.plot.height = 1.8)

plot_data %>%
  mutate(sample_name = fct_relevel(sample_name, subset_samples)) %>% 
  ggplot(aes(x = rank, y = delta_psi, color = sample_name)) +
  geom_line(size = 1) +
  scale_color_manual(values = cbPalette_12[c(5,6,7,12)]) +
  scale_x_continuous(limits = c(1,200), breaks = scales::pretty_breaks(n = 3)) +
  # scale_y_continuous(limits = c(-20, 40), breaks = scales::pretty_breaks(n = 4)) +
  labs(y = "Isoform change (%)\nw.r.t control sgRNA", x = "Intron rank") +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.key.size = unit(0.5, "cm"),
        legend.position = "top",
        legend.direction = "horizontal",
)

# Export source data
plot_data %>%
  mutate(across(where(is.numeric), ~signif(.x, 3))) %>%
  write_csv("../../../source_data/figure_3g.csv")

ggsave("../figures/retained_intron_isoform_change.png", width = 1.75, height = 1.8, units = "in")

plot_data %>% 
  filter(delta_psi >= 10) %>% 
  group_by(sample_name) %>%
  summarize(n = dplyr::n(), .groups = "drop")

