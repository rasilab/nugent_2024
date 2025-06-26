options(warn = -1)
suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(plyranges)
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
))

subset_samples <- c("SF3B5", "SF3B6","AQR", "FLUC")
subset_genes <- c("RPL41", "RPL24")

sample_annotations <- read_csv("../annotations/sample_annotations.csv", show_col_types = F) %>% 
  filter(str_detect(sample_id, "226")) %>%
  select(sample_id, sample_name) %>% 
  mutate(sample_name = toupper(sample_name)) %>% 
  print()

plasmid_sj <- read_tsv("../annotations/plasmids/pHPHS232_pHPHS800_pAS321.ss.tsv", show_col_types = F) %>% 
  mutate(sj_num = row_number()) %>%
  print()

plasmid_junction_counts <- list.files("../data/alignments/",
            recursive = T, full.names = T, pattern = "SJ.out.tab") %>%
  as_tibble_col("file") %>% 
  filter(str_detect(file, "plasmid")) %>%
  mutate(sample_id = str_extract(file, "226p\\d+")) %>% 
  filter(!is.na(sample_id)) %>%
  mutate(data = map(file, readSTARJunctions)) %>%
  select(-file) %>% 
  mutate(data = map(data, as_tibble)) %>% 
  unnest(data) %>%
  left_join(plasmid_sj, by = c("seqnames", "start", "end", "strand")) %>%
  print()


exons <- read_gff2("../data/ensembl/Homo_sapiens.GRCh38.108.exons.gtf") %>%
  print()

sj <- read_tsv("../data/ensembl/Homo_sapiens.GRCh38.108.all.ss.tsv", show_col_types = F) %>% 
  mutate(sj_num = row_number()) %>%
  print()

annotated_sj <- sj %>%
  filter(junction_type == "annotated") %>%
  mutate(row = row_number()) %>%
  GRanges() 

n_introns <- findOverlaps(annotated_sj, annotated_sj, type = "within") %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  filter(query_hits != subject_hits) %>%
  mutate(within_sj_num = annotated_sj[query_hits]$sj_num) %>%
  mutate(within_transcript_id = annotated_sj[query_hits]$transcript_id) %>%
  mutate(intron_start = start(annotated_sj[subject_hits]), intron_end = end(annotated_sj[subject_hits])) %>%
  mutate(within_intron_start = start(annotated_sj[query_hits]), within_intron_end = end(annotated_sj[query_hits])) %>%
  filter(intron_start == within_intron_start | intron_end == within_intron_end) %>%
  group_by(subject_hits, within_transcript_id) %>%
  summarize(n_introns = dplyr::n(), within_sj_nums = list(within_sj_num),.groups = "drop")

n_exons <- findOverlaps(exons, annotated_sj, type = "within") %>%
  as_tibble() %>%
  janitor::clean_names() %>%
  mutate(within_transcript_id = exons[query_hits]$transcript_id) %>%
  group_by(subject_hits, within_transcript_id) %>%
  summarize(n_exons = dplyr::n(), .groups = "drop")

cassette_sj <- inner_join(n_introns, n_exons, by = c("subject_hits", "within_transcript_id")) %>%
  rename(row = subject_hits) %>%
  left_join(as_tibble(annotated_sj), by = "row") %>%
  filter(n_introns == 2, n_exons == 1) %>%
  ungroup() %>%
  select(-row, -n_introns, -n_exons, -n_exon_skip, -junction_type) %>%
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

sj_types <- junction_counts %>% 
  filter(um_reads > 0) %>%
  left_join(cassette_sj %>% select(sj_num, within_transcript_id), by = "sj_num") %>% 
  left_join(cassette_sj  %>% unnest(within_sj_nums) %>% select(within_sj_nums, within_transcript_id), by = c("sj_num" = "within_sj_nums")) %>%
  mutate(sj_type = case_when(
    !is.na(within_transcript_id.x) ~ "skipped_cassette_sj", 
    !is.na(within_transcript_id.y) ~ "included_cassette_sj", 
    TRUE ~ "other"
)) %>% 
  select(sample_name, sj_num, sj_type)

plot_data <- sj_types %>% 
  group_by(sample_name, sj_type) %>%
  count() %>% 
  ungroup() %>%
  group_by(sample_name) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(sj_type = fct_rev(fct_relevel(sj_type, "other")))

options(repr.plot.width = 5, repr.plot.height = 2, repr.matrix.max.rows = 12)

plot_data  %>%
  pivot_wider(names_from = sj_type, values_from = n) 

plot_data %>%
  filter(sj_type != "other") %>%
  ggplot(aes(x = sample_name, y = n, fill = sj_type)) +
  geom_col(position = "fill") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  scale_fill_manual(values = cbPalette_12) +
  labs(y = "Junction fraction", x = "sgRNA")

options(repr.plot.width = 7, repr.plot.height = 2)
plot_data %>%
  filter(sj_type != "other") %>%
  ggplot(aes(x = sample_name, y = n, fill = sj_type)) +
  geom_col(position = "dodge") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  scale_fill_manual(values = cbPalette_12) +
  labs(y = "Num. cassette junctions", x = "sgRNA")

options(repr.plot.width=10, repr.plot.height=6)

plot_data <- junction_counts %>% 
  filter(str_detect(sample_name, "AQR|FLUC|SF3")) %>%
  filter(um_reads > 0) %>% 
  select(sj_num, um_reads, sample_name) %>%  
  complete(sj_num, sample_name, fill = list('um_reads' = 1)) %>% 
  filter(!is.na(sj_num)) %>% 
  group_by(sj_num) %>%
  mutate(um_fluc = um_reads[sample_name == "FLUC"]) %>% 
  ungroup() %>% 
  filter(sample_name != "FLUC", um_reads > 1)

plot_data %>%
  group_by(sample_name) %>% 
  mutate(sample_name = paste0(sample_name, " (N=", dplyr::n(), ")")) %>%
  ungroup() %>% 
  ggplot(aes(x = um_fluc, y = um_reads)) +
  facet_wrap(~sample_name, ncol = 4) +
  geom_point(alpha=0.2) +
  scale_x_log10() +
  scale_y_log10()

options(repr.plot.width=10, repr.plot.height=6)

plot_data %>%
  inner_join(cassette_sj %>% select(sj_num), by = "sj_num") %>%
  group_by(sample_name) %>% 
  mutate(sample_name = paste0(sample_name, " (N=", dplyr::n(), ")")) %>%
  ungroup() %>% 
  ggplot(aes(x = um_fluc, y = um_reads)) +
  facet_wrap(~sample_name, ncol = 4) +
  geom_point(alpha=0.2) +
  scale_x_log10() +
  scale_y_log10()

options(repr.matrix.max.rows = 12, repr.plot.width = 6, repr.plot.height = 4)
plot_data_2 <- plot_data %>% 
  left_join(cassette_sj %>% select(-gene_id), by = "sj_num") %>%
  filter(str_detect(sample_name, "SF3B5|SF3B6")) %>% 
  mutate(lfc = round(log2(um_reads / um_fluc), 2)) %>%
  group_by(sample_name) %>% 
  mutate(lfc = lfc - median(lfc)) %>%
  filter(um_reads > 50, um_fluc > 1) %>%
  arrange(desc(lfc)) %>%
  mutate(rank = seq(1,dplyr::n())) %>%
  ungroup()

plot_data_2 %>%
  filter(rank < 100) %>%
  ggplot(aes(x = rank, y = lfc, color = sample_name)) +
  geom_line(size = 1) +
  scale_color_manual(values = cbPalette_12) +
  geom_text(aes(label = sample_name, y= lfc), x = 15,hjust = 0, vjust = 0,data = plot_data_2 %>% filter(rank == 15))


pso <- cassette_sj %>%
  unnest(within_sj_nums) %>%
  left_join(junction_counts %>% select(sj_num, sample_name, um_reads), by = c("within_sj_nums" = "sj_num"), suffix = c("", ".within")) %>%
  rename(um_reads.within = um_reads) %>%
  filter(um_reads.within >= 100) %>%
  group_by(sample_name, sj_num) %>%
  summarize(um_reads.within = as.integer(mean(um_reads.within)), .groups = "drop") %>%
  left_join(junction_counts %>% select(sj_num, sample_name, um_reads), by = c("sample_name", "sj_num")) %>%
  mutate(um_reads = ifelse(is.na(um_reads), 0, um_reads)) %>%
  mutate(pso = 100 * um_reads / (um_reads.within + um_reads)) %>%
  complete(sample_name, sj_num, fill = list("pso" = NA, "um_reads.within" = NA, "um_reads" = 0)) %>%
  left_join(sj, by = "sj_num") %>%
  print()

options(repr.plot.width = 3, repr.plot.height = 1.6)

plot_data <- pso %>%
  filter(gene_name %in% subset_genes) %>%
  filter(sample_name %in% subset_samples) %>%
  arrange(gene_name, sample_name) %>%
  filter(sj_num %in% c(636761, 1986975)) %>%
  select(sample_name, gene_name, sj_num, pso) %>%
  mutate(sample_name = fct_relevel(sample_name, subset_samples))

plot_data %>%
  ggplot(aes(x = sample_name, y = pso, shape = gene_name, group = gene_name)) +
  geom_point(size = 2) +
  geom_line(color = "#0088aa", size = 1) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) +
  labs(x = "sgRNA", y = "Exon skipped\nisoform (%)") +
  scale_y_continuous(breaks = c(0, 20, 40))

# Export source data
plot_data %>%
  mutate(across(where(is.numeric), ~signif(.x, 3))) %>%
  write_csv("../../../source_data/figure_3i.csv")

ggsave("../figures/rpl24_rpl41_exon_skipped_fraction.png", width = 3, height = 1.6, units = "in")


pso %>%
  group_by(sample_name) %>%
  summarize(mean_pso = mean(pso), n = dplyr::n(), .groups = "drop")

options(repr.plot.width = 2.2, repr.plot.height = 2)
plot_data <- pso %>%
  filter(um_reads > 2) %>%
  group_by(sample_name) %>%
  arrange(desc(pso)) %>%
  mutate(rank = seq(1, dplyr::n())) %>%
  ungroup() %>%
  filter(sample_name %in% c("FLUC", "SF3B5", "SF3B6", "AQR"))

plot_data %>%
  ggplot(aes(x = rank, y = pso, color = sample_name)) +
  geom_line(size = 1) +
  scale_color_manual(values = cbPalette_12) +
  scale_x_continuous(limits = c(1,100)) +
  labs(y = "Skipped isoform (%)", x = "Exon") +
  theme(legend.position = c(0.8, 0.75), legend.title = element_blank(), legend.text = element_text(size = 10), legend.key.size = unit(0.5, "cm"))

options(repr.plot.width = 1.5, repr.plot.height = 2)

plot_data <- pso %>%
  group_by(sj_num) %>%
  mutate(pso_fluc = pso[sample_name == "FLUC"], 
  um_fluc = um_reads[sample_name == "FLUC"],
  um.within_fluc = um_reads.within[sample_name == "FLUC"]
  ) %>% 
  ungroup() %>%
  filter(!is.na(um_reads.within), !is.na(um.within_fluc)) %>%
  mutate(delta_pso = pso - pso_fluc) %>%
  filter(sample_name %in% c("SF3B5", "SF3B6", "AQR")) %>%
  filter(um_reads >= 2, um_fluc >= 2) %>%
  group_by(sample_name) %>%
  arrange(desc(delta_pso)) %>%
  mutate(rank = row_number()) %>%
  ungroup()

plot_data %>%
  ggplot(aes(x = sample_name, y = delta_pso)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(y = "Skipped isoform\nchange (%)", x = "sgRNA")

options(repr.plot.width = 4.2, repr.plot.height = 2)

plot_data %>%
  ggplot(aes(x = um_fluc, y = um_reads)) +
  facet_wrap(~sample_name, ncol = 4) +
  geom_abline(slope = 1, intercept = 0, color = "blue") +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_x_log10() +
  scale_y_log10()
  # theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  # labs(x = "Exon exclusion % (sgFLUC)", y = "Exon exclusion %")

plot_data %>%
  group_by(sample_name) %>%
  nest() %>%
  ungroup() %>%
  mutate(ttest = map(data, ~t.test(.x$pso, .x$pso_fluc, paired = T, alternative = "greater"))) %>%
  mutate(ttest = map(ttest, broom::tidy)) %>%
  mutate(wilcoxtest = map(data, ~wilcox.test(.x$pso, .x$pso_fluc, paired = T, alternative = "greater"))) %>%
  mutate(wilcoxtest = map(wilcoxtest, broom::tidy)) %>%
  select(-data) %>%
  unnest(ttest, wilcoxtest) 


options(repr.plot.width = 1.75, repr.plot.height = 1.8)

plot_data %>%
  mutate(sample_name = fct_relevel(sample_name, subset_samples)) %>% 
  ggplot(aes(x = rank, y = delta_pso, color = sample_name)) +
  geom_line(size = 1) +
  scale_color_manual(values = cbPalette_12[c(5,6,7,12)]) +
  # scale_x_continuous(limits = c(1,100)) +
  scale_y_continuous(breaks = c(0, 30, 60), limits = c(NA, 60)) +
  labs(y = "Isoform change (%)\nw.r.t control sgRNA", x = "Exon rank") +
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

ggsave("../figures/skipped_exon_isoform_change.png", width = 1.75, height = 1.8, units = "in")

plot_data %>% 
  filter(delta_pso >= 10) %>% 
  group_by(sample_name) %>%
  summarize(n = dplyr::n(), .groups = "drop")
