options(warn = -1)
suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(GenomicRanges)
  library(plyranges)
  library(tidyverse)
  library(ggbio)
  library(rasilabRtemplates)
  
})

experiment <- "226"
subset_samples <- c("FLUC")
subset_genes <- c("globin")
subset_transcripts <- c(
  "globin" = c("globin")
)

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

sample_annotations <- read_csv("../annotations/sample_annotations.csv", show_col_types = F) %>% 
  filter(str_detect(sample_id, experiment)) %>%
  select(sample_id, sample_name) %>% 
  mutate(sample_name = toupper(sample_name)) %>% 
  filter(sample_name %in% subset_samples) %>% 
  print()

gtf <- read_gff("../annotations/plasmids/pHPHS232_pHPHS800_pAS321.cleaned.gtf") %>%
  as_tibble() %>%
  mutate(transcript_name = tolower(transcript_id), gene_name = tolower(gene_id)) %>%
  arrange(transcript_name) %>%
  GRanges() %>% 
  filter((transcript_name %in% subset_transcripts) | type == "gene" & gene_name %in% subset_genes) %>%
  print()

gene <- gtf %>%
  filter(type == "transcript") %>%
  mutate(type = "gene") %>% 
  select(gene_name, gene_id)

tx  <- gtf  %>%
  filter(type == "transcript") 

tx_grl <- gtf  %>%
  filter(type == "exon") %>%
  split(.$transcript_name) 

introns <- psetdiff(tx, tx_grl)

gene

tx_with_introns <- introns  %>%
  unlist(use.names = T) %>%
  mutate(transcript_name = names(.), type = "gap") %>%
  unname() %>%
  c(gtf) %>%
  mutate(type = if_else(str_detect(type, "utr"), "utr", as.character(type))) %>%
  mutate(type = tolower(type)) %>%
  filter(type %in% c("exon", "gap")) %>%
  filter(!is.na(gene_name))

options(repr.plot.width = 3, repr.plot.height = 1)

goi <- "globin"
plot_data <- tx_with_introns %>%
  keepSeqlevels(seqlevelsInUse(.)) %>%
  filter(gene_name == goi) %>% 
  split(.$transcript_name) 

start <- gene %>%
  filter(gene_name == goi) %>%
  start()
 
end <- gene %>%
  filter(gene_name == goi) %>%
  end()

p1a <- plot_data %>%
  ggplot() +
  geom_alignment(truncate.gaps = T, fill = "#333333", gap.geom = "chevron") +
  scale_x_continuous(limits = c(start, end)) +
  # scale_y_reverse() +
  theme(
    strip.text.x = element_blank(),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
  ) 

p1a

scan_gr <- ScanBamParam(which = invertStrand(gene))

get_cvg_as_tibble <- function(file, param) {
  readGAlignmentPairs(file, param = param) %>% 
  coverage() %>% 
  GRanges() %>% 
  filter(score > 0) %>% 
  keepSeqlevels(seqlevelsInUse(.)) %>%
  group_by_overlaps(gene) %>%
  as_tibble() %>% 
  mutate(start = map2(start, width, ~seq(.x, .x + .y - 1))) %>% 
  unnest(start) %>% 
  mutate(end = start)
}

cvg <- list.files("../data/alignments/", pattern = "Aligned.out.bam$", recursive = T, full.names = T) %>% 
  as_tibble_col("file") %>% 
  filter(str_detect(file, "plasmid.+226")) %>% 
  mutate(sample_id = str_extract(file, "226p\\d+")) %>% 
  inner_join(sample_annotations, by = "sample_id") %>%
  mutate(cvg = map(file, get_cvg_as_tibble, scan_gr)) %>% 
  select(-file) %>%
  unnest(cvg) 

options(repr.plot.width = 3, repr.plot.height = 1.5)

goi <- "globin"

start <- gene %>%
  filter(gene_name == goi) %>%
  start()
 
end <- gene %>%
  filter(gene_name == goi) %>%
  end()

plot_data <- cvg %>%
  filter(gene_name == goi) %>%
  mutate(sample_name = fct_relevel(sample_name, subset_samples))

label_data <- plot_data %>%  
  mutate(score = max(score)) %>% 
  group_by(sample_name) %>% 
  slice(1)

p2a <- plot_data %>% 
  ggplot(aes(x = start, y = score, fill = sample_name)) +
  facet_wrap(~sample_name, ncol = 1) +
  geom_area() +
  geom_text(aes(label = sample_name), x = end, hjust = 1, vjust  = 1, data = label_data) +
  theme(
    strip.text.x = element_blank(),
    panel.spacing.y = unit(0.1, "lines"),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
  ) +
  guides(fill = F) +
  scale_fill_manual(values = "#333333") +
  scale_x_continuous(limits = c(start, end))

p2a

options(repr.plot.width = 2.5, repr.plot.height = 1.5)
cowplot::plot_grid(p2a,p1a@ggplot, ncol = 1, align = "h", axis = "lr", rel_heights = c(3,1))

# Export source data
plot_data %>%
  mutate(across(where(is.numeric), ~signif(.x, 3))) %>%
  write_csv("../../../source_data/figure_3b.csv")

cowplot::ggsave2("../figures/globin_cvg.png", width = 2.5, height = 1.5)
