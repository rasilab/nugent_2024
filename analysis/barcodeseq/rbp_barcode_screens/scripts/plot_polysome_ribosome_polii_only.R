#!/usr/bin/env Rscript

# Standalone plot: 60S/40S ribosomal subunits and RNA Pol II perturbations

options(warn = -1)

suppressPackageStartupMessages({
  library(tidyverse)
})

cbPalette_12 <- c(
  "#DDCC78", "#CC6677", "#6699CC", "#661100", "#117733", "#999933",
  "#332288", "#AA4499", "#44AA99", "#882255", "#88CCEE", "#999999"
)

# Use rasilabRtemplates theme if available; otherwise fallback
if (requireNamespace("rasilabRtemplates", quietly = TRUE)) {
  theme_set(rasilabRtemplates::theme_rasilab() +
    theme(
      axis.line = element_line(color = "grey"),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 10))
    ))
} else {
  theme_set(theme_bw() +
    theme(
      axis.line = element_line(color = "grey"),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 10))
    ))
}

# Prefer DejaVu Sans for all text
theme_update(text = element_text(family = "DejaVu Sans"))

# Enable DejaVu Sans embedding if available
if (requireNamespace("showtext", quietly = TRUE)) {
  if (requireNamespace("sysfonts", quietly = TRUE)) {
    fams <- tryCatch(sysfonts::font_families(), error = function(e) character())
    if (!("DejaVu Sans" %in% fams)) {
      paths <- c(
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans.ttf"
      )
      path <- paths[file.exists(paths)][1]
      if (!is.na(path)) try(sysfonts::font_add(family = "DejaVu Sans", regular = path), silent = TRUE)
    }
  }
  try(showtext::showtext_auto(), silent = TRUE)
}

fdr_cutoff <- 0.05
set.seed(111)

# Load MaGeCK gene-level results (polysome screens only)
mageck_gene <- list.files("../data/mageck/", full.names = TRUE, recursive = TRUE, pattern = "gene_summary.tsv$") %>%
  as_tibble_col("file") %>%
  filter(str_detect(file, "ntc.+polysome")) %>%
  mutate(sample_name = str_extract(file, "(?<=mageck//).+(?=/mageck.)")) %>%
  separate(sample_name, c("treatment", "control"), sep = "_vs_", remove = FALSE) %>%
  mutate(rep = str_extract(sample_name, "rep1|rep2")) %>%
  mutate(treatment = str_extract(treatment, "heavy|light|mono|supernatant")) %>%
  mutate(control = str_extract(control, "mono|supernatant")) %>%
  filter(!is.na(treatment), !is.na(control), !is.na(rep)) %>%
  mutate(data = map(file, ~read_tsv(.x, show_col_types = FALSE))) %>%
  select(-file) %>%
  unnest(data) %>%
  select(-`neg|lfc`) %>%
  rename(lfc = `pos|lfc`, gene = id) %>%
  mutate(p_value = if_else(lfc > 0, `pos|p-value`, `neg|p-value`)) %>%
  mutate(fdr = if_else(lfc > 0, `pos|fdr`, `neg|fdr`)) %>%
  mutate(goodsgrna = if_else(lfc > 0, `pos|goodsgrna`, `neg|goodsgrna`)) %>%
  mutate(fdr = if_else(is.na(fdr), 1, fdr))

# Prepare data for the focused plot
plot_data <- mageck_gene %>%
  filter(control == "mono", rep == "rep1") %>%
  mutate(annotation = case_when(
    str_detect(gene, "^RPL") ~ "60S ribosome",
    str_detect(gene, "^RPS|RACK1") ~ "40S ribosome",
    str_detect(gene, "^POLR2") ~ "RNA polymerase II",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(annotation)) %>%
  # Invert group order relative to default
  mutate(annotation = factor(annotation, levels = c("RNA polymerase II", "40S ribosome", "60S ribosome"))) %>%
  mutate(sig = if_else(fdr < fdr_cutoff, "sig", "notsig"))

# Plot and save PNG
plot_data %>%
  filter(treatment == "heavy") %>%
  ggplot(aes(x = annotation, y = lfc, fill = annotation, color = annotation)) +
  geom_hline(aes(yintercept = 0), linetype = 2, linewidth = 0.5, color = "grey") +
  geom_point(aes(shape = sig), size = 1, position = position_jitter(height = 0, width = 0.3, seed = 11)) +
  scale_color_manual(values = cbPalette_12[c(6, 2, 7)]) +
  scale_fill_manual(values = cbPalette_12[c(6, 2, 7)]) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_shape_manual(values = c("circle open", "circle")) +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  labs(y = "log2 (heavy \n polysome / monosome)") +
  coord_flip()

devfun <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
ggsave("../figures/polysome_ribosome_polII_groups.pdf", width = 2, height = 2, units = "in", device = devfun, family = "DejaVu Sans")

message("Wrote ../figures/polysome_ribosome_polII_groups.pdf")
