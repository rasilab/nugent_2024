#!/usr/bin/env Rscript
# Extended Data Figure 5c qPCR panel

# Load libraries
options(warn = -1)
suppressPackageStartupMessages({
  library(tidyverse)
  library(rasilabRtemplates)
})

# Set theme
theme_set(theme_rasilab())

# Define color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Create figures directory if it doesn't exist
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}

# Read annotations
annotations <- read_csv("../annotations/sampleannotations.csv", show_col_types = FALSE) %>% 
  mutate(
    Well = as.character(Well),
    Plate = as.character(Plate)
  )

# Read Cq data from Cfx Duet
cq_data <- list.files("../../../../data/qpcr/u937_gcn1_hht/", full.names = TRUE) %>% 
  as_tibble_col("file") %>%
  mutate(data = map(file, ~ read_csv(.x, show_col_types = FALSE) %>% as_tibble())) %>%
  mutate(Plate = str_extract(file, "(?<=plate)\\d")) %>%
  unnest("data") %>%
  select(Cq, Well, Plate) %>%
  rename(CT = Cq) %>%
  inner_join(annotations, by = c("Well", "Plate")) %>%
  filter(!is.na(CT)) %>%
  # Remove outlier
  mutate(platewell = str_c(Plate, Well)) %>%
  filter(platewell != "7F07") %>%
  select(CT, Amplicon, sgRNA, Treatment) %>%
  mutate(Treatment = fct_relevel(Treatment, "DMSO"))

# Normalize expression to GAPDH
norm_ct_inverted <- cq_data %>%
  filter(!is.na(CT)) %>%
  group_by(Amplicon, Treatment, sgRNA) %>%
  summarize(
    std_error = sd(CT),
    n = n(),
    CT = mean(CT),
    .groups = "drop"
  ) %>%
  group_by(Treatment, sgRNA) %>%
  mutate(
    norm_ct = CT - CT[Amplicon == "GAPDH"],
    norm_error = sqrt(std_error^2 + std_error[Amplicon == "GAPDH"]^2)
  ) %>%
  filter(Amplicon != "GAPDH") %>%
  mutate(norm_ct = 10 - norm_ct)

# Export source data
write_csv(norm_ct_inverted, "../../../../source_data/figure_s5c.csv")

# Create plot
p <- norm_ct_inverted %>% 
  ggplot(aes(x = Treatment, y = norm_ct, ymax = norm_ct + norm_error, ymin = norm_ct - norm_error, color = sgRNA)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  facet_wrap(~Amplicon, scales = "free_y") +
  ylab("Relative mRNA level\n(-ΔΔCT)") +
  theme(
    axis.text.x = element_text(),
    axis.line = element_line(color = "grey")
  ) +
  scale_color_manual(values = cbPalette)

# Save figure
ggsave("../figures/egr1_jun_mrna_levels_s5c.pdf", p, height = 2.2, width = 6)

options(warn = 0)