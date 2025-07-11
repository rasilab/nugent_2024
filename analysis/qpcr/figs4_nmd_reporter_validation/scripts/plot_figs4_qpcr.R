#!/usr/bin/env Rscript
# Extended Data Figure 4a qPCR panel

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
annotations <- read_csv("../annotations/sampleannotations.csv", show_col_types = FALSE)

# Read Ct data and join with annotations
cq_data <- read_csv("../../../../data/qpcr/nmd_reporter_validation/cq_values.csv", show_col_types = FALSE) %>%
  select(Cq, Well) %>%
  rename(CT = Cq) %>%
  inner_join(annotations, by = "Well") %>%
  select(CT, Amplicon, Treatment, Reporter) %>%
  mutate(Treatment = fct_relevel(Treatment, "DMSO"))

# Normalize β-globin expression to mCherry
norm_ct_inverted <- cq_data %>%
  filter(!is.na(CT)) %>%
  group_by(Amplicon, Treatment, Reporter) %>%
  summarize(
    std_error = sd(CT),
    n = n(),
    CT = mean(CT),
    .groups = "drop"
  ) %>%
  group_by(Treatment, Reporter) %>%
  mutate(
    norm_ct = CT - CT[Amplicon == "mCherry"],
    norm_error = sqrt(std_error^2 + std_error[Amplicon == "mCherry"]^2)
  ) %>%
  filter(Amplicon == "globin") %>%
  mutate(norm_ct = -5 - norm_ct)

# Export source data
write_csv(norm_ct_inverted, "../../../../source_data/figure_s4a.csv")

# Create plot
p <- norm_ct_inverted %>% 
  ggplot(aes(x = Treatment, y = norm_ct, ymax = norm_ct + norm_error, ymin = norm_ct - norm_error, color = Reporter)) +
  geom_point() +
  geom_errorbar(width = 0.2) +
  ylab("log2 mRNA level\n(-ΔΔCT)") +
  theme(
    axis.text.x = element_text(),
    axis.line = element_line(color = "grey")
  ) +
  scale_x_discrete(guide = guide_axis(angle = 30)) +
  scale_color_manual(values = cbPalette)

# Save figure
ggsave("../figures/mcherry_normalized_ct_values_inverted_s4a.pdf", p, height = 2, width = 3)

options(warn = 0)