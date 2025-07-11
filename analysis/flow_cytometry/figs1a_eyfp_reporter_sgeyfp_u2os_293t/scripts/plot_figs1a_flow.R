#!/usr/bin/env Rscript
# Extended Data Figure 1a flow cytometry panel

# Load libraries
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rasilabRtemplates))
suppressPackageStartupMessages(library(ggridges))

# Set theme
theme_set(theme_rasilab())

# Define color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Create figures directory if it doesn't exist
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}

# Define data folders
fcs_file_folder <- c(
  "../../../../data/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/day1/", 
  "../../../../data/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/day3/", 
  "../../../../data/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/day5/",
  "../../../../data/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/day8/",
  "../../../../data/flow_cytometry/figs1a_eyfp_reporter_sgeyfp_u2os_293t/day11/"
)

# Nice channel names
channels <- c(
  "fitc_a" = "yfp", 
  "pe_texas_rd_a" = "rfp",
  "bv421_a" = "bfp"
)

# Read annotations
annotations <- read_csv("../annotations/sample_annotations.csv", show_col_types = FALSE) %>% 
  mutate(file = as.character(file))

# Read flow data
options(warn = -1)

flowdata <- list.files(fcs_file_folder, full.names = TRUE, recursive = TRUE, pattern = ".fcs") %>% 
  as_tibble_col("file") %>%
  mutate(day = as.integer(str_extract(file, "(?<=day)\\d+"))) %>% 
  mutate(data = map(file, ~ read.FCS(.x, transformation = FALSE, alter.names = TRUE) %>% exprs() %>% as_tibble())) %>%
  mutate(file = str_extract(file, "(?<=events_)[:graph:]+(?=.fcs$)")) %>%
  unnest("data") %>%
  janitor::clean_names() %>%
  rename_with(~ channels[.x], names(channels)) %>%
  select(day, file, ssc_a, bfp, yfp, rfp)

# Join data with annotations
data <- flowdata %>% 
  inner_join(annotations, by = "file")

# Process and save data
plot_data <- data %>%
  pivot_longer(c("bfp", "yfp", "rfp", "ssc_a"), names_to = "channel") %>% 
  filter(channel %in% c("yfp")) %>% 
  mutate(sgRNA = case_when(
    sgRNA == "parent" ~ "Parent",
    sgRNA == "FLUC" ~ "sgCTRL",
    sgRNA == "YFP" ~ "sgEYFP"
  )) %>%
  mutate(channel = case_when(
    channel == "yfp" ~ "EYFP"
  )) %>%
  filter(value >= 1)

# Export source data
write_csv(plot_data, "../../../../source_data/figure_s1a.csv")

# Create plot
p <- plot_data %>%
  ggplot(aes(x = value, y = as.factor(day), fill = sgRNA)) +
  geom_density_ridges(alpha = 0.8) +
  facet_wrap(~cell_line) +
  scale_x_log10(
    limits = c(1, 1e6), 
    breaks = scales::trans_breaks("log10", function(x) 100^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_fill_manual(values = cbPalette[c(1, 3, 2)]) +
  labs(x = "Fluorescence (a.u.)", y = "Days post Cas9", fill = "") +
  theme(
    axis.line = element_line(color = "grey"),
    axis.ticks = element_line(color = "grey")
  )

# Save figure
ggsave("../figures/sgyfp_sgfluc_effects_for_validation_s1a.pdf", p, width = 5.5, height = 2.5)

options(warn = 0)