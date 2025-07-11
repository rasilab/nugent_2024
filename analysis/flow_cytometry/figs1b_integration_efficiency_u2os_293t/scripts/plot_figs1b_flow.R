#!/usr/bin/env Rscript
# Extended Data Figure 1b flow cytometry panel

# Load libraries
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(broom))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rasilabRtemplates))

# Set theme
theme_set(theme_rasilab())

# Create figures directory if it doesn't exist
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}

# Define data folder
fcs_file_folder <- "../../../../data/flow_cytometry/figs1b_integration_efficiency_u2os_293t/"

# Nice channel names
channels <- c(
  "fitc_a" = "yfp", 
  "pe_texas_red_a" = "rfp",
  "bv421_a" = "bfp"
)

# Read annotations
annotations <- read_csv("../annotations/sample_annotations.csv", show_col_types = FALSE) %>% 
  mutate(file = as.character(file))

# Read flow data
flowdata <- list.files(fcs_file_folder, full.names = TRUE, pattern = ".fcs") %>% 
  as_tibble_col("file") %>%
  mutate(data = map(file, ~ read.FCS(.x, transformation = FALSE, alter.names = TRUE) %>% exprs() %>% as_tibble())) %>%
  mutate(file = str_extract(file, "(?<=001_)[:graph:]+(?=.fcs$)")) %>%
  unnest("data") %>%
  janitor::clean_names() %>%
  rename_with(~ channels[.x], names(channels)) %>%
  select(-time, -yfp, -fsc_a)

# Join data with annotations
data <- flowdata %>%
  inner_join(annotations, by = "file") %>%
  mutate(cell_line = fct_relevel(cell_line, "293T", "U2OS"))

# Export source data
write_csv(data, "../../../../source_data/figure_s1b.csv")

# Create plot
options(warn = -1)

p <- data %>%
  mutate(color = if_else(bfp < 950 & rfp > 1200, "red", "black")) %>%
  ggplot(aes(x = bfp, y = rfp, color = color, size = color)) +
  facet_grid(fct_rev(sgrna) ~ cell_line, scales = "fixed", labeller = "label_value") +
  geom_point() +
  scale_size_manual(values = c(0.1, 0.5)) +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(
    axis.line = element_line(color = "grey"), 
    axis.ticks = element_line(color = "grey")
  ) +
  guides(color = "none", size = "none") +
  labs(x = "BFP (a.u.)", y = "RFP (a.u.)")

# Save figure
ggsave("../figures/bfp_vs_rfp_s1b.pdf", p, width = 3, height = 3, units = "in")

# Calculate integration efficiency statistics
top_left <- data %>%
  filter(rfp > 1200 & bfp < 950) %>% 
  group_by(file, cell_line, sgrna) %>%
  count() %>%
  mutate(top_left_percent = n/100) %>%
  arrange(cell_line, sgrna)

print("Integration efficiency (% in top left quadrant):")
print(top_left)

options(warn = 0)