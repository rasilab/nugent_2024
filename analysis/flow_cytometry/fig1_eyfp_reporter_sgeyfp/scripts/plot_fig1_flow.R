# Flow cytometry analysis, Fig. 1b

options(warn = -1)

# Load required libraries
suppressPackageStartupMessages({
  library(flowCore)       # For reading .fcs files
  library(broom)          # For extracting statistical fits
  library(tidyverse)      # For analysis and plotting
  library(rasilabRtemplates) # Lab ggplot2 theme and color-blind palette
})

# Define data folder and channel names
fcs_file_folder <- c(
  "../../../../data/flow_cytometry/fig1_eyfp_reporter_sgeyfp/day2/", 
  "../../../../data/flow_cytometry/fig1_eyfp_reporter_sgeyfp/day5/", 
  "../../../../data/flow_cytometry/fig1_eyfp_reporter_sgeyfp/day7/"
)

channels <- c(
  "fitc_a" = "yfp", 
  "pe_texas_rd_a" = "rfp",
  "bv421_a" = "bfp"
)

if (!file.exists('../figures')){
    dir.create('../figures')
}

# Read in annotations
annotations <- read_csv("../annotations/sampleannotations.csv", show_col_types = F) %>% 
  print()

# Read in flow data
flowdata <- list.files(fcs_file_folder, full.names = TRUE, pattern = ".fcs") %>%
  as_tibble_col("file") %>%
  mutate(
    data = map(file, ~ read.FCS(.x, transformation = FALSE, alter.names = TRUE) %>%
      exprs %>%
      as_tibble)
  ) %>%
  unnest(data) %>%
  janitor::clean_names() %>%
  rename_with(~ channels[.x], names(channels)) %>%
  mutate(
    day = str_extract(file, "(?<=day).") %>% as.integer(),
    file = str_extract(file, "(?<=events_)[:graph:]+(?=.fcs$)") %>% as.integer()
  ) %>%
  select(day, file, ssc_a, bfp, yfp, rfp) %>%
  print()

# Join data with annotations
data <- flowdata %>% 
     inner_join(annotations, by = c("file")) %>% 
     print() 

# Plot data
plot_data <- data %>%
  pivot_longer(c("bfp", "yfp", "rfp", "ssc_a"), names_to = "channel") %>%
  filter(channel == "yfp", value >= 1) %>%
  mutate(
    sgrna = recode(sgrna, "none" = "Background", "fluc" = "sgCTRL", "yfp" = "sgEYFP"),
    channel = recode(channel, "yfp" = "EYFP"),
    value = round(value, 2)
  ) %>%
  write_csv("../../../../source_data/figure_1b.csv")

plot_data %>%
  ggplot(aes(x = value, y = as.factor(day), fill = sgrna)) +
  ggridges::geom_density_ridges(alpha = 0.8) +
  scale_x_log10(
    limits = c(1, 1e5),
    breaks = scales::trans_breaks("log10", function(x) 100^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_fill_manual(values = cbPalette[c(1, 3, 2)]) +
  labs(x = "Fluorescence (a.u.)", y = "Days post Cas9", fill = "") +
  theme(
    axis.line = element_line(color = "grey"),
    axis.ticks = element_line(color = "grey")
  )

ggsave("../figures/sgyfp_sgfluc_effects_for_validation.pdf", width = 3.5, height = 2)