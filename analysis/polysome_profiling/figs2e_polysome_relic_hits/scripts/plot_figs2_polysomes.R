#!/usr/bin/env Rscript
# Extended Data Figure 2e polysome profiling panels

# Load libraries
options(warn = -1)
suppressPackageStartupMessages({
  library(tidyverse)
  library(rasilabRtemplates)
})

# Set theme
theme_set(theme_rasilab() +
  theme(
    axis.line = element_line(color = "grey"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )
)

# Define color palette
cbPalette_12 <- c(
  "#999999", "#CC6677", "#88CCEE", "#661100", "#117733", "#332288", 
  "#999933","#AA4499", "#44AA99", "#882255", "#6699CC", "#DDCC77"
)

# Create figures directory if it doesn't exist
if (!dir.exists("../figures")) {
  dir.create("../figures", recursive = TRUE)
}

# Read annotations
annotations <- read_csv("../annotations/sampleannotations.csv", col_types = cols(.default = "c"))

# Read data
counts <- list.files("../../../../data/polysome_profiling/polysome_relic_hits/", pattern = "\\d\\.csv", full.names = TRUE) %>% 
  enframe("sno", "file") %>% 
  mutate(data = map(file, ~ read_csv(.x, skip = 45, show_col_types = FALSE) %>%
                      select(`Position(mm)`, Absorbance, `Fraction Number`)
                    )) %>%
  unnest(data) %>%
  rename(
    AbsA = Absorbance,
    Fraction = `Fraction Number`,
    Position = `Position(mm)`
  ) %>%
  mutate(file = str_extract(file, "20241104_i287_\\d.csv")) %>%
  left_join(annotations, by = "file") %>%
  filter(!is.na(sample)) %>%
  select(-sno, -file)

# Standardize traces
standardized <- counts %>%
  # remove all initial values < 0
  filter(!(AbsA < 0 & Position < 10)) %>%
  # make baseline equal to 0
  group_by(sample) %>%
  mutate(AbsA = AbsA - min(AbsA)) %>%
  ungroup()

# Calculate P/M ratio
PbyM <- standardized %>%
  group_by(sample) %>%
  # Pick a minimum position just short of 40S and 60S peaks
  summarize(
    M = sum(AbsA[Position > 10 & Position < 28]),
    P = sum(AbsA[Position >= 28 & Position < 80]),
    .groups = "drop"
  ) %>% 
  mutate(PbyM = P/M) %>% 
  mutate(PbyM = PbyM / PbyM[sample == "sgFLUC"]) %>% 
  mutate(PbyM = round(PbyM, 2)) %>% 
  arrange(desc(PbyM))

print("P/M ratios:")
print(PbyM)

# Normalize to approximate value of 80S
normalized_counts <- standardized %>%
  group_by(sample) %>%
  mutate(AbsA = AbsA - min(AbsA)) %>%
  mutate(MonoAbsA = max(AbsA[Position > 10])) %>%
  mutate(AbsA = AbsA/MonoAbsA)

# Prepare data for plotting and export
plot_data <- normalized_counts %>%
  filter(Position > 10 & Position < 80, !sample %in% c("sgRPL5")) %>% 
  left_join(PbyM, by = "sample") %>%
  mutate(sample = paste0(sample, ", ", PbyM))

# Export source data
write_csv(plot_data, "../../../../source_data/figure_s2e.csv")

# Create plot
p <- plot_data %>%
  ggplot(aes(x = Position, y = AbsA, color = sample)) +
  geom_line(linewidth = 0.8) +
  xlab("Position (mm)") +
  ylab("A260 (A.U.)") +
  labs(color = "sgRNA, P / M") +
  scale_y_continuous(limits = c(0, 1.2), breaks = c(0, 0.4, 0.8, 1.2)) +
  scale_color_manual(values = cbPalette_12[c(2, 1, 3, 4, 5)])

# Save figure
ggsave("../figures/normalized_trace_with_PbyM_s2e.pdf", p, width = 6, height = 3, units = "in")

options(warn = 0)