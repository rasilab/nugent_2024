# Plot mRNA vs gDNA fold-change colored by gene essentiality

options(warn = -1)
suppressPackageStartupMessages({
  library(tidyverse)
  library(rasilabRtemplates)
})

theme_set(theme_rasilab() +
  theme(
    axis.line = element_line(color = "grey"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )
)

# colorblind-friendly colors (Okabe-Ito)
essential_colors <- c("non-essential" = "#0072B2", "essential" = "#D55E00")

# read source data
plot_data <- read_csv("../../../../source_data/figure_1d.csv", show_col_types = F)

# read essential gene list
essential_genes <- read_csv("../annotations/DepMapCRISPRInferredCommonEssentials_23Q2.csv",
  show_col_types = F) %>%
  separate(Essentials, c("gene", "num")) %>%
  select(gene) %>%
  mutate(essential = TRUE)

# join and classify
plot_data <- plot_data %>%
  left_join(essential_genes, by = "gene") %>%
  mutate(essential = if_else(is.na(essential), "non-essential", "essential")) %>%
  mutate(day = fct_relevel(day, "day 5", "day 13", "day 21")) %>%
  arrange(desc(essential))  # plot non-essential first, essential on top

# compute correlations per day
correlation <- plot_data %>%
  group_by(day) %>%
  nest() %>%
  mutate(cor = map(data, ~cor.test(~ mRNA + gDNA, data = .x, method = "spearman", exact = F))) %>%
  mutate(cor = map(cor, broom::tidy)) %>%
  select(-data) %>%
  unnest(cor)

plot_data %>%
  ggplot(aes(x = gDNA, y = mRNA, color = essential)) +
  facet_wrap(~day) +
  geom_point(size = 0.5, alpha = 0.3) +
  geom_text(data = correlation, aes(label = paste0("r = ", round(estimate, 2))),
    x = -4.5, y = 2.5, hjust = 0, vjust = 1, size = 3, color = "black", inherit.aes = FALSE) +
  scale_x_continuous(breaks = c(-6, -3, 0, 3), limits = c(-5, 3)) +
  scale_y_continuous(breaks = c(-6, -3, 0, 3), limits = c(-6, 3)) +
  scale_color_manual(values = essential_colors) +
  labs(
    x = "gDNA fold-change (log2, median-centered)",
    y = "mRNA fold-change\n(log2, median-centered)",
    color = "") +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"))

ggsave("../figures/compare_mrna_gdna_foldchange_with_time_essential.png",
  width = 4, height = 2.5, dpi = 96)
ggsave("../figures/compare_mrna_gdna_foldchange_with_time_essential.pdf",
  width = 4, height = 2.5)
