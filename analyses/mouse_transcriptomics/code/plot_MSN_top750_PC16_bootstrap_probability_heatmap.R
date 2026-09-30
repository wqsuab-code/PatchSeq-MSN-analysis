#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(viridisLite)
})

root <- file.path(
  "C:/Users/53461/Documents/Codex/2026-07-13/zh/outputs",
  "05_branch_specific_subtype_mapping"
)
source_dir <- file.path(root, "50_MSN_top750_PC16_gene_bootstrap_500")
out_dir <- file.path(root, "82_MSN_top750_PC16_bootstrap_probability_heatmap")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subtypes <- c(paste0("D1_", 1:8), paste0("D2_", 1:8))
vote_file <- file.path(
  source_dir, "MSN_top750_PC16_gene_bootstrap_500_vote_fractions.csv"
)
per_cell_file <- file.path(
  source_dir, "MSN_top750_PC16_gene_bootstrap_500_per_cell.csv"
)

votes <- fread(vote_file)
per_cell <- fread(per_cell_file)
stopifnot(
  nrow(votes) == 588L,
  uniqueN(votes$cell_id) == 588L,
  all(subtypes %in% names(votes)),
  nrow(per_cell) == 588L
)

votes[, cell_id := as.character(cell_id)]
per_cell[, cell_id := as.character(cell_id)]
votes <- merge(
  votes,
  per_cell[, .(
    cell_id,
    bootstrap_top1,
    bootstrap_top1_fraction,
    bootstrap_top2,
    bootstrap_top2_fraction,
    bootstrap_margin,
    bootstrap_entropy,
    stability_class,
    RPCA_agrees_bootstrap_top1
  )],
  by = "cell_id",
  all.x = TRUE,
  sort = FALSE
)

votes[, confidence_class := fifelse(
  bootstrap_top1_fraction >= 0.80 & bootstrap_margin >= 0.50,
  "High",
  fifelse(
    bootstrap_top1_fraction >= 0.50 & bootstrap_margin >= 0.20,
    "Moderate",
    "Low"
  )
)]
votes[, confidence_class := factor(
  confidence_class,
  levels = c("High", "Moderate", "Low")
)]
votes[, top1_order := match(bootstrap_top1, subtypes)]
setorder(
  votes,
  confidence_class,
  top1_order,
  -bootstrap_top1_fraction,
  cell_id
)
votes[, cell_rank := .I]
votes[, tier_rank := seq_len(.N), by = confidence_class]
votes[, tier_n := .N, by = confidence_class]

long <- melt(
  votes,
  id.vars = c(
    "cell_id", "cell_rank", "tier_rank", "tier_n", "bootstrap_top1",
    "bootstrap_top1_fraction", "bootstrap_top2", "bootstrap_top2_fraction",
    "bootstrap_margin", "bootstrap_entropy", "stability_class",
    "confidence_class"
  ),
  measure.vars = subtypes,
  variable.name = "Subtype",
  value.name = "Bootstrap_probability"
)
long[, Subtype := factor(Subtype, levels = subtypes)]
long[, Bootstrap_probability_plot := Bootstrap_probability]

group_info <- votes[, .(
  n_cells = .N,
  y_min = min(tier_rank),
  y_max = max(tier_rank),
  y_mid = mean(range(tier_rank))
), by = .(confidence_class, bootstrap_top1, top1_order)]
setorder(group_info, confidence_class, top1_order)
group_boundaries <- group_info[, head(.SD, -1), by = confidence_class]
tier_info <- votes[, .(n_cells = .N), by = confidence_class]
setorder(tier_info, confidence_class)
tier_labels <- setNames(
  paste0(
    c("High", "Moderate", "Low"),
    " confidence (n=",
    tier_info$n_cells,
    ")"
  ),
  c("High", "Moderate", "Low")
)

heatmap_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.background = element_rect(fill = "#FFFFFF", colour = NA),
    panel.border = element_rect(
      fill = NA, color = "#222222", linewidth = 0.50 / 2.845
    ),
    panel.grid = element_blank(),
    panel.spacing.y = unit(0.18, "cm"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 55, hjust = 1, vjust = 1, size = 8
    ),
    axis.text.y = element_text(size = 7),
    axis.ticks = element_line(color = "#333333", linewidth = 0.25),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5, color = "#444444"),
    plot.caption = element_text(size = 7, color = "#555555", hjust = 0),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(3.0, "cm"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 0, face = "bold", size = 8.5, hjust = 1
    ),
    plot.margin = margin(8, 8, 8, 8)
  )

prob_colors <- c(
  "#FFFFFF",
  "#FFF2B2",
  "#F28E2B",
  "#D64F67",
  "#A22A88",
  "#661B9A",
  "#0D168B"
)
prob_color_values <- scales::rescale(c(
  0,
  0.08,
  0.25,
  0.50,
  0.65,
  0.80,
  1.00
))
confidence_colors <- c(
  High = "#42B540FF",
  Moderate = "#00468BFF",
  Low = "#ED0000FF"
)

p_heatmap <- ggplot(
  long,
  aes(x = Subtype, y = tier_rank, fill = Bootstrap_probability_plot)
) +
  geom_tile(width = 1, height = 1) +
  geom_rect(
    data = tier_info,
    aes(
      xmin = 8.38, xmax = 8.62,
      ymin = 0.5, ymax = n_cells + 0.5
    ),
    inherit.aes = FALSE,
    fill = "#FFFFFF", color = NA
  ) +
  geom_segment(
    data = tier_info,
    aes(
      x = "Confidence", xend = "Confidence",
      y = 0.5, yend = n_cells + 0.5
    ),
    inherit.aes = FALSE,
    color = "#D9D9D9", linewidth = 3.4
  ) +
  geom_segment(
    data = tier_info,
    aes(
      x = "Confidence", xend = "Confidence",
      y = 0.5, yend = n_cells + 0.5,
      color = confidence_class
    ),
    inherit.aes = FALSE,
    linewidth = 2.6
  ) +
  geom_vline(
    xintercept = c(seq(1.5, 7.5, by = 1), seq(9.5, 15.5, by = 1)),
    color = "#222222", linewidth = 0.50 / 2.845,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = c(8.38, 8.62),
    color = "#222222", linewidth = 0.50 / 2.845,
    linetype = "solid"
  ) +
  scale_fill_gradientn(
    colours = prob_colors,
    values = prob_color_values,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = label_number(accuracy = 0.01),
    oob = squish,
    na.value = "#FFFFFF",
    name = "Bootstrap\nprobability"
  ) +
  scale_color_manual(
    values = confidence_colors,
    drop = FALSE,
    name = "Cell confidence"
  ) +
  scale_x_discrete(
    limits = c(subtypes, "Confidence"),
    breaks = subtypes,
    expand = c(0, 0)
  ) +
  scale_y_reverse(
    breaks = function(x) {
      upper <- floor(max(x, na.rm = TRUE))
      unique(c(1, seq(50, upper, by = 50), upper))
    },
    expand = c(0, 0)
  ) +
  facet_grid(
    rows = vars(confidence_class),
    scales = "free_y",
    space = "free_y",
    switch = "y",
    labeller = as_labeller(tier_labels)
  ) +
  labs(
    title = "MSN Top750 Pearson bootstrap probability heatmap by confidence tier",
    subtitle = paste0(
      "588 stable MSN Patch-seq cells shown in three separate tiers: ",
      "High (n=339), Moderate (n=175), Low (n=74)\n",
      "Within each tier, cells are ordered by bootstrap top1 subtype and probability; 500 replicates"
    ),
    x = "Bootstrap-assigned subtype",
    y = "Cells within confidence tier (ordered rank)",
    caption = paste0(
      "A blank gap with 0.5-pt solid flanking borders separates D1 and D2; 0.5-pt dashed vertical lines separate adjacent subtypes within D1 and D2; ",
      "0.5-pt solid frames delimit the three confidence tiers. ",
      "Each tile is the fraction of bootstrap replicates assigning the cell to that subtype.\n",
      "Probability 0 is white. Confidence annotation uses Lancet colors: High (green), Moderate (blue), Low (red)."
    )
  ) +
  guides(
    fill = guide_colorbar(order = 1),
    color = guide_legend(order = 2, override.aes = list(linewidth = 4))
  ) +
  heatmap_theme

figure_png <- file.path(
  out_dir, "Figure_MSN_top750_PC16_Pearson_bootstrap_probability_heatmap_three_tiers_subtype_dashed_gap_bordered.png"
)
figure_pdf <- file.path(
  out_dir, "Figure_MSN_top750_PC16_Pearson_bootstrap_probability_heatmap_three_tiers_subtype_dashed_gap_bordered.pdf"
)
ggsave(
  figure_png, p_heatmap, width = 10.0, height = 11.5,
  units = "in", dpi = 400, bg = "white"
)
ggsave(
  figure_pdf, p_heatmap, width = 10.0, height = 11.5,
  units = "in", device = cairo_pdf, bg = "white"
)

fwrite(
  votes,
  file.path(out_dir, "MSN_top750_PC16_bootstrap_probability_heatmap_cell_order.csv")
)
fwrite(
  long,
  file.path(out_dir, "MSN_top750_PC16_bootstrap_probability_heatmap_long.csv")
)
fwrite(
  group_info,
  file.path(out_dir, "MSN_top750_PC16_bootstrap_probability_heatmap_group_boundaries.csv")
)

summary <- votes[, .(
  n_cells = .N,
  top1_fraction_mean = mean(bootstrap_top1_fraction),
  top1_fraction_median = median(bootstrap_top1_fraction),
  top1_fraction_q25 = quantile(bootstrap_top1_fraction, 0.25),
  top1_fraction_q75 = quantile(bootstrap_top1_fraction, 0.75),
  margin_median = median(bootstrap_margin),
  entropy_median = median(bootstrap_entropy)
), by = bootstrap_top1]
summary[, subtype_order := match(bootstrap_top1, subtypes)]
setorder(summary, subtype_order)
summary[, subtype_order := NULL]
fwrite(
  summary,
  file.path(out_dir, "MSN_top750_PC16_bootstrap_probability_heatmap_summary.csv")
)

confidence_summary <- votes[, .(
  n_cells = .N,
  fraction = .N / nrow(votes),
  rpca_pearson_concordant = sum(RPCA_agrees_bootstrap_top1),
  rpca_pearson_concordance_rate = mean(RPCA_agrees_bootstrap_top1)
), by = confidence_class]
setorder(confidence_summary, confidence_class)
fwrite(
  confidence_summary,
  file.path(out_dir, "MSN_top750_PC16_bootstrap_confidence_class_summary.csv")
)

readme <- c(
  "MSN Top750 Pearson 500-bootstrap probability heatmap",
  "",
  "Rows are the 588 stable MSN Patch-seq cells; columns are the 16 MSN subtypes.",
  "Cells are drawn in three separate vertically stacked confidence tiers: High, Moderate and Low.",
  "Within each tier, rows are ordered independently by bootstrap top1 subtype and descending top1 probability.",
  "The 34 stable IN cells and all non-MSN/QC-risk cells are excluded.",
  "Each tile is the fraction of 500 bootstrap replicates assigning a cell to a subtype.",
  "Rows are ordered by bootstrap top1 subtype (D1_1-D1_8, D2_1-D2_8), then descending top1 fraction.",
  "A blank vertical gap separates D1 and D2, with a 0.5-pt solid border on each side of the gap.",
  "Adjacent subtypes within D1 and D2 are separated by 0.5-pt vertical dashed lines.",
  "Each confidence tier is enclosed by a solid 0.5-pt frame, with a small gap between tiers.",
  "The right-hand strip uses Lancet colors: High=#42B540FF, Moderate=#00468BFF, Low=#ED0000FF.",
  "High: top1 probability >=0.80 and margin >=0.50; Moderate: probability >=0.50 and margin >=0.20 but not High; Low: all remaining cells.",
  "White is included in the continuous probability gradient: white at 0, pale yellow to orange at low values, red-purple at intermediate values, and dark blue at 1.",
  "Bootstrap probabilities are descriptive stability/support metrics and are not RPCA prediction scores.",
  "Source: MSN_top750_PC16_gene_bootstrap_500_vote_fractions.csv."
)
writeLines(
  readme,
  file.path(out_dir, "README_MSN_top750_PC16_bootstrap_probability_heatmap.txt")
)
writeLines("SUCCESS", file.path(out_dir, "SUCCESS"))

message("Output: ", out_dir)
message("Cells: ", nrow(votes), "; long rows: ", nrow(long))
