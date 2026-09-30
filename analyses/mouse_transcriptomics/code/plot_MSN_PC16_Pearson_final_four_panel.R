options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(circlize)
  library(ggplotify)
})

root <- "C:/Users/53461/Documents/Codex/2026-07-13/zh/outputs/05_branch_specific_subtype_mapping"
package_dir <- file.path(root, "60_MSN_PC16_tSNE_complete_package")
bootstrap_dir <- file.path(root, "50_MSN_top750_PC16_gene_bootstrap_500")
out_dir <- file.path(root, "64_MSN_PC16_Pearson_final_four_panel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subtypes <- c(paste0("D1_", 1:8), paste0("D2_", 1:8))
colors <- fread(file.path(package_dir, "labels", "MSN_fixed_subtype_color_mapping.csv"))
subtype_pal <- setNames(colors$color, colors$Subtype)[subtypes]

coords <- fread(file.path(package_dir, "results", "MSN_PC16_tSNE_4x4_all_coordinates.csv"))
xy <- coords[perplexity == 15 & learning_rate == 100]
bootstrap <- fread(file.path(
  bootstrap_dir,
  "MSN_top750_PC16_gene_bootstrap_500_per_cell.csv"
))

pair <- readRDS(file.path(package_dir, "inputs_prepared", "MSN_top750_dispersion_z_pca_pair.rds"))
if (is.null(pair$ref$Subtype)) pair$ref$Subtype <- pair$ref$orig.ident
ref_labels <- data.table(
  cell_id = colnames(pair$ref),
  Reference_Subtype = as.character(pair$ref$Subtype)
)

ref <- merge(
  xy[source == "Reference"],
  ref_labels,
  by = "cell_id", all.x = TRUE, sort = FALSE
)
qry <- merge(
  xy[source == "Query"],
  bootstrap[, .(
    cell_id,
    RPCA_Subtype = RPCA_top1,
    Pearson_Top750_Subtype = bootstrap_top1,
    Pearson_Top750_Prob = bootstrap_top1_fraction,
    Pearson_Top750_second_Subtype = bootstrap_top2,
    Pearson_Top750_second_Prob = bootstrap_top2_fraction,
    Pearson_Top750_vote_margin = bootstrap_margin,
    RPCA_Pearson_agreement = RPCA_agrees_bootstrap_top1
  )],
  by = "cell_id", all.x = TRUE, sort = FALSE
)

stopifnot(
  nrow(ref) == 20029L,
  nrow(qry) == 588L,
  uniqueN(qry$cell_id) == 588L,
  !anyNA(ref$Reference_Subtype),
  !anyNA(qry$Pearson_Top750_Subtype),
  !anyNA(qry$Pearson_Top750_Prob),
  all(qry$Reference_Subtype %in% c(NA_character_)),
  all(qry$Pearson_Top750_Subtype %in% subtypes),
  sum(qry$RPCA_Pearson_agreement) == 422L
)

ref[, Reference_Subtype := factor(Reference_Subtype, levels = subtypes)]
qry[, Pearson_Top750_Subtype := factor(Pearson_Top750_Subtype, levels = subtypes)]
qry_discordant <- qry[RPCA_Pearson_agreement == FALSE]
qry_consensus <- qry[RPCA_Pearson_agreement == TRUE]

x_lim <- range(xy$tSNE_1, finite = TRUE)
y_lim <- range(xy$tSNE_2, finite = TRUE)

prob_core_colors <- viridisLite::magma(256, direction = 1)
prob_scale_colors <- c(prob_core_colors[1], prob_core_colors, prob_core_colors[256])
prob_scale_values <- c(0, seq(0.25, 0.75, length.out = 256), 1)

common_theme <- theme_void(base_size = 10) +
  theme(
    plot.title = element_text(
      size = 12, face = "bold", hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "none"
  )

p_ref <- ggplot() +
  geom_point(
    data = ref,
    aes(tSNE_1, tSNE_2, fill = Reference_Subtype),
    shape = 21, color = "#111111", stroke = 0, size = 0.42, alpha = 0.82
  ) +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE) +
  coord_equal(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  labs(title = "A  Reference subtypes") +
  common_theme

p_query <- ggplot() +
  geom_point(
    data = ref,
    aes(tSNE_1, tSNE_2),
    color = "#D0D0D0", size = 0.17, alpha = 0.48
  ) +
  geom_point(
    data = qry,
    aes(tSNE_1, tSNE_2, color = Pearson_Top750_Subtype),
    shape = 16, size = 1.45, alpha = 0.98
  ) +
  scale_color_manual(values = subtype_pal, limits = subtypes, drop = FALSE) +
  coord_equal(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  labs(title = "B  Pearson Top750 subtype") +
  common_theme

p_probability <- ggplot() +
  geom_point(
    data = ref,
    aes(tSNE_1, tSNE_2),
    color = "#D0D0D0", size = 0.17, alpha = 0.48
  ) +
  geom_point(
    data = qry_discordant,
    aes(tSNE_1, tSNE_2, fill = Pearson_Top750_Prob),
    shape = 21, color = "#111111", stroke = 0, size = 1.45, alpha = 0.96
  ) +
  geom_point(
    data = qry_consensus,
    aes(tSNE_1, tSNE_2, fill = Pearson_Top750_Prob),
    shape = 21, color = "#111111", stroke = 0.38, size = 1.45, alpha = 0.98
  ) +
  scale_fill_gradientn(
    colours = prob_scale_colors,
    values = prob_scale_values,
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    oob = scales::squish
  ) +
  coord_equal(xlim = x_lim, ylim = y_lim, expand = FALSE) +
  labs(title = "C  Pearson bootstrap top1 Prob.") +
  common_theme

links <- qry[, .N, by = .(
  Pearson = paste0("Pearson_", as.character(Pearson_Top750_Subtype)),
  RPCA = paste0("RPCA_", RPCA_Subtype),
  Pearson_Subtype = as.character(Pearson_Top750_Subtype)
)]
setorder(links, Pearson, RPCA)

pearson_d1 <- paste0("Pearson_D1_", 1:8)
pearson_d2 <- paste0("Pearson_D2_", 1:8)
rpca_d1 <- paste0("RPCA_D1_", 8:1)
rpca_d2 <- paste0("RPCA_D2_", 8:1)
sector_order <- c(pearson_d1, pearson_d2, rpca_d2, rpca_d1)
sector_subtype <- sub("^(Pearson|RPCA)_", "", sector_order)
grid_colors <- setNames(subtype_pal[sector_subtype], sector_order)
ribbon_colors <- adjustcolor(subtype_pal[links$Pearson_Subtype], alpha.f = 0.68)
major_colors <- c(D1 = "#8B1E24", D2 = "#173F73")
gap_after <- rep(1.0, length(sector_order))
gap_after[c(8, 16, 24, 32)] <- c(4.5, 13, 4.5, 13)

draw_chord <- function() {
  circos.clear()
  circos.par(
    start.degree = 92,
    clock.wise = TRUE,
    gap.after = gap_after,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.002, 0.002),
    points.overflow.warning = FALSE
  )
  chordDiagram(
    x = links[, .(Pearson, RPCA, N)],
    order = sector_order,
    grid.col = grid_colors,
    col = ribbon_colors,
    transparency = 0,
    link.border = adjustcolor("#FFFFFF", alpha.f = 0.20),
    link.lwd = 0.15,
    directional = 1,
    direction.type = "diffHeight",
    diffHeight = mm_h(1.2),
    link.arr.type = "big.arrow",
    annotationTrack = "grid",
    preAllocateTracks = list(
      list(track.height = 0.075),
      list(track.height = 0.080)
    )
  )
  circos.trackPlotRegion(
    track.index = 2,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector <- CELL_META$sector.index
      label <- sub("^(Pearson|RPCA)_", "", sector)
      label_facing <- if (sector == "RPCA_D2_8") "reverse.clockwise" else "clockwise"
      label_nice_facing <- sector != "RPCA_D2_8"
      circos.text(
        CELL_META$xcenter, CELL_META$ycenter, label,
        facing = label_facing, niceFacing = label_nice_facing,
        cex = 0.44, font = 2, col = "#222222"
      )
    }
  )
  highlight.sector(
    pearson_d1, track.index = 1,
    col = major_colors["D1"], border = major_colors["D1"], lwd = 0.8,
    text = "Pearson · D1", cex = 0.62, font = 2,
    text.col = "#FFFFFF", facing = "bending.inside"
  )
  highlight.sector(
    pearson_d2, track.index = 1,
    col = major_colors["D2"], border = major_colors["D2"], lwd = 0.8,
    text = "Pearson · D2", cex = 0.62, font = 2,
    text.col = "#FFFFFF", facing = "bending.inside"
  )
  highlight.sector(
    rpca_d2, track.index = 1,
    col = major_colors["D2"], border = major_colors["D2"], lwd = 0.8,
    text = "RPCA · D2", cex = 0.62, font = 2,
    text.col = "#FFFFFF", facing = "bending.inside"
  )
  highlight.sector(
    rpca_d1, track.index = 1,
    col = major_colors["D1"], border = major_colors["D1"], lwd = 0.8,
    text = "RPCA · D1", cex = 0.62, font = 2,
    text.col = "#FFFFFF", facing = "bending.inside"
  )
  title("D  Pearson → RPCA concordance", line = 0.5, cex.main = 1.0, font.main = 2)
  circos.clear()
}

p_chord <- ggplotify::as.ggplot(~{
  par(mar = c(0.3, 0.3, 1.8, 0.3))
  draw_chord()
}) + theme(aspect.ratio = 1)

celltype_key <- data.table(
  Subtype = factor(subtypes, levels = subtypes),
  major = rep(c("D1", "D2"), each = 8),
  subtype_index = rep(1:8, 2)
)
celltype_key[, `:=`(
  x_key = ifelse(major == "D1", 0.00, 1.20),
  x_text = ifelse(major == "D1", 0.18, 1.38),
  y = 9 - subtype_index
)]

p_figure_subtype_legend <- ggplot(celltype_key) +
  geom_point(
    aes(x_key, y, fill = Subtype),
    shape = 22, color = "#00000000", stroke = 0, size = 3.5
  ) +
  geom_text(
    aes(x_text, y, label = as.character(Subtype)),
    hjust = 0, size = 2.50, color = "#222222"
  ) +
  annotate("text", x = 0.38, y = 8.78, label = "D1", fontface = "bold", size = 2.9) +
  annotate("text", x = 1.58, y = 8.78, label = "D2", fontface = "bold", size = 2.9) +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE, guide = "none") +
  scale_x_continuous(limits = c(-0.10, 2.52), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.45, 9.15), expand = c(0, 0)) +
  labs(title = "Cell subtype colors") +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5, margin = margin(b = 1)),
    plot.margin = margin(1, 3, 0, 3),
    plot.background = element_rect(fill = "white", color = NA)
  )

agreement_key <- data.table(
  x = 0,
  y = c(1, 0),
  agreement = c(TRUE, FALSE),
  label = c(
    "Consensus: RPCA = Pearson (n=422)",
    "Discordant: RPCA != Pearson (n=166)"
  )
)
p_figure_agreement_legend <- ggplot() +
  geom_point(
    data = agreement_key[agreement == TRUE], aes(x, y),
    shape = 21, fill = "#BDBDBD", color = "#111111", stroke = 0.55, size = 3.0
  ) +
  geom_point(
    data = agreement_key[agreement == FALSE], aes(x, y),
    shape = 21, fill = "#BDBDBD", color = "#00000000", stroke = 0, size = 3.0
  ) +
  geom_text(
    data = agreement_key, aes(x, y, label = label),
    hjust = 0, nudge_x = 0.18, size = 2.25, color = "#222222"
  ) +
  scale_x_continuous(limits = c(-0.08, 2.82), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.55, 1.55), expand = c(0, 0)) +
  labs(title = "RPCA-Pearson agreement (C)") +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(size = 8.2, face = "bold", hjust = 0.5, margin = margin(b = 1)),
    plot.margin = margin(1, 3, 1, 3),
    plot.background = element_rect(fill = "white", color = NA)
  )

prob_density_fit <- density(
  qry$Pearson_Top750_Prob,
  from = 0, to = 1, n = 512, adjust = 0.15, cut = 0
)
prob_density_data <- data.table(
  probability = prob_density_fit$x,
  cell_count = prob_density_fit$y * nrow(qry) * 0.025
)
prob_reference_lines <- data.table(
  cell_count = c(0, 50, 100, 150),
  label = c("0", "50", "100", "150"),
  vjust = c(-0.35, -0.35, -0.35, 1.10)
)

p_figure_probability_mountain <- ggplot(
  prob_density_data,
  aes(probability, cell_count)
) +
  geom_area(fill = "#BDBDBD", alpha = 0.48, color = NA) +
  geom_hline(
    data = prob_reference_lines,
    aes(yintercept = cell_count),
    color = "#8A8A8A", linewidth = 0.28
  ) +
  geom_line(color = "#333333", linewidth = 0.48, lineend = "round") +
  geom_text(
    data = prob_reference_lines,
    aes(x = 0.012, y = cell_count, label = label, vjust = vjust),
    inherit.aes = FALSE, hjust = 0, size = 1.85, color = "#555555"
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 155), expand = expansion(mult = 0)) +
  labs(
    title = "Pearson bootstrap top1 probability",
    subtitle = "Smoothed query-cell count per 0.025 bin (n=588)",
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.title = element_text(size = 8.2, face = "bold", hjust = 0.5, margin = margin(b = 0)),
    plot.subtitle = element_text(size = 6.8, color = "#555555", hjust = 0.5, margin = margin(b = 0)),
    plot.margin = margin(1, 4, 0, 4),
    plot.background = element_rect(fill = "white", color = NA)
  )

prob_bar_figure_data <- data.table(
  probability = seq(0.0005, 0.9995, length.out = 1000),
  y = 1
)
p_figure_probability_bar <- ggplot(
  prob_bar_figure_data,
  aes(probability, y, fill = probability)
) +
  geom_tile(width = 0.001, height = 0.18) +
  annotate(
    "rect", xmin = 0, xmax = 1, ymin = 0.91, ymax = 1.09,
    fill = NA, color = "#444444", linewidth = 0.25
  ) +
  scale_fill_gradientn(
    colours = prob_scale_colors, values = prob_scale_values,
    limits = c(0, 1), oob = scales::squish, guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = sprintf("%.2f", c(0, 0.25, 0.50, 0.75, 1)),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(limits = c(0.74, 1.26), expand = expansion(mult = 0)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(color = "#444444", linewidth = 0.25),
    axis.text.x = element_text(size = 7.2, color = "#222222"),
    plot.margin = margin(0, 4, 1, 4),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_figure_probability_legend <-
  p_figure_probability_mountain / p_figure_probability_bar +
  plot_layout(heights = c(0.62, 0.38))

p_probability_mountain_no_text <- ggplot(
  prob_density_data,
  aes(probability, cell_count)
) +
  geom_area(fill = "#BDBDBD", alpha = 0.48, color = NA) +
  geom_hline(
    data = prob_reference_lines,
    aes(yintercept = cell_count),
    color = "#8A8A8A", linewidth = 0.28
  ) +
  geom_line(color = "#333333", linewidth = 0.48, lineend = "round") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0, 155), expand = expansion(mult = 0)) +
  theme_void() +
  theme(plot.margin = margin(2, 4, 0, 4))

p_probability_bar_no_text <- ggplot(
  prob_bar_figure_data,
  aes(probability, y, fill = probability)
) +
  geom_tile(width = 0.001, height = 0.18) +
  annotate(
    "rect", xmin = 0, xmax = 1, ymin = 0.91, ymax = 1.09,
    fill = NA, color = "#444444", linewidth = 0.25
  ) +
  scale_fill_gradientn(
    colours = prob_scale_colors, values = prob_scale_values,
    limits = c(0, 1), oob = scales::squish, guide = "none"
  ) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.74, 1.26), expand = expansion(mult = 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 4, 2, 4))

probability_only_with_text <-
  p_figure_probability_mountain / p_figure_probability_bar +
  plot_layout(heights = c(0.68, 0.32))

probability_only_no_text <-
  p_probability_mountain_no_text / p_probability_bar_no_text +
  plot_layout(heights = c(0.72, 0.28))

probability_text_png <- file.path(out_dir, "Legend_Pearson_top1_probability_distribution_with_text.png")
probability_text_pdf <- file.path(out_dir, "Legend_Pearson_top1_probability_distribution_with_text.pdf")
probability_no_text_png <- file.path(out_dir, "Legend_Pearson_top1_probability_distribution_no_text.png")
probability_no_text_pdf <- file.path(out_dir, "Legend_Pearson_top1_probability_distribution_no_text.pdf")

ggsave(probability_text_png, probability_only_with_text, width = 4.8, height = 2.0, dpi = 400, bg = "white")
ggsave(probability_text_pdf, probability_only_with_text, width = 4.8, height = 2.0, device = cairo_pdf, bg = "white")
ggsave(probability_no_text_png, probability_only_no_text, width = 4.8, height = 2.0, dpi = 400, bg = "white")
ggsave(probability_no_text_pdf, probability_only_no_text, width = 4.8, height = 2.0, device = cairo_pdf, bg = "white")

legend_separator <- ggplot() +
  theme_void() +
  theme(plot.background = element_rect(fill = "#111111", color = NA))

right_legend_column <-
  p_figure_subtype_legend /
  legend_separator /
  p_figure_agreement_legend /
  legend_separator /
  p_figure_probability_legend /
  plot_spacer() +
  plot_layout(heights = c(0.42, 0.012, 0.13, 0.012, 0.21, 0.216))

scale_bar_length <- 10
scale_bar_x <- x_lim[1] + 0.055 * diff(x_lim)
scale_bar_y <- y_lim[1] + 0.055 * diff(y_lim)
scale_bar_tick_half <- 0.65

add_tsne_scale_bar <- function(p) {
  p +
    annotate(
      "rect",
      xmin = scale_bar_x - 1.2, xmax = scale_bar_x + scale_bar_length + 1.2,
      ymin = scale_bar_y - 1.4, ymax = scale_bar_y + 3.0,
      fill = adjustcolor("#FFFFFF", alpha.f = 0.82), color = NA
    ) +
    annotate(
      "segment",
      x = scale_bar_x, xend = scale_bar_x + scale_bar_length,
      y = scale_bar_y, yend = scale_bar_y,
      color = "#111111", linewidth = 0.48, lineend = "butt"
    ) +
    annotate(
      "segment",
      x = scale_bar_x, xend = scale_bar_x,
      y = scale_bar_y - scale_bar_tick_half, yend = scale_bar_y + scale_bar_tick_half,
      color = "#111111", linewidth = 0.48, lineend = "butt"
    ) +
    annotate(
      "segment",
      x = scale_bar_x + scale_bar_length, xend = scale_bar_x + scale_bar_length,
      y = scale_bar_y - scale_bar_tick_half, yend = scale_bar_y + scale_bar_tick_half,
      color = "#111111", linewidth = 0.48, lineend = "butt"
    ) +
    annotate(
      "text",
      x = scale_bar_x + scale_bar_length / 2,
      y = scale_bar_y + 1.55,
      label = "10 t-SNE units",
      size = 2.45, color = "#111111"
    )
}

p_ref_scaled <- add_tsne_scale_bar(p_ref)
p_query_scaled <- add_tsne_scale_bar(p_query)
p_probability_scaled <- add_tsne_scale_bar(p_probability)

main_abc_plot <- wrap_plots(
  p_ref_scaled, p_query_scaled, p_probability_scaled,
  ncol = 3, widths = c(1, 1, 1)
) +
  plot_annotation(
    title = "MSN PC16 t-SNE Pearson mapping",
    subtitle = "Perplexity 15 | eta 100 | Pearson Top750 / 721-gene | 500 bootstrap replicates",
    theme = theme(
      plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, color = "#444444", hjust = 0.5)
    )
  )

main_abc_png <- file.path(out_dir, "Figure_MSN_PC16_Pearson_ABC_with_tSNE_scale_bars.png")
main_abc_pdf <- file.path(out_dir, "Figure_MSN_PC16_Pearson_ABC_with_tSNE_scale_bars.pdf")
ggsave(main_abc_png, main_abc_plot, width = 15.6, height = 6.25, dpi = 400, bg = "white")
ggsave(main_abc_pdf, main_abc_plot, width = 15.6, height = 6.25, device = cairo_pdf, bg = "white")
ggsave(
  file.path(out_dir, "Figure_MSN_PC16_Pearson_ABC_with_tSNE_scale_bars_B_Pearson_no_border.png"),
  main_abc_plot, width = 15.6, height = 6.25, dpi = 400, bg = "white"
)
ggsave(
  file.path(out_dir, "Figure_MSN_PC16_Pearson_ABC_with_tSNE_scale_bars_B_Pearson_no_border.pdf"),
  main_abc_plot, width = 15.6, height = 6.25, device = cairo_pdf, bg = "white"
)

standalone_legend <-
  p_figure_subtype_legend /
  legend_separator /
  p_figure_agreement_legend /
  legend_separator /
  p_figure_probability_legend +
  plot_layout(heights = c(0.529, 0.017, 0.182, 0.017, 0.255)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = "#111111", linewidth = 0.45),
      plot.margin = margin(2, 2, 2, 2)
    )
  )

standalone_legend_png <- file.path(out_dir, "Legend_MSN_PC16_Pearson_three_level_compact_tighter_labels.png")
standalone_legend_pdf <- file.path(out_dir, "Legend_MSN_PC16_Pearson_three_level_compact_tighter_labels.pdf")
ggsave(standalone_legend_png, standalone_legend, width = 2.8, height = 4.35, dpi = 400, bg = "white")
ggsave(standalone_legend_pdf, standalone_legend, width = 2.8, height = 4.35, device = cairo_pdf, bg = "white")

final_three_panel <- wrap_plots(
  p_ref, p_query, p_probability, right_legend_column,
  ncol = 4, widths = c(1, 1, 1, 0.48)
) +
  plot_annotation(
    title = "MSN PC16 t-SNE Pearson mapping",
    subtitle = "Perplexity 15 | eta 100 | Pearson Top750 / 721-gene | 500 bootstrap replicates",
    theme = theme(
      plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 9.5, color = "#444444", hjust = 0.5)
    )
  )

three_panel_png <- file.path(out_dir, "Figure_MSN_PC16_Pearson_final_three_panel_right_legend_block_compact.png")
three_panel_pdf <- file.path(out_dir, "Figure_MSN_PC16_Pearson_final_three_panel_right_legend_block_compact.pdf")
ggsave(three_panel_png, final_three_panel, width = 17.26, height = 6.25, dpi = 400, bg = "white")
ggsave(three_panel_pdf, final_three_panel, width = 17.26, height = 6.25, device = cairo_pdf, bg = "white")
ggsave(
  file.path(out_dir, "Figure_MSN_PC16_Pearson_final_three_panel_B_Pearson_no_border.png"),
  final_three_panel, width = 17.26, height = 6.25, dpi = 400, bg = "white"
)
ggsave(
  file.path(out_dir, "Figure_MSN_PC16_Pearson_final_three_panel_B_Pearson_no_border.pdf"),
  final_three_panel, width = 17.26, height = 6.25, device = cairo_pdf, bg = "white"
)

png(
  file.path(out_dir, "Panel_D_Pearson_primary_RPCA_chord.png"),
  width = 6, height = 6, units = "in", res = 600, bg = "white"
)
par(mar = c(0.3, 0.3, 1.8, 0.3))
draw_chord()
dev.off()
cairo_pdf(
  file.path(out_dir, "Panel_D_Pearson_primary_RPCA_chord.pdf"),
  width = 6, height = 6, bg = "white"
)
par(mar = c(0.3, 0.3, 1.8, 0.3))
draw_chord()
dev.off()

p_subtype_legend <- ggplot(
  data.table(Pearson_Top750_Subtype = factor(subtypes, levels = subtypes)),
  aes(Pearson_Top750_Subtype, 1, fill = Pearson_Top750_Subtype)
) +
  geom_point(shape = 22, color = "#00000000", alpha = 0, show.legend = TRUE) +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE) +
  labs(fill = NULL) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9, color = "#222222"),
    legend.key.height = unit(0.48, "cm"),
    legend.key.width = unit(0.48, "cm"),
    legend.spacing.x = unit(0.35, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(fill = guide_legend(
    ncol = 2, byrow = FALSE,
    override.aes = list(shape = 22, color = NA, alpha = 1, stroke = 0, size = 5)
  ))

subtype_legend <- cowplot::get_legend(p_subtype_legend)
cowplot::save_plot(
  file.path(out_dir, "Legend_MSN_Pearson_Top750_subtype_two_column.png"),
  subtype_legend, base_width = 4.1, base_height = 3.1, dpi = 400, bg = "white"
)
ggsave(
  file.path(out_dir, "Legend_MSN_Pearson_Top750_subtype_two_column.pdf"),
  subtype_legend, width = 4.1, height = 3.1, device = cairo_pdf, bg = "white"
)

p_prob_histogram <- ggplot(qry, aes(Pearson_Top750_Prob, fill = after_stat(x))) +
  geom_histogram(binwidth = 0.025, boundary = 0, color = "#333333", linewidth = 0.28) +
  scale_fill_gradientn(
    colours = prob_scale_colors, values = prob_scale_values,
    limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1),
    oob = scales::squish, guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Pearson bootstrap top1 probability", x = NULL, y = "Query cells") +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(5, 12, 0, 12)
  )

prob_bar_data <- data.table(probability = seq(0.0005, 0.9995, length.out = 1000), y = 1)
p_prob_bar <- ggplot(prob_bar_data, aes(probability, y, fill = probability)) +
  geom_tile(width = 0.001, height = 0.55 / 3) +
  annotate(
    "rect", xmin = 0, xmax = 1, ymin = 1 - 0.55 / 6, ymax = 1 + 0.55 / 6,
    fill = NA, color = "#444444", linewidth = 0.25
  ) +
  scale_fill_gradientn(
    colours = prob_scale_colors, values = prob_scale_values,
    limits = c(0, 1), oob = scales::squish, guide = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = sprintf("%.2f", c(0, 0.25, 0.50, 0.75, 1)),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(limits = c(0.6, 1.4), expand = expansion(mult = 0)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.line = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(), axis.ticks.x = element_line(color = "#444444"),
    axis.text.x = element_text(size = 9, color = "#222222"),
    plot.margin = margin(0, 12, 5, 12)
  )

prob_legend <- p_prob_histogram / p_prob_bar + plot_layout(heights = c(2.2, 0.8))
cowplot::save_plot(
  file.path(out_dir, "Legend_Pearson_Top750_bootstrap_top1_probability_magma.png"),
  prob_legend, base_width = 7.2, base_height = 3.0, dpi = 400, bg = "white"
)
ggsave(
  file.path(out_dir, "Legend_Pearson_Top750_bootstrap_top1_probability_magma.pdf"),
  prob_legend, width = 7.2, height = 3.0, device = cairo_pdf, bg = "white"
)

fwrite(qry, file.path(out_dir, "MSN_PC16_per15_eta100_Pearson_Top750_bootstrap_plot_data.csv"))
fwrite(links, file.path(out_dir, "MSN_Pearson_Top750_RPCA_subtype_chord_links.csv"))
fwrite(
  qry[, .(
    n_query = .N,
    n_RPCA_Pearson_agreement = sum(RPCA_Pearson_agreement),
    n_RPCA_Pearson_discordant = sum(!RPCA_Pearson_agreement),
    agreement_fraction = mean(RPCA_Pearson_agreement),
    top1_probability_mean = mean(Pearson_Top750_Prob),
    top1_probability_median = median(Pearson_Top750_Prob),
    top1_probability_q25 = quantile(Pearson_Top750_Prob, 0.25),
    top1_probability_q75 = quantile(Pearson_Top750_Prob, 0.75)
  )],
  file.path(out_dir, "MSN_PC16_Pearson_final_four_panel_summary.csv")
)

readme <- c(
  "MSN PC16 Pearson final three-panel figure with integrated legends",
  "",
  "A: Reference subtype map.",
  "B: Query labels are Pearson Top750 / 721-gene 500-bootstrap modal top1 labels.",
  "C: Query feature is the 500-bootstrap top1 vote probability.",
  "Black border: RPCA and Pearson top1 agree; no border: discordant.",
  "Right legend column: D1 and D2 subtype-color columns, RPCA-Pearson agreement border styles, and probability mapping.",
  "The probability legend includes a continuous kernel-density mountain scaled to smoothed counts per 0.025 probability interval.",
  "t-SNE: PC1-16, perplexity 15, eta 100.",
  "The requested final layout is exported as a single A-B-C row.",
  "Formal RPCA-primary labels were not modified."
)
writeLines(readme, file.path(out_dir, "README_MSN_PC16_Pearson_final_three_panel.txt"))
writeLines("SUCCESS", file.path(out_dir, "SUCCESS"))

message("Output: ", out_dir)
message("Query cells: ", nrow(qry))
message("RPCA-Pearson agreement: ", sum(qry$RPCA_Pearson_agreement))
