options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
})

root <- "C:/Users/53461/Documents/Codex/2026-07-13/zh/outputs/05_branch_specific_subtype_mapping"
package_dir <- file.path(root, "60_MSN_PC16_tSNE_complete_package")
bootstrap_dir <- file.path(root, "50_MSN_top750_PC16_gene_bootstrap_500")
raw_input_file <- file.path(root, "01_MSN_input", "MSN_raw_branch_counts_and_metadata.rds")
out_dir <- file.path(root, "65_MSN_RefQuery_normalized_stacked_violin")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

subtypes <- c(paste0("D1_", 1:8), paste0("D2_", 1:8))
genes <- c(
  "Ppp1r1b", "Drd1", "Drd2", "Resp18", "Onecut2",
  "Tcerg1l", "Penk", "Nts", "Igfbp4", "Itm2a",
  "Ncald", "Gng2", "Calcr", "Tac1", "Lamp5",
  "Mfge8", "Oprk1", "Pde1a", "Stard5", "Cpne4"
)

colors <- fread(file.path(package_dir, "labels", "MSN_fixed_subtype_color_mapping.csv"))
subtype_pal <- setNames(colors$color, colors$Subtype)[subtypes]

pair <- readRDS(raw_input_file)
bootstrap <- fread(file.path(
  bootstrap_dir,
  "MSN_top750_PC16_gene_bootstrap_500_per_cell.csv"
))

query_label_map <- setNames(bootstrap$bootstrap_top1, bootstrap$cell_id)
ref_types <- as.character(pair$ref_meta$Subtype)
query_types <- unname(query_label_map[colnames(pair$query_counts)])

stopifnot(
  all(genes %in% rownames(pair$ref_counts)),
  all(genes %in% rownames(pair$query_counts)),
  all(ref_types %in% subtypes),
  all(query_types %in% subtypes),
  !anyNA(query_types)
)

# Reference and Query are independently library-size normalized from the full
# 23,950-gene raw matrices, followed by log1p transformation.
# For display, each gene is divided by its combined Reference+Query maximum;
# the original log-normalized maximum is printed at the right of that row.
normalize_requested <- function(counts, genes) {
  library_size <- Matrix::colSums(counts)
  scale_factor <- 10000 / pmax(library_size, 1)
  normalized <- counts[genes, , drop = FALSE] %*%
    Matrix::Diagonal(x = scale_factor)
  normalized@x <- log1p(normalized@x)
  as.matrix(normalized)
}
ref_expr <- normalize_requested(pair$ref_counts, genes)
query_expr <- normalize_requested(pair$query_counts, genes)
expr <- cbind(ref_expr, query_expr)
cell_types <- c(ref_types, query_types)
sources <- c(rep("Reference", ncol(ref_expr)), rep("Query", ncol(query_expr)))

gene_max <- apply(expr, 1, max, na.rm = TRUE)
gene_max[!is.finite(gene_max) | gene_max <= 0] <- 1

plot_dt <- data.table(
  Gene = rep(genes, times = ncol(expr)),
  Subtype = rep(cell_types, each = length(genes)),
  Source = rep(sources, each = length(genes)),
  Expression = as.vector(expr)
)
plot_dt[, NormalizedExpression := Expression / gene_max[Gene]]
plot_dt[, `:=`(
  Gene = factor(Gene, levels = genes),
  Subtype = factor(Subtype, levels = subtypes),
  SubtypeIndex = match(as.character(Subtype), subtypes),
  Source = factor(Source, levels = c("Reference", "Query"))
)]
plot_dt[, PairX := SubtypeIndex + fifelse(Source == "Reference", -0.20, 0.20)]

max_dt <- data.table(
  Gene = factor(genes, levels = genes),
  MaxExpression = as.numeric(gene_max[genes]),
  MaxLabel = sprintf("%.2f", as.numeric(gene_max[genes])),
  x = 16.78,
  y = 0.56
)

count_dt <- data.table(
  Subtype = factor(cell_types, levels = subtypes)
)[, .N, by = Subtype]
count_dt[, SubtypeIndex := match(as.character(Subtype), subtypes)]

source_count_dt <- data.table(
  Subtype = factor(cell_types, levels = subtypes),
  Source = factor(sources, levels = c("Reference", "Query"))
)[, .(CellCount = .N), by = .(Source, Subtype)]
source_count_dt <- merge(
  CJ(
    Source = factor(c("Reference", "Query"), levels = c("Reference", "Query")),
    Subtype = factor(subtypes, levels = subtypes),
    unique = TRUE
  ),
  source_count_dt,
  by = c("Source", "Subtype"),
  all.x = TRUE,
  sort = FALSE
)
source_count_dt[is.na(CellCount), CellCount := 0L]
source_count_dt[, SourceTotal := sum(CellCount), by = Source]
source_count_dt[, `:=`(
  Percent = 100 * CellCount / SourceTotal,
  SubtypeIndex = match(as.character(Subtype), subtypes)
)]
source_count_dt[, PairX := SubtypeIndex + fifelse(Source == "Reference", -0.20, 0.20)]

strip_dt <- data.table(
  Subtype = factor(subtypes, levels = subtypes),
  SubtypeIndex = 1:16,
  y = 1
)

p_top <- ggplot(strip_dt, aes(SubtypeIndex, y, fill = Subtype)) +
  geom_tile(width = 0.98, height = 0.68) +
  annotate("text", x = 16.78, y = 1, label = "Max", size = 1.45, color = "#555555") +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE, guide = "none") +
  scale_x_continuous(limits = c(0.45, 17.15), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.55, 1.45), expand = expansion(mult = 0)) +
  theme_void() +
  theme(plot.margin = margin(1, 2, 0, 40, unit = "pt"))

p_violin <- ggplot(
  plot_dt,
  aes(
    x = SubtypeIndex,
    y = NormalizedExpression,
    group = Subtype,
    fill = Subtype
  )
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.15) +
  geom_violin(
    width = 0.88,
    scale = "width",
    trim = TRUE,
    adjust = 0.72,
    color = NA,
    linewidth = 0,
    alpha = 0.98
  ) +
  geom_text(
    data = max_dt,
    aes(x = x, y = y, label = MaxLabel),
    inherit.aes = FALSE,
    hjust = 0.5,
    size = 1.35,
    color = "#666666"
  ) +
  facet_grid(rows = vars(Gene), switch = "y") +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE, guide = "none") +
  scale_x_continuous(
    limits = c(0.45, 17.15),
    breaks = 1:16,
    labels = subtypes,
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(limits = c(0, 1.12), expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0, "pt"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      size = 4.2, color = "#333333",
      angle = 90, vjust = 0.5, hjust = 1,
      margin = margin(t = 1)
    ),
    axis.ticks.x = element_line(color = "#555555", linewidth = 0.25),
    axis.ticks.length.x = unit(1.5, "pt"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 0, hjust = 1,
      size = 5.0, color = "#222222",
      margin = margin(r = 3)
    ),
    plot.margin = margin(0, 2, 0, 2, unit = "pt")
  )

half_violin_dt <- plot_dt[, {
  values <- NormalizedExpression[is.finite(NormalizedExpression)]
  y_grid <- seq(0, 1, length.out = 81)
  if (length(values) < 2L || diff(range(values)) < 1e-8) {
    center_y <- if (length(values)) values[1] else 0
    density_y <- dnorm(y_grid, mean = center_y, sd = 0.025)
  } else {
    density_y <- density(
      values,
      from = 0, to = 1, n = length(y_grid),
      adjust = 0.72, cut = 0
    )$y
  }
  half_width <- if (max(density_y) > 0) 0.41 * density_y / max(density_y) else rep(0, length(y_grid))
  direction <- if (as.character(Source[1]) == "Reference") -1 else 1
  center_x <- SubtypeIndex[1]
  outer_x <- center_x + direction * half_width
  .(
    x = c(outer_x, center_x, center_x),
    y = c(y_grid, 1, 0),
    Vertex = seq_len(length(y_grid) + 2L)
  )
}, by = .(Gene, Subtype, SubtypeIndex, Source)]
half_violin_dt[, Source := factor(Source, levels = c("Reference", "Query"))]

p_top_split <- ggplot(strip_dt, aes(SubtypeIndex, y, fill = Subtype)) +
  geom_tile(width = 0.98, height = 0.68) +
  annotate(
    "text", x = 8.5, y = 1.58,
    label = "Reference ← | → Query",
    size = 1.55, fontface = "bold", color = "#444444"
  ) +
  annotate("text", x = 16.78, y = 1, label = "Max", size = 1.45, color = "#555555") +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE, guide = "none") +
  scale_x_continuous(limits = c(0.45, 17.15), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.55, 1.75), expand = expansion(mult = 0)) +
  theme_void() +
  theme(plot.margin = margin(1, 2, 0, 40, unit = "pt"))

p_violin_split <- ggplot(
  half_violin_dt,
  aes(
    x = x, y = y,
    group = interaction(Subtype, Source, drop = TRUE),
    fill = Subtype,
    alpha = Source
  )
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.15) +
  geom_polygon(color = NA, linewidth = 0) +
  geom_text(
    data = max_dt,
    aes(x = x, y = y, label = MaxLabel),
    inherit.aes = FALSE,
    hjust = 0.5,
    size = 1.35,
    color = "#666666"
  ) +
  facet_grid(rows = vars(Gene), switch = "y") +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE, guide = "none") +
  scale_alpha_manual(values = c(Reference = 0.98, Query = 0.56), guide = "none") +
  scale_x_continuous(
    limits = c(0.45, 17.15),
    breaks = 1:16,
    labels = subtypes,
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(limits = c(0, 1.12), expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0, "pt"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      size = 4.2, color = "#333333",
      angle = 90, vjust = 0.5, hjust = 1,
      margin = margin(t = 1)
    ),
    axis.ticks.x = element_line(color = "#555555", linewidth = 0.25),
    axis.ticks.length.x = unit(1.5, "pt"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 0, hjust = 1,
      size = 5.0, color = "#222222",
      margin = margin(r = 3)
    ),
    plot.margin = margin(0, 2, 0, 2, unit = "pt")
  )

p_top_paired <- ggplot(strip_dt, aes(SubtypeIndex, y, fill = Subtype)) +
  geom_tile(width = 0.98, height = 0.40) +
  annotate("text", x = 16.88, y = 1, label = "Max", size = 1.45, color = "#555555") +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE, guide = "none") +
  scale_x_continuous(limits = c(0.45, 17.25), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.77, 1.23), expand = expansion(mult = 0)) +
  theme_void() +
  theme(plot.margin = margin(1, 2, 0, 40, unit = "pt"))

p_violin_paired <- ggplot(
  plot_dt,
  aes(
    x = PairX,
    y = NormalizedExpression,
    group = interaction(Subtype, Source, drop = TRUE),
    fill = Subtype
  )
) +
  geom_hline(yintercept = 0, color = "#D0D0D0", linewidth = 0.15) +
  geom_violin(
    width = 0.36,
    scale = "width",
    trim = TRUE,
    adjust = 0.72,
    color = NA,
    linewidth = 0,
    alpha = 0.5
  ) +
  geom_point(
    data = plot_dt[Source == "Query"],
    aes(x = PairX, y = NormalizedExpression),
    inherit.aes = FALSE,
    position = position_jitter(width = 0.105, height = 0, seed = 20260718),
    shape = 16,
    size = 0.18,
    alpha = 0.20,
    color = "#202020"
  ) +
  geom_text(
    data = max_dt,
    aes(x = 16.88, y = y, label = MaxLabel),
    inherit.aes = FALSE,
    hjust = 0.5,
    size = 1.35,
    color = "#666666"
  ) +
  facet_grid(rows = vars(Gene), switch = "y") +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE, guide = "none") +
  scale_x_continuous(
    limits = c(0.45, 17.25),
    breaks = 1:16,
    labels = subtypes,
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(limits = c(0, 1.12), expand = expansion(mult = 0)) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.y = unit(0, "pt"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      size = 4.2, color = "#333333",
      angle = 90, vjust = 0.5, hjust = 1,
      margin = margin(t = 1)
    ),
    axis.ticks.x = element_line(color = "#555555", linewidth = 0.25),
    axis.ticks.length.x = unit(1.5, "pt"),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(
      angle = 0, hjust = 1,
      size = 5.0, color = "#222222",
      margin = margin(r = 3)
    ),
    plot.margin = margin(0, 2, 0, 2, unit = "pt")
  )

percent_y_max <- max(source_count_dt$Percent) * 1.28
p_percent <- ggplot(
  source_count_dt,
  aes(x = PairX, y = Percent, fill = Subtype)
) +
  geom_col(width = 0.34, color = NA) +
  geom_text(
    aes(y = Percent + 0.35, label = CellCount),
    angle = 90, hjust = 0, vjust = 0.5,
    size = 1.02, color = "#444444"
  ) +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE, guide = "none") +
  scale_x_continuous(
    limits = c(0.45, 17.25),
    breaks = 1:16,
    labels = subtypes,
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, percent_y_max),
    breaks = pretty(c(0, max(source_count_dt$Percent)), n = 3),
    labels = function(x) paste0(format(x, trim = TRUE), "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "% MSN") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#E2E2E2", linewidth = 0.18),
    axis.title.y = element_text(size = 4.2, color = "#444444", margin = margin(r = 2)),
    axis.text.y = element_text(size = 3.8, color = "#555555"),
    axis.ticks.y = element_line(color = "#777777", linewidth = 0.2),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 2, 1, 2, unit = "pt")
  )

p_counts <- ggplot(count_dt, aes(SubtypeIndex, N, fill = Subtype)) +
  geom_col(width = 0.72, color = NA) +
  geom_text(aes(label = N), vjust = -0.20, size = 1.12, color = "#555555") +
  scale_fill_manual(values = subtype_pal, limits = subtypes, drop = FALSE, guide = "none") +
  scale_x_continuous(limits = c(0.45, 17.15), expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.26))) +
  theme_void() +
  theme(plot.margin = margin(0, 2, 1, 40, unit = "pt"))

final_plot <- p_top / p_violin / p_counts +
  plot_layout(heights = c(0.045, 0.875, 0.080))

final_split_plot <- p_top_split / p_violin_split / p_counts +
  plot_layout(heights = c(0.060, 0.860, 0.080))

final_paired_plot <- p_top_paired / p_violin_paired / p_percent +
  plot_layout(heights = c(0.04000, 0.84375, 0.11625))

png_out <- file.path(out_dir, "Figure_MSN_RefQuery_normalized_stacked_violin_W6.3_H4.png")
pdf_out <- file.path(out_dir, "Figure_MSN_RefQuery_normalized_stacked_violin_W6.3_H4.pdf")
ggsave(png_out, final_plot, width = 6.3, height = 4.0, dpi = 600, bg = "white")
ggsave(pdf_out, final_plot, width = 6.3, height = 4.0, device = cairo_pdf, bg = "white")

split_png_out <- file.path(out_dir, "Figure_MSN_RefQuery_split_normalized_stacked_violin_W6.3_H4.png")
split_pdf_out <- file.path(out_dir, "Figure_MSN_RefQuery_split_normalized_stacked_violin_W6.3_H4.pdf")
ggsave(split_png_out, final_split_plot, width = 6.3, height = 4.0, dpi = 600, bg = "white")
ggsave(split_pdf_out, final_split_plot, width = 6.3, height = 4.0, device = cairo_pdf, bg = "white")

paired_png_out <- file.path(out_dir, "Figure_MSN_RefQuery_paired_normalized_stacked_violin_D2_8_Cpne4_QueryPointsAlpha0p2_NoTopText_W6.3_H2.67.png")
paired_pdf_out <- file.path(out_dir, "Figure_MSN_RefQuery_paired_normalized_stacked_violin_D2_8_Cpne4_QueryPointsAlpha0p2_NoTopText_W6.3_H2.67.pdf")
ggsave(paired_png_out, final_paired_plot, width = 6.3, height = 8 / 3, dpi = 600, bg = "white")
ggsave(paired_pdf_out, final_paired_plot, width = 6.3, height = 8 / 3, device = cairo_pdf, bg = "white")

fwrite(
  max_dt[, .(Gene = as.character(Gene), MaxExpression)],
  file.path(out_dir, "MSN_RefQuery_normalized_stacked_violin_gene_max.csv")
)
fwrite(
  count_dt[, .(Subtype = as.character(Subtype), CellCount = N)],
  file.path(out_dir, "MSN_RefQuery_normalized_stacked_violin_subtype_counts.csv")
)
fwrite(
  source_count_dt[, .(
    Source = as.character(Source),
    Subtype = as.character(Subtype),
    CellCount,
    SourceTotal,
    Percent
  )],
  file.path(out_dir, "MSN_RefQuery_paired_subtype_counts_and_percentages.csv")
)

readme <- c(
  "MSN Reference+Query normalized stacked violin plot",
  "",
  "Original/split canvas: width 6.3 inches; height 4.0 inches.",
  "Current paired canvas: width 6.3 inches; height 2.667 inches (two-thirds of the prior height).",
  "Reference and Query are independently normalized from the full 23,950-gene raw count matrices as log1p(count / library size * 10,000).",
  "Query subtype labels: Pearson Top750 / 721-gene / 500-bootstrap modal top1.",
  "Subtype order: D1_1-D1_8, then D2_1-D2_8.",
  "Gene order follows the supplied example from top to bottom.",
  "Ddc was replaced by Cpne4 as the final D2_8-supporting row.",
  "Within each gene, display values are divided by the combined Reference+Query maximum.",
  "The right-side number is the original combined log-normalized maximum for that gene.",
  "Split version: Reference is the left half; Query is the right half; Query uses lighter alpha.",
  "Bottom bars in the original and split versions report total Reference+Query cell counts by subtype.",
  "Paired version: each subtype has two full violins, Reference left and Query right, with identical subtype color and opacity.",
  "Reference and Query violin fills use alpha = 0.5; the top color strip remains opaque.",
  "Individual Query cells are overlaid as small dark points on the right-hand violins only, with alpha = 0.2.",
  "The centered Reference/Query text above the subtype color strip is omitted.",
  "The top subtype color strip is compressed to reduce unused vertical space.",
  "The paired-version bottom percentage panel is approximately 0.31 inches high, half of its prior physical height.",
  "Paired-version bottom bars show subtype percentage within the total MSN population of each source separately; labels give exact cell counts."
)
writeLines(readme, file.path(out_dir, "README_MSN_RefQuery_normalized_stacked_violin.txt"))

message("Output: ", out_dir)
message("Reference cells: ", ncol(ref_expr))
message("Query cells: ", ncol(query_expr))
message("Genes: ", length(genes))
