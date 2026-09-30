#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
  library(ggplot2)
})

ref_file <- file.path(
  "C:/Users/53461/Documents/Codex/2026-07-13/zh/work/05_input",
  "05_Codex_branch_specific_mapping_input/00_original_inputs",
  "YZ_all_final_fixed_RNA_clean.rds"
)
query_file <- file.path(
  "C:/Users/53461/Documents/Codex/2026-07-13/zh/work/05_input",
  "05_Codex_branch_specific_mapping_input/00_original_inputs",
  "query_geneBody_symbol_raw_QC_no39cells_noO9273_Batch2_noO9268_Batch2_n641.rds"
)
pearson_qc_file <- "C:/Users/53461/Downloads/Top15_risk_x040_x050_QC_table.csv"
npc5_qc_file <- "C:/Users/53461/Downloads/Top15CleanMarker_MaxPC_NPC5_QC_table.csv"
out_dir <- file.path("outputs", "T-type_QC_RiskCells")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes <- c("Olig1", "Aqp4", "Siglech", "C1qb", "Snap25")
ref_obj <- readRDS(ref_file)
query_obj <- readRDS(query_file)

ref_meta <- ref_obj[[]]
ref_cells <- rownames(ref_meta)[as.character(ref_meta$CellType) %in% c("D1", "D2", "IN")]

ref_counts <- LayerData(ref_obj, assay = "RNA", layer = "counts")[, ref_cells,
  drop = FALSE]
query_counts <- LayerData(query_obj, assay = "RNA", layer = "counts")

stopifnot(
  all(genes %in% rownames(ref_counts)),
  all(genes %in% rownames(query_counts))
)

normalize_genes <- function(counts, genes, library_size) {
  cell_ids <- colnames(counts)
  scale_factor <- 10000 / pmax(library_size, 1)
  x <- counts[genes, , drop = FALSE] %*%
    Matrix::Diagonal(x = scale_factor)
  dimnames(x) <- list(genes, cell_ids)
  x@x <- log1p(x@x)
  as.matrix(x)
}

ref_library_size <- as.numeric(ref_meta[colnames(ref_counts), "nCount_RNA"])
query_meta <- query_obj[[]]
query_library_size <- as.numeric(query_meta[colnames(query_counts), "nCount_RNA"])

ref_expr <- normalize_genes(ref_counts, genes, ref_library_size)
query_expr <- normalize_genes(query_counts, genes, query_library_size)

rm(ref_obj, query_obj, ref_counts, query_counts)
gc(verbose = FALSE)

pearson_qc <- read.csv(pearson_qc_file, check.names = FALSE,
  stringsAsFactors = FALSE)
npc5_qc <- read.csv(npc5_qc_file, check.names = FALSE,
  stringsAsFactors = FALSE)

qc <- merge(
  npc5_qc[, c("Query_cell", "risk_class")],
  pearson_qc[, c("Query_cell", "risk_class")],
  by = "Query_cell", suffixes = c("_NPC5", "_Pearson"), all = FALSE
)
qc$QC_level <- ifelse(
  qc$risk_class_NPC5 == "retained" & qc$risk_class_Pearson == "retained",
  "Patch-seq (Good)", "Patch-seq (Risk)"
)
query_group <- setNames(qc$QC_level, qc$Query_cell)[colnames(query_expr)]

stopifnot(
  !anyNA(query_group),
  sum(query_group == "Patch-seq (Good)") == 622L,
  sum(query_group == "Patch-seq (Risk)") == 19L
)

ref_dt <- data.table(
  Gene = rep(genes, times = ncol(ref_expr)),
  Expression = as.vector(ref_expr),
  Group = "Ref (Chen et al. 2021, neuron)"
)
query_dt <- data.table(
  Gene = rep(genes, times = ncol(query_expr)),
  Expression = as.vector(query_expr),
  Group = rep(query_group, each = length(genes))
)
plot_dt <- rbind(ref_dt, query_dt)

group_levels <- c(
  "Ref (Chen et al. 2021, neuron)",
  "Patch-seq (Good)",
  "Patch-seq (Risk)"
)
plot_dt[, Gene := factor(Gene, levels = genes)]
plot_dt[, Group := factor(Group, levels = group_levels)]
group_offset <- c(
  "Ref (Chen et al. 2021, neuron)" = -0.25,
  "Patch-seq (Good)" = 0,
  "Patch-seq (Risk)" = 0.25
)
plot_dt[, GeneIndex := match(as.character(Gene), genes)]
plot_dt[, X := GeneIndex + unname(group_offset[as.character(Group)])]

p <- ggplot(plot_dt, aes(X, Expression, color = Group)) +
  geom_point(
    data = plot_dt[Group == "Ref (Chen et al. 2021, neuron)"],
    position = position_jitter(width = 0.065, height = 0, seed = 20260720),
    size = 1, alpha = 0.1, shape = 16
  ) +
  geom_point(
    data = plot_dt[Group != "Ref (Chen et al. 2021, neuron)"],
    position = position_jitter(width = 0.065, height = 0, seed = 20260721),
    size = 1, alpha = 0.5, shape = 16
  ) +
  scale_color_manual(
    values = c(
      "Ref (Chen et al. 2021, neuron)" = "#0000FF",
      "Patch-seq (Good)" = "#1F1F1F",
      "Patch-seq (Risk)" = "#D62728"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = seq_along(genes),
    labels = genes,
    expand = expansion(add = 0.48)
  ) +
  scale_y_continuous(
    breaks = seq(0, 8, by = 2),
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  labs(
    x = NULL,
    y = "Gene expression (log-normalized)",
    color = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    axis.title = element_text(size = 7.2, color = "black"),
    axis.text.y = element_text(size = 6.0, color = "black"),
    axis.text.x = element_text(size = 6.0, color = "black",
      angle = 90, hjust = 1, vjust = 0.5),
    axis.line = element_line(linewidth = 0.3528, color = "black"),
    axis.ticks = element_line(linewidth = 0.3528, color = "black"),
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.99),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 5.2),
    legend.key.width = unit(5, "pt"),
    legend.key.height = unit(5, "pt"),
    legend.spacing.y = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.background = element_rect(fill = "white", color = NA),
    aspect.ratio = 1 / 1.4,
    plot.margin = margin(2, 6, 2, 2, unit = "pt")
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 1.5, alpha = c(0.1, 0.5, 0.5))))

png_file <- file.path(out_dir,
  "Figure_marker_expression_Ref_Good_Risk_three_columns.png")
pdf_file <- file.path(out_dir,
  "Figure_marker_expression_Ref_Good_Risk_three_columns.pdf")
data_file <- file.path(out_dir,
  "Figure_marker_expression_Ref_Good_Risk_three_columns_plot_data.csv")

ggsave(png_file, p, width = 2.8, height = 2, units = "in", dpi = 600,
  bg = "white")
ggsave(pdf_file, p, width = 2.8, height = 2, units = "in",
  device = cairo_pdf, bg = "white")
fwrite(plot_dt, data_file)

message("Reference cells: ", ncol(ref_expr))
message("Query Good: ", sum(query_group == "Patch-seq (Good)"))
message("Query Risk: ", sum(query_group == "Patch-seq (Risk)"))
message("PNG: ", normalizePath(png_file))
message("PDF: ", normalizePath(pdf_file))
