#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(circlize)
})

options(stringsAsFactors = FALSE)

pca_file <- file.path(
  "macaque_m", "m18_adaptive_pca126", "04_pca_scores.csv"
)
assignment_file <- file.path(
  "macaque_m", "m18_tempfreeze_NPC5_HCK4_res2.3",
  "01_temp_frozen_assignments_126.csv"
)
output_dir <- file.path("outputs", "Macaque_M4_circular_dendrogram")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_frozen_csv <- function(path) {
  x <- read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  names(x)[grepl("cell_label$", names(x))] <- "cell_label"
  x
}

scores <- read_frozen_csv(pca_file)
labels <- read_frozen_csv(assignment_file)
pc_names <- paste0("PC", 1:5)

stopifnot(
  nrow(scores) == 126L,
  nrow(labels) == 126L,
  !anyDuplicated(scores$cell_label),
  !anyDuplicated(labels$cell_label),
  setequal(scores$cell_label, labels$cell_label),
  all(pc_names %in% names(scores)),
  all(c("HC_K4", "GC_merged_K4", "concordant") %in% names(labels))
)

# Preserve the frozen PCA row order.  Neither the HC nor GC labels participate
# in tree construction or leaf ordering.
labels <- labels[match(scores$cell_label, labels$cell_label), , drop = FALSE]
pcs <- as.matrix(scores[, pc_names, drop = FALSE])
storage.mode(pcs) <- "double"
rownames(pcs) <- scores$cell_label
stopifnot(!anyNA(pcs), !anyNA(labels$cell_label))

hc_class <- setNames(paste0("M", labels$HC_K4), labels$cell_label)
gc_class <- setNames(paste0("M", labels$GC_merged_K4), labels$cell_label)
consensus <- setNames(as.logical(labels$concordant), labels$cell_label)
palette <- c(
  M1 = "#1F77B4", M2 = "#D9A400",
  M3 = "#8C564B", M4 = "#E377C2"
)
stopifnot(all(hc_class %in% names(palette)), all(gc_class %in% names(palette)))

hc_fit <- hclust(dist(pcs, method = "euclidean"), method = "ward.D2")
dend <- as.dendrogram(hc_fit)
ordered_ids <- labels(dend)
n_cells <- length(ordered_ids)
tree_height <- attr(dend, "height")

# edgePar belongs to the child edge.  Parent height therefore makes central
# trunks thick and tapers continuously toward terminal branches.
style_edges <- function(node, parent_height, max_height,
                        min_pt = 0.18, max_pt = 2.00, gamma = 0.72) {
  edge_par <- attr(node, "edgePar")
  if (is.null(edge_par)) edge_par <- list()
  height_fraction <- if (max_height > 0) parent_height / max_height else 0
  width_pt <- min_pt + (max_pt - min_pt) * height_fraction^gamma
  edge_par$col <- "#707070"
  edge_par$lwd <- width_pt / 0.75
  attr(node, "edgePar") <- edge_par

  if (!is.leaf(node)) {
    node_height <- attr(node, "height")
    for (i in seq_along(node)) {
      node[[i]] <- style_edges(node[[i]], node_height, max_height,
                               min_pt, max_pt, gamma)
    }
  }
  node
}
dend <- style_edges(dend, tree_height, tree_height)

draw_ring <- function(colours) {
  circos.trackPlotRegion(
    ylim = c(0, 1), track.height = 0.045, bg.border = NA,
    panel.fun = function(x, y) {
      for (i in seq_len(n_cells)) {
        circos.rect(i - 1, 0, i, 1, col = colours[i], border = NA)
      }
    }
  )
}

draw_core <- function() {
  on.exit(circos.clear(), add = TRUE)
  par(mar = c(0, 0, 0, 0), family = "Arial", xpd = NA)
  circos.clear()
  circos.par(
    start.degree = 90,
    gap.degree = 0,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(0.002, 0.002),
    circle.margin = c(0.015, 0.015, 0.015, 0.015),
    points.overflow.warning = FALSE
  )
  circos.initialize(factors = "All_cells", xlim = c(0, n_cells))

  # Tracks are created from outside inward: merged GC outer, HC inner.
  draw_ring(unname(palette[gc_class[ordered_ids]]))
  draw_ring(unname(palette[hc_class[ordered_ids]]))
  circos.trackPlotRegion(
    ylim = c(0, tree_height), track.height = 0.82, bg.border = NA,
    panel.fun = function(x, y) {
      circos.dendrogram(dend, facing = "outside", max_height = tree_height)
    }
  )
}

clean_png <- file.path(output_dir, "Macaque_M4_circular_HCinner_GCouter_clean_1in.png")
clean_pdf <- file.path(output_dir, "Macaque_M4_circular_HCinner_GCouter_clean_1in.pdf")
annot_png <- file.path(output_dir, "Macaque_M4_circular_HCinner_GCouter_annotated.png")
annot_pdf <- file.path(output_dir, "Macaque_M4_circular_HCinner_GCouter_annotated.pdf")

open_png <- function(path, width_in, height_in) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = width_in, height = height_in, units = "in",
                  res = 600, background = "white")
  } else {
    png(path, width = width_in, height = height_in, units = "in", res = 600,
        type = "cairo", bg = "white")
  }
}

open_png(clean_png, 1, 1)
draw_core()
dev.off()

cairo_pdf(clean_pdf, width = 1, height = 1, family = "Arial", bg = "white")
draw_core()
dev.off()

# The annotated canvas is expanded around an unchanged physical 1 x 1 inch
# core region, so the circular tree retains the clean panel's actual diameter.
annot_w <- 1.82
annot_h <- 1.16
draw_annotated <- function() {
  core_x0 <- 0.02 / annot_w
  core_x1 <- 1.02 / annot_w
  core_y0 <- 0.02 / annot_h
  core_y1 <- 1.02 / annot_h
  par(fig = c(core_x0, core_x1, core_y0, core_y1), mar = c(0, 0, 0, 0),
      family = "Arial")
  draw_core()

  par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE,
      family = "Arial", xpd = NA)
  plot.new()
  text(0.28, 0.965, "Macaque M1-M4 morphology hierarchy",
       cex = 4 / 12, font = 2)
  text(0.735, 0.74, "Inner: HC\nOuter: merged GC",
       cex = 4 / 12, adj = c(0, 0.5))
  text(0.735, 0.60,
       sprintf("Agreement: %d/%d (%.2f%%)", sum(consensus), n_cells,
               100 * mean(consensus)),
       cex = 4 / 12, adj = c(0, 0.5))
  legend(
    x = 0.735, y = 0.51, legend = names(palette), fill = unname(palette),
    border = NA, bty = "n", cex = 4 / 12, x.intersp = 0.45,
    y.intersp = 0.82, text.width = 0.06
  )
}

open_png(annot_png, annot_w, annot_h)
draw_annotated()
dev.off()

cairo_pdf(annot_pdf, width = annot_w, height = annot_h,
          family = "Arial", bg = "white")
draw_annotated()
dev.off()

leaf_table <- data.frame(
  circular_leaf_order = seq_len(n_cells),
  cell_label = ordered_ids,
  HC_M_class = unname(hc_class[ordered_ids]),
  merged_GC_M_class = unname(gc_class[ordered_ids]),
  consensus = unname(consensus[ordered_ids]),
  check.names = FALSE
)
write.csv(
  leaf_table,
  file.path(output_dir, "Macaque_M4_circular_leaf_order.csv"),
  row.names = FALSE
)

consensus_counts <- table(hc_class[consensus])
all_hc_counts <- table(hc_class)
all_gc_counts <- table(gc_class)
writeLines(
  c(
    "Macaque M1-M4 circular dendrogram audit",
    "========================================",
    "Tree input: 126 frozen Macaque MSN morphology-complete cells",
    "PCA input: frozen PC1-PC5 scores",
    "Clustering: Euclidean distance; Ward.D2; full HC tree",
    "Leaf order: labels(as.dendrogram(hclust(...))); no class-based reordering or rotation",
    sprintf("HC counts (all 126): %s", paste(names(all_hc_counts), all_hc_counts, collapse = ", ")),
    sprintf("Merged GC counts (all 126): %s", paste(names(all_gc_counts), all_gc_counts, collapse = ", ")),
    sprintf("Consensus HC counts (117): %s", paste(names(consensus_counts), consensus_counts, collapse = ", ")),
    sprintf("HC-GC agreement: %d/%d (%.2f%%)", sum(consensus), n_cells, 100 * mean(consensus)),
    "Ring order: inner HC; outer merged GC",
    "Palette: M1 #1F77B4; M2 #D9A400; M3 #8C564B; M4 #E377C2",
    "Reference E files were used for visual geometry only; no E features or E PCA were read."
  ),
  file.path(output_dir, "Macaque_M4_circular_audit.txt")
)

file.copy(
  normalizePath(file.path("scripts", "plot_macaque_M4_circular_dendrogram.R"),
                winslash = "/", mustWork = TRUE),
  file.path(output_dir, "plot_macaque_M4_circular_dendrogram.R"),
  overwrite = TRUE
)

cat(normalizePath(output_dir, winslash = "/"), "\n")
