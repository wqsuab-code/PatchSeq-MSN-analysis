#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

input_file <- file.path(
  "outputs", "morph_qc", "seurat_npc4_res2_tsne_separation_optimization",
  "ec_hc_gc_consensus_clean", "06_cellwise_EC_HC_GC_consensus_215cells.csv"
)
output_dir <- file.path(
  "outputs", "morph_qc", "triple_consensus189_M1_M4_D1D2_100pct_stacked"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- fread(input_file)
dat <- dat[EC_HC_GC_consensus == TRUE]
stopifnot(nrow(dat) == 189L)

gc_to_m <- c(`S-1` = "M1", `S-3` = "M2", `S-2` = "M3", `S-4` = "M4")
dat[, M_class := unname(gc_to_m[Merged_class])]
stopifnot(!anyNA(dat$M_class))

plot_dt <- rbindlist(list(
  dat[, .(Panel = M_class, D1D2)]
))
# The main figure shows only cells with a stable D1 or D2 assignment. Percentages
# are renormalized within each bar after excluding D1/D2-unstable cells.
plot_dt <- plot_dt[D1D2 %in% c("D1", "D2")]
panel_levels <- c("M1", "M2", "M3", "M4")
d_levels <- c("D1", "D2")
plot_dt[, Panel := factor(Panel, levels = panel_levels)]
plot_dt[, D1D2 := factor(D1D2, levels = d_levels)]

counts <- plot_dt[, .N, by = .(Panel, D1D2)]
counts <- CJ(
  Panel = factor(panel_levels, levels = panel_levels),
  D1D2 = factor(d_levels, levels = d_levels), unique = TRUE
)[counts, on = .(Panel, D1D2)]
counts[is.na(N), N := 0L]
counts[, Total := sum(N), by = Panel]
counts[, Percent := 100 * N / Total]
setorder(counts, Panel, D1D2)

palette <- c(
  D1 = "#E41A1C",
  D2 = "#377EB8"
)

make_plot <- function(show_text = TRUE) {
  p <- ggplot(counts, aes(x = Panel, y = Percent, fill = D1D2)) +
    geom_col(
      width = 0.82, colour = "white", linewidth = 0.28,
      position = position_stack(reverse = TRUE)
    ) +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
    # Bar centres are 1:4 and bar width is 0.82, so the outer bar edges are
    # 0.59 and 4.41. Limits 0.426 and 4.574 leave exactly 0.164 category units
    # (=20% of the 0.82 bar width) outside the first and last bars.
    scale_x_discrete(expand = expansion(add = c(0, 0))) +
    coord_cartesian(xlim = c(0.426, 4.574), clip = "off") +
    theme_classic(base_size = 7) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.32),
      axis.ticks = element_line(colour = "black", linewidth = 0.30),
      axis.ticks.length = grid::unit(1.6, "pt"),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    labs(x = NULL, y = NULL, fill = NULL)
  if (show_text) {
    p <- p +
      scale_y_continuous(
        limits = c(0, 100), breaks = c(0, 50, 100),
        labels = c("0", "50", "100"), expand = expansion(mult = c(0, 0))
      ) +
      labs(y = "Cells (%)") +
      theme(
        axis.text.x = element_text(colour = "black", size = 6.5),
        axis.text.y = element_text(colour = "black", size = 6),
        axis.title.y = element_text(colour = "black", size = 6.5, margin = margin(r = 3)),
        legend.position = "top",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.text = element_text(size = 5.3, colour = "black"),
        legend.key.size = grid::unit(6.5, "pt"),
        legend.spacing.x = grid::unit(2, "pt"),
        legend.margin = margin(0, 0, 1, 0)
      )
  } else {
    p <- p + theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line.y = element_blank()
    )
  }
  p
}

p_labeled <- make_plot(TRUE)
p_clean <- make_plot(FALSE)

save_pair <- function(stem, plot, width, height) {
  ggsave(
    file.path(output_dir, paste0(stem, ".png")), plot,
    width = width, height = height, units = "in", dpi = 600, bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "in", device = cairo_pdf, bg = "white"
  )
}

save_pair("01_M1_M4_D1D2_100pct_stacked_labeled", p_labeled, 1.75, 1.55)
save_pair("02_M1_M4_D1D2_100pct_stacked_clean", p_clean, 1.20, 1.20)

fwrite(counts, file.path(output_dir, "03_M1_M4_D1D2_counts_percentages.csv"))
fwrite(dat[, .(MSN_unique_ID, M_class, D1D2)], file.path(output_dir, "04_cell_assignments_189.csv"))
writeLines(c(
  "100% stacked D1/D2 composition plot for EC-HC-GC triple-consensus Morph cells with stable D1/D2 assignments.",
  "Bars from left to right: M1, M2, M3, M4.",
  "D1: red (#E41A1C); D2: blue (#377EB8).",
  "Five D1/D2-unstable cells are excluded; each bar is renormalized to stable D1 + D2 = 100%.",
  "PNG resolution: 600 dpi; vector PDF also supplied."
), file.path(output_dir, "05_methods.txt"), useBytes = TRUE)

print(counts)
cat("Output:", normalizePath(output_dir, winslash = "/"), "\n")
