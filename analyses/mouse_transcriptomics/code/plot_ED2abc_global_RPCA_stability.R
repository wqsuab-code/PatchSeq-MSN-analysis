#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggalluvial)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
bundle_dir <- dirname(dirname(script_path))
global_dir <- file.path(bundle_dir, "data", "global")
out_dir <- file.path(bundle_dir, "figures", "extended_data_2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scan <- fread(file.path(global_dir, "stage02B_RPCA_scan_summary.csv"))
pred <- fread(file.path(global_dir, "stage02B_RPCA_all_predictions.csv"))
fine <- fread(file.path(global_dir, "03A_final_RPCA_CellType_summary.csv"))
root <- fread(file.path(global_dir, "03B_final_RPCA_Root_summary.csv"))
final_celltype <- fread(file.path(
  global_dir, "03A_final_RPCA_per_cell_CellType_classification.csv"
))[, .(cell_id, Final_consensus_CellType)]

npc_levels <- c(10L, 12L, 14L, 16L, 18L, 20L, 22L, 24L, 26L, 28L, 30L)
stopifnot(
  identical(scan$nPC, npc_levels),
  nrow(pred) == 641L * length(npc_levels),
  identical(sort(unique(pred$nPC)), npc_levels),
  range(scan$anchor_count) == c(654L, 1003L)
)

major_colors <- c(D1 = "#E31A1C", D2 = "#377EB8", IN = "#00A63C")
gray_variable <- "#BDBDBD"
alluvium_colors <- c(major_colors, Ambiguous = gray_variable)

theme_panel <- theme_classic(base_size = 7.2) +
  theme(
    axis.title = element_text(size = 7.2, color = "black"),
    axis.text = element_text(size = 6.3, color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    legend.text = element_text(size = 6.2),
    legend.title = element_blank(),
    plot.margin = margin(3, 4, 3, 4, unit = "pt")
  )

scan[, npc_label := factor(paste0("PC", nPC), levels = paste0("PC", rev(npc_levels)))]
p_a <- ggplot(scan, aes(anchor_count, npc_label)) +
  geom_col(width = 0.72, fill = "#D9D9D9", color = "black", linewidth = 0.35) +
  scale_x_continuous(
    limits = c(0, 1050), breaks = c(0, 500, 1000),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(x = "Anchor count", y = NULL, tag = "a") +
  theme_panel +
  theme(panel.grid = element_blank())

pred[, nPC_factor := factor(nPC, levels = npc_levels)]
pred[, predicted_CellType := factor(
  predicted_CellType,
  levels = c("D1", "D2", "IN")
)]
pred <- merge(pred, final_celltype, by = "cell_id", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(pred$Final_consensus_CellType))

p_b <- ggplot(
  pred,
  aes(x = nPC_factor, stratum = predicted_CellType, alluvium = cell_id, y = 1)
) +
  geom_alluvium(
    aes(fill = Final_consensus_CellType),
    color = NA, alpha = 0.40, width = 0.20, knot.pos = 0.45
  ) +
  geom_stratum(
    aes(fill = predicted_CellType),
    width = 0.42, color = "white", linewidth = 0.15
  ) +
  scale_fill_manual(values = alluvium_colors, drop = FALSE) +
  scale_x_discrete(labels = as.character(npc_levels), expand = c(0.02, 0.02)) +
  labs(x = "Number of PCs", y = "Patch-seq cells", fill = NULL, tag = "b") +
  theme_panel +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 5.8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(linewidth = 0.35)
  )

fine_counts <- data.table(
  level = factor(c("D1", "D2", "IN"), levels = c("D1", "D2", "IN")),
  stable = c(272L, 279L, 34L),
  variable = c(28L, 18L, 10L)
)
root_counts <- data.table(
  level = factor(c("MSN", "IN"), levels = c("MSN", "IN")),
  stable = c(588L, 34L),
  variable = c(9L, 10L)
)

stopifnot(
  sum(fine_counts$stable) == 585L,
  sum(fine_counts$stable + fine_counts$variable) == 641L,
  sum(root_counts$stable) == 622L,
  sum(root_counts$stable + root_counts$variable) == 641L
)

fine_long <- melt(
  fine_counts, id.vars = "level", variable.name = "status", value.name = "n"
)
fine_long[, fill_key := fifelse(status == "stable", as.character(level), "Variable")]
fine_labels <- fine_counts[, .(
  level, y = stable + variable + 12,
  label = sprintf("%d\n(%d)", stable, variable)
)]

root_long <- melt(
  root_counts, id.vars = "level", variable.name = "status", value.name = "n"
)
root_long[, fill_key := fifelse(
  status == "stable" & level == "MSN", "MSN",
  fifelse(status == "stable" & level == "IN", "IN", "Variable")
)]
root_labels <- root_counts[, .(
  level, y = stable + variable + 25,
  label = sprintf("%d\n(%d)", stable, variable)
)]

p_c1 <- ggplot(fine_long, aes(level, n, fill = fill_key)) +
  geom_col(
    width = 0.72, color = "black", linewidth = 0.25,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = fine_labels, aes(level, y, label = label),
    inherit.aes = FALSE, size = 2.2, lineheight = 0.85
  ) +
  scale_fill_manual(values = c(major_colors, Variable = gray_variable), guide = "none") +
  scale_y_continuous(limits = c(0, 335), breaks = c(0, 100, 200, 300), expand = c(0, 0)) +
  labs(x = NULL, y = "Cell counts", tag = "c") +
  theme_panel

p_c2 <- ggplot(root_long, aes(level, n, fill = fill_key)) +
  geom_col(
    width = 0.72, color = "black", linewidth = 0.25,
    position = position_stack(reverse = TRUE)
  ) +
  geom_text(
    data = root_labels, aes(level, y, label = label),
    inherit.aes = FALSE, size = 2.2, lineheight = 0.85
  ) +
  scale_fill_manual(
    values = c(MSN = "#8B4A6B", IN = major_colors[["IN"]], Variable = gray_variable),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(0, 650), breaks = c(0, 200, 400, 600), expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_panel

p_c <- p_c1 + p_c2 + plot_layout(widths = c(1.05, 1))

combined <- p_a + p_b + p_c +
  plot_layout(widths = c(0.78, 1.55, 1.28)) +
  plot_annotation(
    title = "Major cell-type mapping (RPCA)",
    theme = theme(
      plot.title = element_text(size = 7.5, hjust = 0.5, margin = margin(b = 3)),
      plot.tag = element_text(size = 8, face = "bold")
    )
  )

ggsave(
  file.path(out_dir, "ED2abc_global_RPCA_stability_correct_nPC.png"),
  combined, width = 7.2, height = 2.25, dpi = 600, bg = "white"
)
ggsave(
  file.path(out_dir, "ED2abc_global_RPCA_stability_correct_nPC.pdf"),
  combined, width = 7.2, height = 2.25, device = cairo_pdf, bg = "white"
)

fwrite(
  rbind(
    fine_counts[, .(classification = "CellType", level = as.character(level), stable, variable)],
    root_counts[, .(classification = "Root", level = as.character(level), stable, variable)]
  ),
  file.path(global_dir, "ED2c_stable_variable_counts.csv")
)

message("Wrote corrected ED Fig. 2a-c with actual nPC values: ", paste(npc_levels, collapse = ", "))
