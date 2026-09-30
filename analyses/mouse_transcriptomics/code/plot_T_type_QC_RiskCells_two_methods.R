#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

input_dir <- "C:/Users/53461/Downloads"
output_dir <- file.path("outputs", "T-type_QC_RiskCells")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_qc <- function(path, method, npc = NA_integer_) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  x$Method <- method
  if (!"NPC" %in% names(x)) x$NPC <- npc
  x$QC_level <- ifelse(x$risk_class == "retained", "Good", "Risk")
  x$QC_level <- factor(x$QC_level, levels = c("Good", "Risk"))
  x$risk_class <- factor(x$risk_class,
    levels = c("retained", "low-risk", "high-risk"))
  x
}

maxpc_files <- c(
  `5` = "Top15CleanMarker_MaxPC_NPC5_QC_table.csv",
  `10` = "Top15CleanMarker_MaxPC_NPC10_QC_table.csv",
  `15` = "Top15CleanMarker_MaxPC_NPC15_QC_table.csv"
)

maxpc <- do.call(rbind, lapply(names(maxpc_files), function(npc) {
  read_qc(file.path(input_dir, maxpc_files[[npc]]),
    method = "Top15CleanMarker MaxPC", npc = as.integer(npc))
}))
maxpc$NPC_label <- factor(paste0("NPC = ", maxpc$NPC),
  levels = paste0("NPC = ", c(5, 10, 15)))

risk_threshold <- read_qc(
  file.path(input_dir, "Top15_risk_x040_x050_QC_table.csv"),
  method = "Top15 risk x040/x050"
)

class_colors <- c(
  retained = "#3B82F6",
  `low-risk` = "#F59E0B",
  `high-risk` = "#DC2626"
)

base_scatter <- function(dat, title) {
  ggplot(dat, aes(plot_corr_neural, plot_corr_nonneural, color = risk_class)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
      linewidth = 0.45, color = "#6B7280") +
    geom_point(alpha = 0.78, size = 1.65) +
    scale_color_manual(values = class_colors, drop = FALSE) +
    coord_equal() +
    labs(
      title = title,
      subtitle = "Good = retained; Risk = low-risk + high-risk",
      x = "Best neural correlation",
      y = "Best non-neural correlation",
      color = "Original QC class"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "#4B5563"),
      strip.background = element_rect(fill = "#EEF2F7", color = "#9CA3AF"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

p_maxpc <- base_scatter(maxpc, "Method 1: Top15CleanMarker MaxPC") +
  facet_grid(QC_level ~ NPC_label, scales = "fixed")

p_risk <- base_scatter(risk_threshold, "Method 2: Top15 risk x040/x050") +
  facet_wrap(~QC_level, nrow = 1, scales = "fixed") +
  geom_text_repel(
    data = subset(risk_threshold, QC_level == "Risk"),
    aes(label = Query_cell), size = 2.6, max.overlaps = Inf,
    box.padding = 0.25, point.padding = 0.15, min.segment.length = 0,
    show.legend = FALSE, seed = 25
  )

summary_counts <- rbind(
  transform(as.data.frame(table(maxpc$NPC_label, maxpc$QC_level)),
    Method = "Top15CleanMarker MaxPC", Setting = Var1,
    QC_level = Var2, Count = Freq)[, c("Method", "Setting", "QC_level", "Count")],
  transform(as.data.frame(table(risk_threshold$QC_level)),
    Method = "Top15 risk x040/x050", Setting = "x040/x050",
    QC_level = Var1, Count = Freq)[, c("Method", "Setting", "QC_level", "Count")]
)
summary_counts$QC_level <- factor(summary_counts$QC_level,
  levels = c("Good", "Risk"))
summary_counts$Percent <- ave(summary_counts$Count,
  interaction(summary_counts$Method, summary_counts$Setting),
  FUN = function(z) 100 * z / sum(z))

p_summary <- ggplot(summary_counts,
    aes(Setting, Percent, fill = QC_level)) +
  geom_col(width = 0.68, color = "white") +
  geom_text(data = subset(summary_counts, QC_level == "Good"),
    aes(y = 53, label = sprintf("%d\n(%.1f%%)", Count, Percent)),
    size = 3.4, color = "white", fontface = "bold") +
  geom_text(data = subset(summary_counts, QC_level == "Risk"),
    aes(y = 6.5, label = sprintf("Risk: %d (%.1f%%)", Count, Percent)),
    size = 3.15, color = "#991B1B", fontface = "bold") +
  facet_wrap(~Method, scales = "free_x") +
  scale_fill_manual(values = c(Good = "#3B82F6", Risk = "#DC2626")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(title = "Good and Risk proportions across the two QC methods",
    x = NULL, y = "Cells (%)", fill = "QC level") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "#EEF2F7"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(output_dir, "Figure_1_Top15CleanMarker_MaxPC_Good_Risk.png"),
  p_maxpc, width = 12, height = 8, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "Figure_1_Top15CleanMarker_MaxPC_Good_Risk.pdf"),
  p_maxpc, width = 12, height = 8, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure_2_Top15_risk_x040_x050_Good_Risk.png"),
  p_risk, width = 11, height = 5.6, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "Figure_2_Top15_risk_x040_x050_Good_Risk.pdf"),
  p_risk, width = 11, height = 5.6, device = cairo_pdf)
ggsave(file.path(output_dir, "Figure_3_TwoMethods_Good_Risk_summary.png"),
  p_summary, width = 10, height = 5.8, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "Figure_3_TwoMethods_Good_Risk_summary.pdf"),
  p_summary, width = 10, height = 5.8, device = cairo_pdf)

write.csv(summary_counts,
  file.path(output_dir, "TwoMethods_Good_Risk_counts.csv"), row.names = FALSE)

message("Saved figures and summary to: ", normalizePath(output_dir))
