#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

input_file <- file.path(
  "outputs", "e_type_qc", "parallel_preprocessing",
  "Ephys_QCpass_494_raw_25features.csv"
)
output_dir <- file.path(
  "outputs", "e_type_qc", "pca_shift_sum_final18"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

features <- c(
  "E_Holding.MP..mV.",
  "E_Input.resistance..MOhm.",
  "E_Membrane.time.constant..ms.",
  "E_Rheobase..pA.",
  "E_Sag.ratio",
  "E_Sag.time..s.",
  "E_AP.threshold..mV.",
  "E_AP.amplitude..mV.",
  "E_AP.width..ms.",
  "E_Upstroke.to.downstroke.ratio",
  "E_Afterhyperpolarization..mV.",
  "E_Max.number.of.APs",
  "E_Latency..ms.",
  "E_Latency....20pA.current..ms.",
  "E_ISI.adaptation.index",
  "E_ISI.coefficient.of.variation",
  "E_AP.amplitude.adaptation.index",
  "E_AP.coefficient.of.variation"
)

raw_table <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
missing_features <- setdiff(features, names(raw_table))
if (length(missing_features) > 0) {
  stop("Missing final-panel features: ", paste(missing_features, collapse = ", "))
}
if (nrow(raw_table) != 494) {
  stop("Expected 494 QC-pass cells but found ", nrow(raw_table), ".")
}
if (anyDuplicated(raw_table$MSN_unique_ID)) {
  stop("MSN_unique_ID is not unique.")
}

raw_matrix <- as.matrix(raw_table[, features, drop = FALSE])
storage.mode(raw_matrix) <- "double"
rownames(raw_matrix) <- raw_table$MSN_unique_ID
if (anyNA(raw_matrix) || any(!is.finite(raw_matrix))) {
  stop("Final 18-feature matrix contains missing or non-finite values.")
}

# Transcriptome-like branch: feature-wise translation to zero, column-sum
# normalization to 10,000, log1p transformation, and feature-wise z scoring.
feature_minimum <- apply(raw_matrix, 2, min)
shifted_matrix <- sweep(raw_matrix, 2, feature_minimum, FUN = "-")
shifted_sum <- colSums(shifted_matrix)
if (any(shifted_sum <= 0 | !is.finite(shifted_sum))) {
  stop("At least one shifted feature has a non-positive or invalid column sum.")
}
normalized_matrix <- sweep(shifted_matrix, 2, shifted_sum, FUN = "/") * 10000
log1p_matrix <- log1p(normalized_matrix)
z_matrix <- scale(log1p_matrix, center = TRUE, scale = TRUE)
if (anyNA(z_matrix) || any(!is.finite(z_matrix))) {
  stop("Final log1p-z-scored PCA matrix contains missing or non-finite values.")
}

# Data are already feature-centered and scaled by the preceding z-score step.
pca_fit <- prcomp(z_matrix, center = FALSE, scale. = FALSE)
variance_percent <- 100 * pca_fit$sdev^2 / sum(pca_fit$sdev^2)
cumulative_percent <- cumsum(variance_percent)
pc_names <- paste0("PC", seq_along(variance_percent))

variance_table <- data.frame(
  PC = pc_names,
  PC_number = seq_along(variance_percent),
  Individual_variance_percent = variance_percent,
  Cumulative_variance_percent = cumulative_percent,
  stringsAsFactors = FALSE
)

scores_table <- data.frame(
  MSN_unique_ID = raw_table$MSN_unique_ID,
  pca_fit$x,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
loadings_table <- data.frame(
  Feature = rownames(pca_fit$rotation),
  pca_fit$rotation,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
zscore_table <- data.frame(
  MSN_unique_ID = raw_table$MSN_unique_ID,
  z_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
parameter_table <- data.frame(
  Feature = features,
  Original_minimum = feature_minimum[features],
  Shifted_column_sum = shifted_sum[features],
  Normalized_column_sum = colSums(normalized_matrix)[features],
  Log1p_mean = colMeans(log1p_matrix)[features],
  Log1p_sd = apply(log1p_matrix, 2, sd)[features],
  Final_zscore_mean = colMeans(z_matrix)[features],
  Final_zscore_sd = apply(z_matrix, 2, sd)[features],
  stringsAsFactors = FALSE
)

write.csv(
  variance_table,
  file.path(output_dir, "ShiftSumNorm10000_Log1p_Zscore_Final18_PCA_variance.csv"),
  row.names = FALSE
)
write.csv(
  scores_table,
  file.path(output_dir, "ShiftSumNorm10000_Log1p_Zscore_Final18_PCA_scores.csv"),
  row.names = FALSE
)
write.csv(
  loadings_table,
  file.path(output_dir, "ShiftSumNorm10000_Log1p_Zscore_Final18_PCA_loadings.csv"),
  row.names = FALSE
)
write.csv(
  zscore_table,
  file.path(output_dir, "ShiftSumNorm10000_Log1p_Zscore_Final18_PCA_input.csv"),
  row.names = FALSE
)
write.csv(
  parameter_table,
  file.path(output_dir, "ShiftSumNorm10000_Log1p_Zscore_Final18_parameters.csv"),
  row.names = FALSE
)

display_n <- nrow(variance_table)
plot_data <- variance_table[seq_len(display_n), , drop = FALSE]
plot_data$PC_factor <- factor(
  plot_data$PC_number,
  levels = plot_data$PC_number
)

first_pc_at <- function(threshold) {
  index <- which(cumulative_percent >= threshold)[1]
  if (length(index) == 0 || is.na(index)) NA_integer_ else index
}
pc70 <- first_pc_at(70)
pc80 <- first_pc_at(80)
pc90 <- first_pc_at(90)

variance_plot <- ggplot(
  plot_data,
  aes(x = PC_factor, y = Cumulative_variance_percent, group = 1)
) +
  geom_col(
    aes(y = Individual_variance_percent),
    width = 0.62,
    fill = "#B9DCE8",
    alpha = 0.78,
    color = NA
  ) +
  geom_hline(
    yintercept = 70,
    color = "grey35",
    linewidth = 0.32,
    linetype = "dotted"
  ) +
  geom_hline(
    yintercept = 80,
    color = "#8B0000",
    linewidth = 0.36,
    linetype = "dashed"
  ) +
  geom_line(color = "#2C6FA9", linewidth = 0.55) +
  geom_point(shape = 21, size = 1.35, stroke = 0.25,
             fill = "#8B0000", color = "#8B0000") +
  geom_text(
    aes(label = sprintf("%.1f", Cumulative_variance_percent)),
    vjust = -0.72,
    size = 1.48,
    fontface = "bold"
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, by = 20),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "PCA variance explained",
    subtitle = "Bars: individual | Line: cumulative (70%, 80%)",
    x = "Principal components",
    y = "Variance explained (%)"
  ) +
  theme_classic(base_size = 6) +
  theme(
    panel.grid.major = element_line(color = "grey87", linewidth = 0.22),
    panel.grid.minor = element_line(color = "grey94", linewidth = 0.18),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.40),
    axis.line = element_blank(),
    axis.text = element_text(color = "black", size = 5.1),
    axis.title = element_text(color = "black", size = 5.8),
    plot.title = element_text(face = "bold", size = 7.2, hjust = 0),
    plot.subtitle = element_text(size = 5.0, hjust = 0, margin = margin(b = 1.5)),
    plot.margin = margin(3, 3, 2.5, 3)
  )

png_file <- file.path(
  output_dir,
  "ShiftSumNorm10000_Log1p_Zscore_Final18_PCA_Variance_W3_H2.png"
)
pdf_file <- file.path(
  output_dir,
  "ShiftSumNorm10000_Log1p_Zscore_Final18_PCA_Variance_W3_H2.pdf"
)
ggsave(
  png_file,
  variance_plot,
  width = 3,
  height = 2,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  pdf_file,
  variance_plot,
  width = 3,
  height = 2,
  units = "in",
  useDingbats = FALSE,
  bg = "white"
)

log_file <- file.path(output_dir, "ShiftSumNorm10000_Log1p_Zscore_Final18_PCA_run_log.txt")
log_lines <- c(
  sprintf("Cells: %d", nrow(z_matrix)),
  sprintf("Features: %d", ncol(z_matrix)),
  "Preprocessing: feature minimum subtraction; divide by shifted feature sum; multiply by 10000; log1p; feature-wise z score.",
  "PCA: prcomp(center=FALSE, scale.=FALSE) on the already z-scored matrix.",
  sprintf("PCs required for >=70%% cumulative variance: %s", pc70),
  sprintf("PCs required for >=80%% cumulative variance: %s", pc80),
  sprintf("PCs required for >=90%% cumulative variance: %s", pc90),
  sprintf(
    "Maximum absolute normalized column-sum error: %.12g",
    max(abs(colSums(normalized_matrix) - 10000))
  ),
  sprintf(
    "Maximum absolute final z-score column mean: %.12g",
    max(abs(colMeans(z_matrix)))
  ),
  sprintf(
    "Maximum absolute final z-score SD-minus-one: %.12g",
    max(abs(apply(z_matrix, 2, sd) - 1))
  ),
  sprintf("PNG: %s", normalizePath(png_file, winslash = "/", mustWork = TRUE)),
  sprintf("PDF: %s", normalizePath(pdf_file, winslash = "/", mustWork = TRUE))
)
writeLines(log_lines, log_file, useBytes = TRUE)

cat(paste(log_lines, collapse = "\n"), "\n")
