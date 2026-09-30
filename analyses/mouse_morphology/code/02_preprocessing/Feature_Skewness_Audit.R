#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(ggplot2)
})

input_file <- file.path(".codex-work", "e_type_qc", "stage1_filtered.json")
output_dir <- file.path("outputs", "e_type_qc", "feature_redundancy_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

png_file <- file.path(output_dir, "Feature_Skewness_Audit.png")
pdf_file <- file.path(output_dir, "Feature_Skewness_Audit.pdf")
csv_file <- file.path(output_dir, "Feature_Skewness_Audit_values.csv")

payload <- fromJSON(input_file, simplifyVector = FALSE)
headers <- unlist(payload$analysis_headers)
matrix_data <- do.call(rbind, lapply(payload$analysis_rows, function(row) {
  vapply(seq_along(headers), function(j) {
    if (is.null(row[[j]])) NA_character_ else as.character(row[[j]])
  }, character(1))
}))
dat <- as.data.frame(matrix_data, stringsAsFactors = FALSE, check.names = FALSE)
names(dat) <- headers

# Fixed top-to-bottom order supplied by the user. Do not sort by skewness.
feature_order <- c(
  "E_Holding.MP..mV.",
  "E_Fitted.MP..mV.",
  "E_Input.resistance..MOhm.",
  "E_Membrane.time.constant..ms.",
  "E_AP.threshold..mV.",
  "E_AP.amplitude..mV.",
  "E_AP.width..ms.",
  "E_Upstroke.to.downstroke.ratio",
  "E_Afterhyperpolarization..mV.",
  "E_ISI.adaptation.index",
  "E_Max.number.of.APs",
  "E_Rheobase..pA.",
  "E_Sag.ratio",
  "E_Latency..ms.",
  "E_Latency....20pA.current..ms.",
  "E_ISI.Fano.factor",
  "E_ISI.coefficient.of.variation",
  "E_ISI.average.adaptation.index",
  "E_Rebound..mV.",
  "E_Sag.time..s.",
  "E_Sag.area..mV.s.",
  "E_AP.amplitude.adaptation.index",
  "E_AP.amplitude.average.adaptation.index",
  "E_AP.Fano.factor",
  "E_AP.coefficient.of.variation"
)

missing_features <- setdiff(feature_order, names(dat))
if (length(missing_features) > 0) {
  stop("Missing requested features: ", paste(missing_features, collapse = ", "))
}

qc_features <- c(
  "E_Holding.MP..mV.", "E_Input.resistance..MOhm.",
  "E_AP.amplitude..mV.", "E_Max.number.of.APs",
  "E_AP.width..ms.", "E_Rheobase..pA."
)
numeric_features <- unique(c(feature_order, qc_features))
for (feature in numeric_features) {
  dat[[feature]] <- suppressWarnings(as.numeric(dat[[feature]]))
}

qc_pass <- with(
  dat,
  E_Holding.MP..mV. < -55 &
    E_Input.resistance..MOhm. >= 100 &
    E_Input.resistance..MOhm. <= 1000 &
    E_AP.amplitude..mV. > 40 &
    E_Max.number.of.APs > 2 &
    E_AP.width..ms. <= 3 &
    E_Rheobase..pA. >= 10 &
    E_Rheobase..pA. <= 300
)
qc_pass[is.na(qc_pass)] <- FALSE
metric_df <- dat[qc_pass, feature_order, drop = FALSE]

# Fisher-Pearson moment coefficient, matching the preceding skewness QC plots.
moment_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  standard_deviation <- sd(x)
  if (!is.finite(standard_deviation) || standard_deviation == 0) return(NA_real_)
  mean(((x - mean(x)) / standard_deviation)^3)
}

skewness_values <- vapply(metric_df, moment_skewness, numeric(1))
plot_data <- data.frame(
  Order = seq_along(feature_order),
  Feature = feature_order,
  Display_label = sub("^E_", "", feature_order),
  Skewness = unname(skewness_values[feature_order]),
  stringsAsFactors = FALSE
)
plot_data$Skewness_class <- ifelse(
  abs(plot_data$Skewness) > 1,
  "|Skewness| > 1",
  "|Skewness| <= 1"
)

# ggplot draws the first factor level at the bottom, so reverse the levels to
# preserve the requested list as the visual top-to-bottom order.
plot_data$Display_label <- factor(
  plot_data$Display_label,
  levels = rev(plot_data$Display_label)
)

skewness_plot <- ggplot(
  plot_data,
  aes(x = Skewness, y = Display_label, fill = Skewness_class)
) +
  geom_col(
    width = 0.76,
    colour = "#202020",
    linewidth = 0.18
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    colour = "#B2182B",
    linetype = "dashed",
    linewidth = 0.28
  ) +
  scale_fill_manual(
    values = c(
      "|Skewness| > 1" = "#D98989",
      "|Skewness| <= 1" = "#6FA8CF"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    title = sprintf(
      "Feature Skewness Audit (N=%d, %d Metrics)",
      nrow(metric_df), ncol(metric_df)
    ),
    subtitle = "Red dashed lines indicate skewness thresholds of -1 and +1",
    x = "Skewness",
    y = NULL
  ) +
  theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#E6E6E6", linewidth = 0.16),
    panel.grid.major.x = element_line(colour = "#EFEFEF", linewidth = 0.16),
    panel.border = element_rect(colour = "#333333", fill = NA, linewidth = 0.18),
    axis.ticks = element_line(colour = "#333333", linewidth = 0.18),
    axis.text.y = element_text(size = 5.6, colour = "black"),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.title.x = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
    plot.subtitle = element_text(size = 6.5, hjust = 0.5),
    plot.margin = margin(3, 5, 3, 3)
  )

ggsave(
  png_file, skewness_plot,
  width = 9, height = 3, units = "in",
  dpi = 600, bg = "white"
)
ggsave(
  pdf_file, skewness_plot,
  width = 9, height = 3, units = "in",
  useDingbats = FALSE, bg = "white"
)

# Restore plain character labels before writing the audit table.
plot_data$Display_label <- as.character(plot_data$Display_label)
plot_data <- plot_data[order(plot_data$Order), , drop = FALSE]
write.csv(plot_data, csv_file, row.names = FALSE)

cat(sprintf(
  "Completed Feature Skewness Audit: N=%d, metrics=%d.\n",
  nrow(metric_df), ncol(metric_df)
))
cat(normalizePath(png_file, winslash = "/", mustWork = TRUE), "\n")
cat(normalizePath(pdf_file, winslash = "/", mustWork = TRUE), "\n")
cat(normalizePath(csv_file, winslash = "/", mustWork = TRUE), "\n")
