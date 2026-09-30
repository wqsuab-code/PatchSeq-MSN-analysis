#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(corrplot)
})

input_file <- file.path(".codex-work", "e_type_qc", "stage1_filtered.json")
output_dir <- file.path("outputs", "e_type_qc", "feature_redundancy_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

map_file <- file.path(output_dir, "Feature_Abbreviation_Map.csv")
connection_thresholds <- c(0.60, 0.75)
threshold_tags <- sprintf("%03d", round(connection_thresholds * 100))

payload <- fromJSON(input_file, simplifyVector = FALSE)
headers <- unlist(payload$analysis_headers)
matrix_data <- do.call(rbind, lapply(payload$analysis_rows, function(row) {
  vapply(seq_along(headers), function(j) {
    if (is.null(row[[j]])) NA_character_ else as.character(row[[j]])
  }, character(1))
}))
dat <- as.data.frame(matrix_data, stringsAsFactors = FALSE, check.names = FALSE)
names(dat) <- headers

# Fixed feature order shared by the skewness and correlation panels.
features <- c(
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

abbreviations <- c(
  "Vhold",
  "Vfit",
  "Rin",
  "Tau_m",
  "AP_thr",
  "AP_amp",
  "AP_width",
  "UD_ratio",
  "AHP",
  "ISI_adapt",
  "Max_APs",
  "Rheobase",
  "Sag_ratio",
  "Latency",
  "Latency_-20pA",
  "ISI_Fano",
  "ISI_CV",
  "ISI_adapt_avg",
  "Rebound",
  "Sag_time",
  "Sag_area",
  "AP_amp_adapt",
  "AP_amp_adapt_avg",
  "AP_Fano",
  "AP_CV"
)

if (length(features) != length(abbreviations)) {
  stop("Feature and abbreviation vectors have different lengths.")
}
missing_features <- setdiff(features, names(dat))
if (length(missing_features) > 0) {
  stop("Missing requested features: ", paste(missing_features, collapse = ", "))
}

abbreviation_map <- data.frame(
  Order = seq_along(features),
  Full_feature = features,
  Abbreviation = abbreviations,
  stringsAsFactors = FALSE
)
write.csv(abbreviation_map, map_file, row.names = FALSE)

qc_features <- c(
  "E_Holding.MP..mV.", "E_Input.resistance..MOhm.",
  "E_AP.amplitude..mV.", "E_Max.number.of.APs",
  "E_AP.width..ms.", "E_Rheobase..pA."
)
for (feature in unique(c(features, qc_features))) {
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
metric_df <- dat[qc_pass, features, drop = FALSE]

moment_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  standard_deviation <- sd(x)
  if (!is.finite(standard_deviation) || standard_deviation == 0) return(NA_real_)
  mean(((x - mean(x)) / standard_deviation)^3)
}

skewness_values <- vapply(metric_df, moment_skewness, numeric(1))
corr_matrix <- cor(metric_df, method = "pearson", use = "pairwise.complete.obs")
corr_matrix[is.na(corr_matrix)] <- 0
diag(corr_matrix) <- 1
rownames(corr_matrix) <- abbreviations
colnames(corr_matrix) <- abbreviations

correlation_palette <- colorRampPalette(
  c(
    "#8E002B", "#D7301F", "#F4A582", "#F7F7F7",
    "#92C5DE", "#4393C3", "#053061"
  )
)(201)

plot_combined <- function(
  output_file,
  connection_threshold,
  file_type = c("png", "pdf")
) {
  file_type <- match.arg(file_type)
  connection_index <- which(
    upper.tri(corr_matrix, diag = FALSE) &
      abs(corr_matrix) > connection_threshold,
    arr.ind = TRUE
  )
  if (file_type == "png") {
    png(
      output_file,
      width = 3900,
      height = 2400,
      res = 300,
      bg = "white"
    )
  } else {
    pdf(
      output_file,
      width = 13,
      height = 8,
      useDingbats = FALSE,
      bg = "white"
    )
  }
  on.exit(dev.off(), add = TRUE)

  # Three coordinated regions: upper-triangle correlation matrix, one shared
  # centered label column, and absolute-skewness lollipops.
  layout(
    matrix(c(1, 2, 3), nrow = 1),
    widths = c(0.60, 0.14, 0.26)
  )

  # Left: upper-triangle Pearson correlation matrix. Axis labels are suppressed
  # here because the shared row labels are drawn in the middle region.
  par(bg = "white", xpd = NA)
  corr_result <- corrplot::corrplot(
    corr_matrix,
    method = "ellipse",
    type = "upper",
    order = "original",
    diag = TRUE,
    col = correlation_palette,
    tl.pos = "n",
    cl.pos = "b",
    col.lim = c(-1, 1),
    cl.length = 11,
    cl.cex = 0.72,
    addgrid.col = "grey82",
    outline = FALSE,
    mar = c(3, 1, 6, 1)
  )
  number_positions <- corr_result$corrPos
  text(
    number_positions$x,
    number_positions$y,
    labels = sprintf("%.2f", number_positions$corr),
    cex = 0.38,
    font = 2,
    col = "black"
  )

  # Curved redundancy links occupy the otherwise empty lower triangle. Their
  # endpoints sit just outside the diagonal cells, avoiding the ellipses and
  # coefficient labels. Line width is proportional to |r|.
  if (nrow(connection_index) > 0) {
    feature_count <- length(abbreviations)
    for (k in seq_len(nrow(connection_index))) {
      first_index <- connection_index[k, 1]
      second_index <- connection_index[k, 2]
      correlation_value <- corr_matrix[first_index, second_index]

      x_start <- first_index - 0.34
      y_start <- feature_count + 1 - first_index - 0.34
      x_end <- second_index - 0.34
      y_end <- feature_count + 1 - second_index - 0.34

      separation <- abs(second_index - first_index)
      # Use a deliberately generous outward bend so the redundancy links are
      # visually distinct from the diagonal cells and from one another.
      bend <- 1.00 + 0.24 * separation
      x_control <- (x_start + x_end) / 2 - bend
      y_control <- (y_start + y_end) / 2 - bend

      curve_t <- seq(0, 1, length.out = 120)
      curve_x <-
        (1 - curve_t)^2 * x_start +
        2 * (1 - curve_t) * curve_t * x_control +
        curve_t^2 * x_end
      curve_y <-
        (1 - curve_t)^2 * y_start +
        2 * (1 - curve_t) * curve_t * y_control +
        curve_t^2 * y_end

      strength_scaled <- min(
        max((abs(correlation_value) - 0.60) / 0.40, 0),
        1
      )
      connection_alpha <- 0.28 + 0.50 * strength_scaled
      connection_color <- if (correlation_value >= 0) {
        adjustcolor("#2166AC", alpha.f = connection_alpha)
      } else {
        adjustcolor("#B2182B", alpha.f = connection_alpha)
      }
      connection_width <- 0.65 + 2.35 * strength_scaled

      lines(
        curve_x,
        curve_y,
        col = connection_color,
        lwd = connection_width,
        lty = if (correlation_value >= 0) 1 else 2
      )
      points(
        c(x_start, x_end),
        c(y_start, y_end),
        pch = 16,
        cex = 0.28,
        col = connection_color
      )
    }
  }

  # Abbreviated column labels across the upper edge of the matrix.
  text(
    x = seq_along(abbreviations),
    y = length(abbreviations) + 0.65,
    labels = abbreviations,
    srt = 48,
    adj = c(0, 0),
    cex = 0.62,
    col = "black",
    xpd = NA
  )
  title(
    main = sprintf(
      "Pearson Correlation (links: |r| > %.2f; %d pairs)",
      connection_threshold,
      nrow(connection_index)
    ),
    font.main = 2,
    cex.main = 1.05,
    line = 2
  )
  correlation_usr <- par("usr")
  y_positions <- rev(seq_along(features))

  # Middle: the only row-label column used by both panels. Matching top/bottom
  # margins and the copied corrplot y-range guarantee row-center alignment.
  par(
    bg = "white",
    mar = c(3, 0, 6, 0),
    xpd = NA
  )
  plot(
    NA_real_, NA_real_,
    xlim = c(0, 1),
    ylim = correlation_usr[3:4],
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  text(
    x = 0.5,
    y = y_positions,
    labels = abbreviations,
    adj = c(0.5, 0.5),
    cex = 0.72,
    col = "black"
  )

  # Right: lollipop length encodes absolute skewness. Color retains the signed
  # direction, with separate negative/positive saturation scales so the most
  # negative value is deep blue, zero is gray, and the most positive is red.
  par(
    bg = "white",
    mar = c(3, 0.5, 6, 2),
    xpd = FALSE
  )
  absolute_skewness <- abs(skewness_values)
  skew_max <- max(c(absolute_skewness, 1), na.rm = TRUE)
  skew_xlim <- c(0, skew_max * 1.08)
  plot(
    NA_real_, NA_real_,
    xlim = skew_xlim,
    ylim = correlation_usr[3:4],
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  abline(h = y_positions, col = "grey92", lwd = 0.6)
  abline(v = 0, col = "grey45", lwd = 0.9)
  abline(v = 1, col = "#B2182B", lty = 2, lwd = 1)

  negative_palette <- colorRampPalette(c("#BDBDBD", "#2166AC"))(101)
  positive_palette <- colorRampPalette(c("#BDBDBD", "#B2182B"))(101)
  negative_limit <- max(abs(skewness_values[skewness_values < 0]), na.rm = TRUE)
  positive_limit <- max(skewness_values[skewness_values > 0], na.rm = TRUE)
  if (!is.finite(negative_limit) || negative_limit == 0) negative_limit <- 1
  if (!is.finite(positive_limit) || positive_limit == 0) positive_limit <- 1

  lollipop_colors <- vapply(skewness_values, function(value) {
    if (!is.finite(value) || value == 0) return("#BDBDBD")
    if (value < 0) {
      index <- 1 + round(100 * min(abs(value) / negative_limit, 1))
      return(negative_palette[[index]])
    }
    index <- 1 + round(100 * min(value / positive_limit, 1))
    positive_palette[[index]]
  }, character(1))

  segments(
    x0 = 0,
    y0 = y_positions,
    x1 = absolute_skewness,
    y1 = y_positions,
    col = lollipop_colors,
    lwd = 4.6
  )
  points(
    absolute_skewness,
    y_positions,
    pch = 21,
    cex = 1.20,
    bg = lollipop_colors,
    col = "#202020",
    lwd = 0.95
  )
  axis(1, cex.axis = 0.78, lwd = 0.7, lwd.ticks = 0.7)
  box(col = "grey35", lwd = 0.7)
  mtext("Absolute skewness", side = 1, line = 1.8, cex = 0.88)
  title(
    main = sprintf("Feature Skewness (N=%d)", nrow(metric_df)),
    font.main = 2,
    cex.main = 1.05,
    line = 2
  )
  mtext(
    "Blue: negative   Gray: near zero   Red: positive",
    side = 3,
    line = 0.55,
    cex = 0.65,
    col = "grey20"
  )
}

for (threshold_index in seq_along(connection_thresholds)) {
  connection_threshold <- connection_thresholds[[threshold_index]]
  threshold_tag <- threshold_tags[[threshold_index]]
  file_stem <- sprintf(
    "Feature_Skewness_Correlation_Combined_r_gt_%s",
    threshold_tag
  )
  png_file <- file.path(output_dir, paste0(file_stem, ".png"))
  pdf_file <- file.path(output_dir, paste0(file_stem, ".pdf"))
  connections_file <- file.path(
    output_dir,
    sprintf("Feature_HighCorrelation_Connections_r_gt_%s.csv", threshold_tag)
  )

  connection_index <- which(
    upper.tri(corr_matrix, diag = FALSE) &
      abs(corr_matrix) > connection_threshold,
    arr.ind = TRUE
  )
  connection_table <- data.frame(
    Feature_1 = features[connection_index[, 1]],
    Feature_2 = features[connection_index[, 2]],
    Abbreviation_1 = abbreviations[connection_index[, 1]],
    Abbreviation_2 = abbreviations[connection_index[, 2]],
    Pearson_r = corr_matrix[connection_index],
    Abs_r = abs(corr_matrix[connection_index]),
    stringsAsFactors = FALSE
  )
  connection_table <- connection_table[
    order(connection_table$Abs_r, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  write.csv(connection_table, connections_file, row.names = FALSE)

  plot_combined(png_file, connection_threshold, "png")
  plot_combined(pdf_file, connection_threshold, "pdf")

  cat(sprintf(
    paste0(
      "Completed combined figure: N=%d, metrics=%d, ",
      "|r| > %.2f links=%d.\n"
    ),
    nrow(metric_df), ncol(metric_df),
    connection_threshold, nrow(connection_table)
  ))
  cat(normalizePath(png_file, winslash = "/", mustWork = TRUE), "\n")
  cat(normalizePath(pdf_file, winslash = "/", mustWork = TRUE), "\n")
  cat(normalizePath(connections_file, winslash = "/", mustWork = TRUE), "\n")
}
cat(normalizePath(map_file, winslash = "/", mustWork = TRUE), "\n")
