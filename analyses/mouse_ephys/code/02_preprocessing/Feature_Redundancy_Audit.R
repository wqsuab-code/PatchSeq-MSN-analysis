#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
  library(corrplot)
})

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else
  ".codex-work/e_type_qc/stage1_filtered.json"
output_dir <- if (length(args) >= 2) args[[2]] else
  "outputs/e_type_qc/feature_redundancy_audit"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

output_paths <- c(
  png = file.path(output_dir, "Feature_Redundancy_Audit.png"),
  pdf = file.path(output_dir, "Feature_Redundancy_Audit.pdf"),
  matrix_csv = file.path(output_dir, "Feature_Redundancy_Correlation_Matrix.csv"),
  pairs_csv = file.path(output_dir, "Feature_Redundancy_HighCorrelation_Pairs.csv"),
  log = file.path(output_dir, "Feature_Redundancy_Audit_run_log.txt")
)

log_messages <- character()
log_line <- function(...) {
  msg <- paste0(...)
  log_messages <<- c(log_messages, msg)
  cat(msg, "\n", sep = "")
}

fmt_list <- function(x) {
  if (length(x) == 0) "None" else paste(x, collapse = "; ")
}

# -----------------------------------------------------------------------------
# Read the current cleaned main cohort (549 cells)
# -----------------------------------------------------------------------------
if (!file.exists(input_file)) stop("Input file does not exist: ", input_file)

p <- fromJSON(input_file, simplifyVector = FALSE)
headers <- unlist(p$analysis_headers)
mat <- do.call(rbind, lapply(p$analysis_rows, function(r) {
  vapply(seq_along(headers), function(j) {
    if (is.null(r[[j]])) NA_character_ else as.character(r[[j]])
  }, character(1))
}))
dat <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
names(dat) <- headers
original_n <- nrow(dat)

# E-prefixed columns that are identifiers or vector-valued protocol traces.
non_numeric_or_id_e <- c("E_cellName", "E_current", "E_spike_count")

# Previously defined protocol/metadata fields and QC-rejected features. These
# remain in the source data but are not inputs to the E-type feature audit.
excluded_protocol_or_metadata <- c(
  "E_Vholding", "E_Input.step", "E_Initial.ampltitude"
)
excluded_prior_qc <- c(
  "E_Afterdepolarization..mV.", "E_Burstiness", "E_Wildness",
  "E_Rebound.number.of.APs"
)

all_e_columns <- headers[startsWith(headers, "E_")]
candidate_metrics <- setdiff(
  all_e_columns,
  c(non_numeric_or_id_e, excluded_protocol_or_metadata, excluded_prior_qc)
)
original_candidate_n <- length(candidate_metrics)

# -----------------------------------------------------------------------------
# Force numeric type and remove unusable columns
# -----------------------------------------------------------------------------
conversion_failed <- character()
for (nm in candidate_metrics) {
  original_nonmissing <- !is.na(dat[[nm]]) & nzchar(trimws(dat[[nm]]))
  converted <- suppressWarnings(as.numeric(dat[[nm]]))
  if (any(original_nonmissing & is.na(converted))) {
    conversion_failed <- c(conversion_failed, nm)
  } else {
    dat[[nm]] <- converted
  }
}

candidate_metrics <- setdiff(candidate_metrics, conversion_failed)

all_na_metrics <- candidate_metrics[
  vapply(dat[candidate_metrics], function(x) all(is.na(x)), logical(1))
]
candidate_metrics <- setdiff(candidate_metrics, all_na_metrics)

# -----------------------------------------------------------------------------
# Apply the current six-feature QC definition
# -----------------------------------------------------------------------------
qc_required <- c(
  "E_Holding.MP..mV.", "E_Input.resistance..MOhm.",
  "E_AP.amplitude..mV.", "E_Max.number.of.APs", "E_AP.width..ms.",
  "E_Rheobase..pA."
)
if (!all(qc_required %in% names(dat))) {
  stop("Missing QC columns: ", paste(setdiff(qc_required, names(dat)), collapse = ", "))
}
for (nm in qc_required) dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))

qc_pass <- with(dat,
  E_Holding.MP..mV. < -55 &
  E_Input.resistance..MOhm. >= 100 & E_Input.resistance..MOhm. <= 1000 &
  E_AP.amplitude..mV. > 40 &
  E_Max.number.of.APs > 2 &
  E_AP.width..ms. <= 3 &
  E_Rheobase..pA. >= 10 & E_Rheobase..pA. <= 300
)
qc_pass[is.na(qc_pass)] <- FALSE

metric_df <- dat[qc_pass, candidate_metrics, drop = FALSE]
final_n <- nrow(metric_df)

# -----------------------------------------------------------------------------
# Remove zero-variance and near-zero-variance metrics
# caret-compatible near-zero-variance rule: frequency ratio >= 19 and no more
# than 10% unique values. Zero-variance metrics are reported separately.
# -----------------------------------------------------------------------------
zero_variance <- names(metric_df)[vapply(metric_df, function(x) {
  ux <- unique(x[!is.na(x)])
  length(ux) <= 1
}, logical(1))]
metric_df <- metric_df[setdiff(names(metric_df), zero_variance)]

is_near_zero_variance <- function(x, freq_cut = 19, unique_cut = 10) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(TRUE)
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) <= 1) return(TRUE)
  freq_ratio <- as.numeric(tab[[1]] / tab[[2]])
  pct_unique <- 100 * length(tab) / length(x)
  freq_ratio >= freq_cut && pct_unique <= unique_cut
}

near_zero_variance <- names(metric_df)[vapply(
  metric_df, is_near_zero_variance, logical(1)
)]
metric_df <- metric_df[setdiff(names(metric_df), near_zero_variance)]
final_metric_n <- ncol(metric_df)

if (final_metric_n < 2) stop("Fewer than two usable metrics remain.")

# -----------------------------------------------------------------------------
# Pearson correlation matrix
# -----------------------------------------------------------------------------
corr_matrix <- cor(
  metric_df,
  method = "pearson",
  use = "pairwise.complete.obs"
)
diag(corr_matrix) <- 1

na_idx <- which(is.na(corr_matrix) & lower.tri(corr_matrix), arr.ind = TRUE)
if (nrow(na_idx) > 0) {
  na_pairs <- apply(na_idx, 1, function(z) {
    paste0(rownames(corr_matrix)[z[[1]]], " <-> ", colnames(corr_matrix)[z[[2]]])
  })
  corr_matrix[is.na(corr_matrix)] <- 0
  diag(corr_matrix) <- 1
} else {
  na_pairs <- character()
}

write.csv(corr_matrix, output_paths[["matrix_csv"]], row.names = TRUE)

lower_idx <- which(lower.tri(corr_matrix, diag = FALSE), arr.ind = TRUE)
pair_table <- data.frame(
  Metric_1 = rownames(corr_matrix)[lower_idx[, 1]],
  Metric_2 = colnames(corr_matrix)[lower_idx[, 2]],
  Pearson_r = corr_matrix[lower_idx],
  stringsAsFactors = FALSE
)
pair_table$Abs_r <- abs(pair_table$Pearson_r)
high_pairs <- pair_table[pair_table$Abs_r >= 0.80, , drop = FALSE]
high_pairs <- high_pairs[order(high_pairs$Abs_r, decreasing = TRUE), , drop = FALSE]
row.names(high_pairs) <- NULL
write.csv(high_pairs, output_paths[["pairs_csv"]], row.names = FALSE)

# -----------------------------------------------------------------------------
# Corrplot figure
# -----------------------------------------------------------------------------
plot_ellipse_correlation <- function(
    corr_matrix,
    sample_n,
    output_file,
    file_type = c("png", "pdf")) {
  file_type <- match.arg(file_type)

  if (file_type == "png") {
    png(
      filename = output_file,
      width = 2600,
      height = 2400,
      res = 300,
      bg = "white"
    )
  } else {
    pdf(
      file = output_file,
      width = 8.7,
      height = 8,
      useDingbats = FALSE,
      bg = "white"
    )
  }
  on.exit(dev.off(), add = TRUE)

  correlation_palette <- colorRampPalette(
    c(
      "#8E002B", "#D7301F", "#F4A582", "#F7F7F7",
      "#92C5DE", "#4393C3", "#053061"
    )
  )(201)

  # Shorten labels only for plotting. The source data, exported correlation
  # matrix, and high-correlation table retain the original E_ column names.
  display_corr_matrix <- corr_matrix
  rownames(display_corr_matrix) <- sub("^E_", "", rownames(corr_matrix))
  colnames(display_corr_matrix) <- sub("^E_", "", colnames(corr_matrix))

  par(bg = "white", xpd = NA)
  corr_result <- corrplot::corrplot(
    display_corr_matrix,
    method = "ellipse",
    type = "lower",
    order = "original",
    diag = TRUE,
    col = correlation_palette,
    # "ld" places each column label directly on the corresponding diagonal
    # cell edge; "lt" would put every label on one horizontal top line.
    tl.pos = "ld",
    tl.col = "black",
    tl.cex = 0.58,
    tl.srt = 48,
    tl.offset = 0.35,
    cl.pos = "b",
    col.lim = c(-1, 1),
    cl.length = 11,
    cl.cex = 0.70,
    addgrid.col = "grey82",
    outline = FALSE,
    mar = c(2.5, 1, 4.5, 1)
  )

  # Add centered Pearson r values to every displayed cell, including 1.00 on
  # the diagonal. Black labels match the reference figure and remain legible
  # because strong correlations are represented by narrow ellipses.
  number_pos <- corr_result$corrPos
  text(
    x = number_pos$x,
    y = number_pos$y,
    labels = sprintf("%.2f", number_pos$corr),
    cex = 0.42,
    font = 2,
    col = "black"
  )

  title(
    main = sprintf(
      "Feature Redundancy Audit (N=%d, %d Metrics)",
      sample_n, ncol(corr_matrix)
    ),
    font.main = 2,
    cex.main = 1.05,
    line = 2
  )

  invisible(corr_result)
}

plot_ellipse_correlation(
  corr_matrix = corr_matrix,
  sample_n = nrow(metric_df),
  output_file = output_paths[["png"]],
  file_type = "png"
)

plot_ellipse_correlation(
  corr_matrix = corr_matrix,
  sample_n = nrow(metric_df),
  output_file = output_paths[["pdf"]],
  file_type = "pdf"
)

# -----------------------------------------------------------------------------
# Run log
# -----------------------------------------------------------------------------
log_line("Feature Redundancy Audit")
log_line("Run time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
log_line("Input file: ", normalizePath(input_file, winslash = "/", mustWork = TRUE))
log_line("Original sample count: ", original_n)
log_line("Final included sample count (QC Pass): ", final_n)
log_line("All E-prefixed source columns: ", length(all_e_columns))
log_line("Original candidate metric count: ", original_candidate_n)
log_line("Removed non-numeric/ID/vector E columns: ", fmt_list(non_numeric_or_id_e[non_numeric_or_id_e %in% all_e_columns]))
log_line("Removed protocol/metadata E columns: ", fmt_list(excluded_protocol_or_metadata[excluded_protocol_or_metadata %in% all_e_columns]))
log_line("Previously QC-excluded E metrics: ", fmt_list(excluded_prior_qc[excluded_prior_qc %in% all_e_columns]))
log_line("Failed numeric conversion: ", fmt_list(conversion_failed))
log_line("Removed all-NA metrics: ", fmt_list(all_na_metrics))
log_line("Removed zero-variance metrics: ", fmt_list(zero_variance))
log_line("Removed near-zero-variance metrics: ", fmt_list(near_zero_variance))
log_line("Final metric count: ", final_metric_n)
log_line("Residual NA correlation pairs replaced by 0: ", fmt_list(na_pairs))
log_line("|r| >= 0.80 metric-pair count: ", nrow(high_pairs))
log_line("Output files:")
for (path in output_paths) {
  log_line("  ", normalizePath(path, winslash = "/", mustWork = FALSE))
}

writeLines(log_messages, output_paths[["log"]], useBytes = TRUE)

cat("Completed Feature Redundancy Audit.\n")
