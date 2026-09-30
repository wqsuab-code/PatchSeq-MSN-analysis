#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

input_file <- file.path(".codex-work", "e_type_qc", "stage1_filtered.json")
output_dir <- file.path("outputs", "e_type_qc", "parallel_preprocessing")
yj_dir <- file.path(output_dir, "01_yeo_johnson")
sum_dir <- file.path(output_dir, "02_shift_sum_norm_10000")
dir.create(yj_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sum_dir, recursive = TRUE, showWarnings = FALSE)

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

payload <- fromJSON(input_file, simplifyVector = FALSE)
headers <- unlist(payload$analysis_headers)
matrix_data <- do.call(rbind, lapply(payload$analysis_rows, function(row) {
  vapply(seq_along(headers), function(j) {
    if (is.null(row[[j]])) NA_character_ else as.character(row[[j]])
  }, character(1))
}))
dat <- as.data.frame(matrix_data, stringsAsFactors = FALSE, check.names = FALSE)
names(dat) <- headers

missing_features <- setdiff(features, names(dat))
if (length(missing_features) > 0) {
  stop("Missing requested features: ", paste(missing_features, collapse = ", "))
}
for (feature in features) {
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

qc_dat <- dat[qc_pass, , drop = FALSE]
raw_df <- qc_dat[, features, drop = FALSE]
raw_matrix <- as.matrix(raw_df)
storage.mode(raw_matrix) <- "double"

if (nrow(raw_matrix) != 494) {
  stop("Expected 494 QC-pass cells but found ", nrow(raw_matrix), ".")
}
if (anyNA(raw_matrix) || any(!is.finite(raw_matrix))) {
  stop("The 494-cell PCA input contains missing or non-finite values.")
}
if (anyDuplicated(qc_dat$MSN_unique_ID)) {
  stop("MSN_unique_ID is not unique in the QC-pass cohort.")
}
rownames(raw_matrix) <- qc_dat$MSN_unique_ID

feature_sd <- apply(raw_matrix, 2, sd)
if (any(!is.finite(feature_sd) | feature_sd == 0)) {
  stop(
    "Zero-variance features detected: ",
    paste(names(feature_sd)[!is.finite(feature_sd) | feature_sd == 0], collapse = ", ")
  )
}

metadata_columns <- setdiff(names(qc_dat), features)
metadata_df <- qc_dat[, metadata_columns, drop = FALSE]
metadata_df$Ephys_QC_status <- "Pass"

id_table <- data.frame(
  MSN_unique_ID = qc_dat$MSN_unique_ID,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write.csv(
  metadata_df,
  file.path(output_dir, "Ephys_QCpass_494_metadata.csv"),
  row.names = FALSE
)
write.csv(
  cbind(id_table, as.data.frame(raw_matrix, check.names = FALSE)),
  file.path(output_dir, "Ephys_QCpass_494_raw_25features.csv"),
  row.names = FALSE
)

# Branch 1: bounded, feature-wise Yeo-Johnson maximum-likelihood fitting,
# followed by Z-scoring. Boundary optima are capped at +/-5 instead of silently
# reverting to the identity transform; this guarantees that all 25 metrics are
# processed while avoiding numerically extreme lambda estimates.
yj_transform <- function(x, lambda, eps = 1e-8) {
  transformed <- numeric(length(x))
  nonnegative <- x >= 0
  if (abs(lambda) < eps) {
    transformed[nonnegative] <- log1p(x[nonnegative])
  } else {
    transformed[nonnegative] <-
      expm1(lambda * log1p(x[nonnegative])) / lambda
  }
  if (abs(lambda - 2) < eps) {
    transformed[!nonnegative] <- -log1p(-x[!nonnegative])
  } else {
    transformed[!nonnegative] <-
      -expm1((2 - lambda) * log1p(-x[!nonnegative])) / (2 - lambda)
  }
  transformed
}

yj_log_likelihood <- function(lambda, x) {
  transformed <- yj_transform(x, lambda)
  variance_mle <- mean((transformed - mean(transformed))^2)
  if (!is.finite(variance_mle) || variance_mle <= 0) return(-Inf)
  jacobian_constant <- sum(sign(x) * log1p(abs(x)))
  -0.5 * length(x) * log(variance_mle) +
    (lambda - 1) * jacobian_constant
}

yj_limits <- c(-5, 5)
yj_fit <- lapply(features, function(feature) {
  values <- raw_matrix[, feature]
  fit <- optimize(
    yj_log_likelihood,
    interval = yj_limits,
    maximum = TRUE,
    x = values,
    tol = 1e-6
  )
  lambda <- fit$maximum
  boundary_tolerance <- 0.001
  boundary_capped <-
    abs(lambda - yj_limits[1]) <= boundary_tolerance ||
    abs(lambda - yj_limits[2]) <= boundary_tolerance
  list(
    lambda = lambda,
    boundary_capped = boundary_capped,
    log_likelihood = fit$objective
  )
})
names(yj_fit) <- features

yj_lambda_table <- data.frame(
  Feature = features,
  Lambda = vapply(yj_fit, `[[`, numeric(1), "lambda"),
  Boundary_capped = vapply(yj_fit, `[[`, logical(1), "boundary_capped"),
  Search_lower = yj_limits[1],
  Search_upper = yj_limits[2],
  Log_likelihood = vapply(yj_fit, `[[`, numeric(1), "log_likelihood"),
  stringsAsFactors = FALSE
)

yj_matrix <- vapply(features, function(feature) {
  yj_transform(raw_matrix[, feature], yj_fit[[feature]]$lambda)
}, numeric(nrow(raw_matrix)))
colnames(yj_matrix) <- features
rownames(yj_matrix) <- qc_dat$MSN_unique_ID
yj_z_matrix <- scale(yj_matrix, center = TRUE, scale = TRUE)

write.csv(
  cbind(id_table, as.data.frame(yj_matrix, check.names = FALSE)),
  file.path(yj_dir, "Ephys_YeoJohnson_analysis_table.csv"),
  row.names = FALSE
)
write.csv(
  cbind(id_table, as.data.frame(yj_z_matrix, check.names = FALSE)),
  file.path(yj_dir, "Ephys_YeoJohnson_Zscore_PCA_input.csv"),
  row.names = FALSE
)
write.csv(
  cbind(metadata_df, as.data.frame(yj_matrix, check.names = FALSE)),
  file.path(yj_dir, "Ephys_YeoJohnson_with_metadata.csv"),
  row.names = FALSE
)
write.csv(
  yj_lambda_table,
  file.path(yj_dir, "Ephys_YeoJohnson_lambda.csv"),
  row.names = FALSE
)

# Branch 2: exact requested translation and column-sum normalization.
feature_minimum <- apply(raw_matrix, 2, min)
shift_matrix <- sweep(raw_matrix, 2, feature_minimum, FUN = "-")
shifted_column_sum <- colSums(shift_matrix)
if (any(!is.finite(shifted_column_sum) | shifted_column_sum <= 0)) {
  stop(
    "Cannot sum-normalize features: ",
    paste(
      names(shifted_column_sum)[
        !is.finite(shifted_column_sum) | shifted_column_sum <= 0
      ],
      collapse = ", "
    )
  )
}
sum_scale_factor <- 10000 / shifted_column_sum
sum_norm_matrix <- sweep(shift_matrix, 2, shifted_column_sum, FUN = "/") * 10000
sum_log1p_matrix <- log1p(sum_norm_matrix)
sum_log1p_z_matrix <- scale(sum_log1p_matrix, center = TRUE, scale = TRUE)

sum_parameter_table <- data.frame(
  Feature = features,
  Original_minimum = feature_minimum[features],
  Shifted_column_sum = shifted_column_sum[features],
  Scale_factor_10000_over_sum = sum_scale_factor[features],
  Final_column_sum = colSums(sum_norm_matrix)[features],
  stringsAsFactors = FALSE
)

write.csv(
  cbind(id_table, as.data.frame(shift_matrix, check.names = FALSE)),
  file.path(sum_dir, "Ephys_Shifted_to_zero_analysis_table.csv"),
  row.names = FALSE
)
write.csv(
  cbind(id_table, as.data.frame(sum_norm_matrix, check.names = FALSE)),
  file.path(sum_dir, "Ephys_ShiftSumNorm10000_intermediate.csv"),
  row.names = FALSE
)
write.csv(
  cbind(id_table, as.data.frame(sum_log1p_matrix, check.names = FALSE)),
  file.path(sum_dir, "Ephys_ShiftSumNorm10000_Log1p_intermediate.csv"),
  row.names = FALSE
)
write.csv(
  cbind(id_table, as.data.frame(sum_log1p_z_matrix, check.names = FALSE)),
  file.path(sum_dir, "Ephys_ShiftSumNorm10000_Log1p_Zscore_PCA_input.csv"),
  row.names = FALSE
)
# Keep the generic branch filename synchronized with the corrected final PCA
# input so earlier links cannot silently point to the incomplete two-step data.
write.csv(
  cbind(id_table, as.data.frame(sum_log1p_z_matrix, check.names = FALSE)),
  file.path(sum_dir, "Ephys_ShiftSumNorm10000_PCA_input.csv"),
  row.names = FALSE
)
write.csv(
  cbind(metadata_df, as.data.frame(sum_log1p_z_matrix, check.names = FALSE)),
  file.path(sum_dir, "Ephys_ShiftSumNorm10000_with_metadata.csv"),
  row.names = FALSE
)
write.csv(
  sum_parameter_table,
  file.path(sum_dir, "Ephys_ShiftSumNorm10000_parameters.csv"),
  row.names = FALSE
)

moment_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3 || sd(x) == 0) return(NA_real_)
  mean(((x - mean(x)) / sd(x))^3)
}

diagnostics <- data.frame(
  Feature = features,
  Raw_min = apply(raw_matrix, 2, min)[features],
  Raw_max = apply(raw_matrix, 2, max)[features],
  Raw_mean = colMeans(raw_matrix)[features],
  Raw_SD = apply(raw_matrix, 2, sd)[features],
  Raw_skewness = apply(raw_matrix, 2, moment_skewness)[features],
  YJ_lambda = yj_lambda_table$Lambda[match(features, yj_lambda_table$Feature)],
  YJ_mean = colMeans(yj_matrix)[features],
  YJ_SD = apply(yj_matrix, 2, sd)[features],
  YJ_skewness = apply(yj_matrix, 2, moment_skewness)[features],
  ShiftSum_min = apply(sum_norm_matrix, 2, min)[features],
  ShiftSum_max = apply(sum_norm_matrix, 2, max)[features],
  ShiftSum_mean = colMeans(sum_norm_matrix)[features],
  ShiftSum_SD = apply(sum_norm_matrix, 2, sd)[features],
  ShiftSum_skewness = apply(sum_norm_matrix, 2, moment_skewness)[features],
  ShiftSum_Log1p_mean = colMeans(sum_log1p_matrix)[features],
  ShiftSum_Log1p_SD = apply(sum_log1p_matrix, 2, sd)[features],
  ShiftSum_Log1p_skewness = apply(
    sum_log1p_matrix, 2, moment_skewness
  )[features],
  stringsAsFactors = FALSE
)
write.csv(
  diagnostics,
  file.path(output_dir, "Ephys_Parallel_Preprocessing_Feature_Diagnostics.csv"),
  row.names = FALSE
)

raw_z_matrix <- scale(raw_matrix, center = TRUE, scale = TRUE)
linear_stage_z_matrix <- scale(sum_norm_matrix, center = TRUE, scale = TRUE)
max_linear_zscore_difference <- max(abs(raw_z_matrix - linear_stage_z_matrix))
raw_corr <- cor(raw_matrix, method = "pearson")
sum_corr <- cor(sum_norm_matrix, method = "pearson")
log_corr <- cor(sum_log1p_matrix, method = "pearson")
max_corr_difference <- max(abs(raw_corr - sum_corr))
max_log_corr_change <- max(abs(raw_corr - log_corr))
max_skewness_difference <- max(abs(
  diagnostics$Raw_skewness - diagnostics$ShiftSum_skewness
))
max_log_skewness_change <- max(abs(
  diagnostics$Raw_skewness - diagnostics$ShiftSum_Log1p_skewness
))
column_sum_error <- max(abs(colSums(sum_norm_matrix) - 10000))
final_z_mean_error <- max(abs(colMeans(sum_log1p_z_matrix)))
final_z_sd_error <- max(abs(apply(sum_log1p_z_matrix, 2, sd) - 1))

bundle_file <- file.path(output_dir, "Ephys_Parallel_Preprocessing_bundle.rds")
saveRDS(
  list(
    cell_id = qc_dat$MSN_unique_ID,
    feature_order = features,
    metadata = metadata_df,
    raw = raw_matrix,
    yeo_johnson = yj_matrix,
    yeo_johnson_zscore = yj_z_matrix,
    yeo_johnson_lambda = yj_lambda_table,
    shifted_to_zero = shift_matrix,
    shift_sum_norm_10000 = sum_norm_matrix,
    shift_sum_norm_10000_log1p = sum_log1p_matrix,
    shift_sum_norm_10000_log1p_zscore = sum_log1p_z_matrix,
    shift_sum_norm_parameters = sum_parameter_table,
    diagnostics = diagnostics
  ),
  bundle_file
)

verification_file <- file.path(
  output_dir,
  "Ephys_Parallel_Preprocessing_verification.txt"
)
verification_lines <- c(
  sprintf("Input cells before six-feature QC: %d", nrow(dat)),
  sprintf("QC-pass cells retained: %d", nrow(raw_matrix)),
  sprintf("Numeric E-features: %d", ncol(raw_matrix)),
  sprintf("Missing/non-finite values in PCA input: %d", sum(!is.finite(raw_matrix))),
  sprintf("Unique MSN_unique_ID values: %d", length(unique(qc_dat$MSN_unique_ID))),
  "",
  paste0(
    "Branch 1: bounded per-feature Yeo-Johnson maximum-likelihood ",
    "transformation (lambda in [-5, 5]) followed by Z-score."
  ),
  sprintf(
    "Features with boundary-capped lambda: %d of %d",
    sum(yj_lambda_table$Boundary_capped),
    nrow(yj_lambda_table)
  ),
  sprintf(
    "Median absolute skewness: raw=%.6f; Yeo-Johnson=%.6f",
    median(abs(diagnostics$Raw_skewness)),
    median(abs(diagnostics$YJ_skewness))
  ),
  "",
  paste0(
    "Branch 2: per-feature shift to zero, divide by shifted column sum, ",
    "multiply by 10000, apply log1p, then Z-score each feature."
  ),
  sprintf("Maximum absolute final column-sum error: %.12g", column_sum_error),
  sprintf(
    "Maximum raw vs pre-log ShiftSum Pearson-correlation difference: %.12g",
    max_corr_difference
  ),
  sprintf(
    "Maximum raw vs pre-log ShiftSum skewness difference: %.12g",
    max_skewness_difference
  ),
  sprintf(
    "Maximum raw-Zscore vs pre-log ShiftSum-Zscore difference: %.12g",
    max_linear_zscore_difference
  ),
  sprintf(
    "Maximum raw vs post-log Pearson-correlation change: %.12g",
    max_log_corr_change
  ),
  sprintf(
    "Maximum raw vs post-log skewness change: %.12g",
    max_log_skewness_change
  ),
  sprintf(
    "Final log1p-Zscore maximum absolute column mean: %.12g",
    final_z_mean_error
  ),
  sprintf(
    "Final log1p-Zscore maximum absolute SD-minus-one: %.12g",
    final_z_sd_error
  ),
  "",
  paste0(
    "Interpretation: the shift and sum-normalization stages alone are positive ",
    "affine transformations and preserve skewness and Pearson correlation. ",
    "The subsequent log1p stage is nonlinear and changes both distributional ",
    "shape and correlation structure; feature-wise Z-scoring is therefore ",
    "required before PCA."
  ),
  paste0(
    "For parallel PCA, use Ephys_YeoJohnson_Zscore_PCA_input.csv for ",
    "Branch 1 and Ephys_ShiftSumNorm10000_Log1p_Zscore_PCA_input.csv for ",
    "the article-matched Branch 2."
  )
)
writeLines(verification_lines, verification_file, useBytes = TRUE)

cat(paste(verification_lines, collapse = "\n"), "\n")
cat("\nOutputs:\n")
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
cat(normalizePath(bundle_file, winslash = "/", mustWork = TRUE), "\n")
cat(normalizePath(verification_file, winslash = "/", mustWork = TRUE), "\n")
