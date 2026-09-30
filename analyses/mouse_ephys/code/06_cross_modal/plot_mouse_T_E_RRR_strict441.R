#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(clue)
})

set.seed(777)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript so the module root can be resolved.")
script_file <- sub("^--file=", "", script_arg[[1]])
root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
default_out <- file.path(root, "results", "04_T_E_RRR", "generated")
out <- Sys.getenv("MOUSE_E_RRR_OUT", unset = default_out)
source_site_data <- file.path(root, "interactive", "data")
site_data <- Sys.getenv("MOUSE_E_RRR_SITE_OUT", unset = source_site_data)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(site_data, recursive = TRUE, showWarnings = FALSE)

stability_file <- file.path(root, "data", "02_frozen_input", "e_stability_all493.csv")
e_file <- file.path(root, "data", "02_frozen_input", "NPC3_HC5_Raw_Final18_Features.csv")
rna_rds <- Sys.getenv("MOUSE_E_RNA_RDS", unset = "")
if (!nzchar(rna_rds) || !file.exists(rna_rds)) {
  stop("Set MOUSE_E_RNA_RDS to the external RNA RDS listed in data/00_source_manifest/release_assets.csv")
}

e_names <- c(
  "E_Holding.MP..mV." = "Holding MP", "E_Input.resistance..MOhm." = "Input resistance",
  "E_Membrane.time.constant..ms." = "Membrane tau", "E_Rheobase..pA." = "Rheobase",
  "E_Sag.ratio" = "Sag ratio", "E_Sag.time..s." = "Sag time",
  "E_AP.threshold..mV." = "AP threshold", "E_AP.amplitude..mV." = "AP amplitude",
  "E_AP.width..ms." = "AP width", "E_Upstroke.to.downstroke.ratio" = "Up/down ratio",
  "E_Afterhyperpolarization..mV." = "AHP", "E_Max.number.of.APs" = "Max APs",
  "E_Latency..ms." = "Latency", "E_Latency....20pA.current..ms." = "Latency (+20 pA)",
  "E_ISI.adaptation.index" = "ISI adaptation", "E_ISI.coefficient.of.variation" = "ISI CV",
  "E_AP.amplitude.adaptation.index" = "AP amp adaptation",
  "E_AP.coefficient.of.variation" = "AP CV"
)
d_colors <- c(D1 = "#ff379b", D2 = "#00eeb3")
e_colors <- c(E1 = "#4E79A7", E2 = "#F28E2B", E3 = "#59A14F", E4 = "#2AA6B8", E5 = "#B07AA1")

truthy <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")

stability <- fread(stability_file)
strict <- stability[
  truthy(HC_GC_consensus) &
    truthy(major_class_agreement) &
    T_identity %in% c("D1", "D2"),
  .(MSN_unique_ID, D1_D2 = T_identity, E_class)
]
e <- fread(e_file)
dat <- merge(strict, e, by = "MSN_unique_ID", sort = FALSE)
setorderv(dat, "MSN_unique_ID")
stopifnot(
  nrow(dat) == 441L,
  dat[D1_D2 == "D1", .N] == 209L,
  dat[D1_D2 == "D2", .N] == 232L,
  all(dat$E_class %in% names(e_colors))
)

# Transcriptomic input: library-size 10,000, log1p, 1,000 nontechnical HVGs,
# gene z-score, then 20 PCs, all refit in the 441 strict T-stable cells.
raw <- readRDS(rna_rds)
counts <- raw$query_counts[, dat$MSN_unique_ID, drop = FALSE]
stopifnot(inherits(counts, "Matrix"), all(Matrix::colSums(counts) > 0))
log_norm <- log1p(counts %*% Diagonal(x = 10000 / Matrix::colSums(counts)))
gene_mean <- Matrix::rowMeans(log_norm)
gene_var <- pmax(Matrix::rowMeans(log_norm ^ 2) - gene_mean ^ 2, 0)
genes <- rownames(log_norm)
technical <- grepl("^(mt-|Mt-|Rpl|Rps|Gm[0-9]|[0-9].*Rik$|Malat1$|Xist$)", genes)
eligible <- which(!technical & is.finite(gene_var) & gene_var > 0)
eligible <- eligible[order(gene_var[eligible], decreasing = TRUE)][seq_len(1000)]
gene_z <- t(scale(t(as.matrix(log_norm[eligible, , drop = FALSE])), center = TRUE, scale = TRUE))
gene_z <- gene_z[apply(gene_z, 1, function(x) all(is.finite(x))), , drop = FALSE]
rna_pca <- prcomp(t(gene_z), center = FALSE, scale. = FALSE, rank. = 20)
X <- rna_pca$x[, 1:20, drop = FALSE]

# Frozen 18-feature E panel; transformation refit within the 441-cell cohort.
features <- names(e_names)
Y0 <- as.matrix(dat[, ..features])
mins <- apply(Y0, 2, min)
shifted <- sweep(Y0, 2, mins, "-")
norm <- sweep(shifted, 2, colSums(shifted), "/") * 10000
Y <- scale(log1p(norm), center = TRUE, scale = TRUE)

# Rank-3 reduced-rank regression, T -> E.
Xc <- scale(X, center = TRUE, scale = FALSE)
Yc <- scale(Y, center = TRUE, scale = FALSE)
B <- qr.solve(Xc, Yc)
fit <- Xc %*% B
sv <- svd(fit, nu = 0, nv = 3)
V <- sv$v[, 1:3, drop = FALSE]
for (k in 1:3) {
  anchor <- which.max(abs(V[, k]))
  if (V[anchor, k] < 0) V[, k] <- -V[, k]
}
Tproj <- fit %*% V
Eproj <- Yc %*% V
Tload <- cor(t(gene_z), Tproj)
Eload <- cor(Y, Eproj)
colnames(Tproj) <- colnames(Eproj) <- colnames(Tload) <- colnames(Eload) <- paste0("Component", 1:3)
rownames(Tproj) <- rownames(Eproj) <- dat$MSN_unique_ID

scores <- data.table(
  MSN_unique_ID = dat$MSN_unique_ID,
  D1_D2 = dat$D1_D2,
  E_class = dat$E_class,
  T_Component1 = Tproj[, 1], T_Component2 = Tproj[, 2], T_Component3 = Tproj[, 3],
  E_Component1 = Eproj[, 1], E_Component2 = Eproj[, 2], E_Component3 = Eproj[, 3]
)
genes_out <- data.table(Gene = rownames(Tload), Tload)
features_out <- data.table(Feature = features, Display = unname(e_names), Eload)
pca_out <- data.table(PC = 1:20, variance_fraction = rna_pca$sdev[1:20]^2 / sum(rna_pca$sdev^2))
counts_out <- as.data.table(table(dat$D1_D2, dat$E_class))
r2 <- 1 - sum((Yc - fit %*% V %*% t(V))^2) / sum(Yc^2)

fwrite(scores, file.path(out, "RRR_cell_scores_n441.csv"))
fwrite(genes_out, file.path(out, "RRR_gene_correlation_loadings_n441.csv"))
fwrite(features_out, file.path(out, "RRR_Efeature_correlation_loadings_n441.csv"))
fwrite(pca_out, file.path(out, "transcriptomic_PCA_variance_n441.csv"))
fwrite(counts_out, file.path(out, "D1D2_by_Eclass_counts_n441.csv"))

writeLines(c(
  "Mouse NAc transcriptomic-to-electrophysiology reduced-rank regression.",
  "Cohort: 441 strict T-stable GC-HC E-class consensus cells (D1=209, D2=232).",
  "T: library-size normalization to 10,000, log1p, 1,000 nontechnical HVGs, gene z-score, 20 PCs.",
  "E: 18 frozen core E-features, minimum shift, column-sum normalization to 10,000, log1p, feature z-score.",
  "RRR: pooled T-to-E model, rank=3; refit directly in the strict 441-cell cohort.",
  sprintf("Descriptive in-sample multivariate R2: %.4f.", r2),
  "D1/D2 ellipses in the interactive editor are descriptive 90% bivariate-normal regions.",
  "Display layout is controlled by Mouse_T-E_RRR_layout_09172026.json."
), file.path(out, "README_strict441.txt"))

# Update interactive editor data to the newly refit strict-441 RRR tables.
fwrite(scores, file.path(site_data, "rrr_scores.csv"))
fwrite(genes_out, file.path(site_data, "rrr_gene_loadings.csv"))
fwrite(features_out, file.path(site_data, "rrr_feature_loadings.csv"))

# Preserve explicit GC-merged names for the interactive two-version selector.
fwrite(scores, file.path(site_data, "rrr_scores_gc_merged.csv"))
fwrite(genes_out, file.path(site_data, "rrr_gene_loadings_gc_merged.csv"))
fwrite(features_out, file.path(site_data, "rrr_feature_loadings_gc_merged.csv"))

# Exploratory GC-unmerged analysis. The fixed raw graph clustering has K=11;
# its best independently scanned Ward.D2 comparison is HC K=13. Raw GC labels
# are matched one-to-one to HC labels by the Hungarian solution, and only
# matched cells with a strict D1/D2 transcriptomic identity enter the refit.
raw_assignment_file <- file.path(root, "legacy", "sensitivity", "GCrawK11_HCK13_consensus_audit.csv")
raw_assignment <- fread(raw_assignment_file)

strict_t <- stability[
  truthy(major_class_agreement) & T_identity %in% c("D1", "D2"),
  .(MSN_unique_ID, D1_D2 = T_identity)
]
strict_raw <- merge(
  raw_assignment[raw_GC_HC_consensus == TRUE], strict_t,
  by = "MSN_unique_ID", sort = FALSE
)
strict_raw[, E_class := sub("^S", "GC", Seurat_raw_cluster)]
raw_dat <- merge(strict_raw, e, by = "MSN_unique_ID", sort = FALSE)
setorderv(raw_dat, "MSN_unique_ID")

fit_rrr_version <- function(version_dat) {
  version_counts <- raw$query_counts[, version_dat$MSN_unique_ID, drop = FALSE]
  version_log_norm <- log1p(version_counts %*% Diagonal(x = 10000 / Matrix::colSums(version_counts)))
  version_gene_mean <- Matrix::rowMeans(version_log_norm)
  version_gene_var <- pmax(Matrix::rowMeans(version_log_norm ^ 2) - version_gene_mean ^ 2, 0)
  version_genes <- rownames(version_log_norm)
  version_technical <- grepl("^(mt-|Mt-|Rpl|Rps|Gm[0-9]|[0-9].*Rik$|Malat1$|Xist$)", version_genes)
  version_eligible <- which(!version_technical & is.finite(version_gene_var) & version_gene_var > 0)
  version_eligible <- version_eligible[order(version_gene_var[version_eligible], decreasing = TRUE)][seq_len(1000)]
  version_gene_z <- t(scale(t(as.matrix(version_log_norm[version_eligible, , drop = FALSE])), center = TRUE, scale = TRUE))
  version_gene_z <- version_gene_z[apply(version_gene_z, 1, function(x) all(is.finite(x))), , drop = FALSE]
  version_pca <- prcomp(t(version_gene_z), center = FALSE, scale. = FALSE, rank. = 20)
  version_x <- version_pca$x[, 1:20, drop = FALSE]

  version_y0 <- as.matrix(version_dat[, ..features])
  version_mins <- apply(version_y0, 2, min)
  version_shifted <- sweep(version_y0, 2, version_mins, "-")
  version_norm <- sweep(version_shifted, 2, colSums(version_shifted), "/") * 10000
  version_y <- scale(log1p(version_norm), center = TRUE, scale = TRUE)
  version_xc <- scale(version_x, center = TRUE, scale = FALSE)
  version_yc <- scale(version_y, center = TRUE, scale = FALSE)
  version_b <- qr.solve(version_xc, version_yc)
  version_fit <- version_xc %*% version_b
  version_sv <- svd(version_fit, nu = 0, nv = 3)
  version_v <- version_sv$v[, 1:3, drop = FALSE]
  for (component in 1:3) {
    anchor <- which.max(abs(version_v[, component]))
    if (version_v[anchor, component] < 0) version_v[, component] <- -version_v[, component]
  }
  version_tproj <- version_fit %*% version_v
  version_eproj <- version_yc %*% version_v
  version_tload <- cor(t(version_gene_z), version_tproj)
  version_eload <- cor(version_y, version_eproj)
  colnames(version_tproj) <- colnames(version_eproj) <- colnames(version_tload) <- colnames(version_eload) <- paste0("Component", 1:3)
  list(
    scores = data.table(
      MSN_unique_ID = version_dat$MSN_unique_ID,
      D1_D2 = version_dat$D1_D2,
      E_class = version_dat$E_class,
      T_Component1 = version_tproj[, 1], T_Component2 = version_tproj[, 2], T_Component3 = version_tproj[, 3],
      E_Component1 = version_eproj[, 1], E_Component2 = version_eproj[, 2], E_Component3 = version_eproj[, 3]
    ),
    genes = data.table(Gene = rownames(version_tload), version_tload),
    efeatures = data.table(Feature = features, Display = unname(e_names), version_eload),
    r2 = 1 - sum((version_yc - version_fit %*% version_v %*% t(version_v))^2) / sum(version_yc^2)
  )
}

raw_fit <- fit_rrr_version(raw_dat)
fwrite(raw_fit$scores, file.path(site_data, "rrr_scores_gc_unmerged.csv"))
fwrite(raw_fit$genes, file.path(site_data, "rrr_gene_loadings_gc_unmerged.csv"))
fwrite(raw_fit$efeatures, file.path(site_data, "rrr_feature_loadings_gc_unmerged.csv"))
fwrite(raw_fit$scores, file.path(out, "RRR_cell_scores_GCrawK11_HCK13.csv"))
fwrite(raw_fit$genes, file.path(out, "RRR_gene_correlation_loadings_GCrawK11_HCK13.csv"))
fwrite(raw_fit$efeatures, file.path(out, "RRR_Efeature_correlation_loadings_GCrawK11_HCK13.csv"))
fwrite(raw_assignment, file.path(out, "GCrawK11_HCK13_consensus_audit.csv"))

version_manifest <- data.table(
  version_key = c("gc_merged", "gc_unmerged"),
  label = c("GC merged (K=5; primary)", "GC unmerged (raw K=11; exploratory)"),
  GC_K = c(5L, 11L), HC_K = c(5L, 13L),
  agreement_n = c(450L, sum(raw_assignment$raw_GC_HC_consensus)),
  agreement_denominator = c(493L, nrow(raw_assignment)),
  strict_rrr_n = c(nrow(dat), nrow(raw_dat)),
  D1_n = c(dat[D1_D2 == "D1", .N], raw_dat[D1_D2 == "D1", .N]),
  D2_n = c(dat[D1_D2 == "D2", .N], raw_dat[D1_D2 == "D2", .N]),
  descriptive_R2 = c(r2, raw_fit$r2),
  class_field = c("E1-E5", "GC1-GC11")
)
fwrite(version_manifest, file.path(site_data, "rrr_gc_versions_manifest.csv"))
fwrite(version_manifest, file.path(out, "RRR_GC_versions_manifest.csv"))

cat("strict441 T->E RRR complete\n")
cat("out:", out, "\n")
cat(sprintf("in-sample R2: %.4f\n", r2))
print(table(dat$D1_D2, dat$E_class))
cat(sprintf("GC-unmerged strict cohort: n=%d; D1=%d; D2=%d; in-sample R2=%.4f\n",
            nrow(raw_dat), raw_dat[D1_D2 == "D1", .N], raw_dat[D1_D2 == "D2", .N], raw_fit$r2))
print(table(raw_dat$D1_D2, raw_dat$E_class))
