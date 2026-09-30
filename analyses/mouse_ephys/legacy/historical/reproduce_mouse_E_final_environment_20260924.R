#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Rtsne)
  library(mclust)
})

options(stringsAsFactors = FALSE)
set.seed(20260722L)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
bundle <- file.path(root, "outputs", "Mouse_E_frozen_analysis_bundle_20260923")
out <- file.path(root, "outputs", "Mouse_E_final_performed_using_20260924", "01_core_taxonomy_R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

input_dir <- file.path(bundle, "01_frozen_inputs", "electrophysiology")
raw_file <- file.path(input_dir, "Ephys_QCpass_494_raw_25features.csv")
active_z_file <- file.path(input_dir, "Active_Final18_PCA_input.csv")
active_pc_file <- file.path(input_dir, "NPC3_Frozen_PC_Scores_HC5.csv")
active_hc_file <- file.path(input_dir, "Active_Frozen_HC5_assignments.csv")
source_file <- file.path(input_dir, "Active_Frozen_source_assignments.csv")
class_file <- file.path(bundle, "02_frozen_analysis", "classification", "HC_GC_only_cell_assignments.csv")
working_rds <- file.path(
  root, "outputs", "e_type_qc", "fixed_npc3_res1.5_merged_k5_tsne",
  "Fixed_NPC3_Res1.5_Merged_K5_tSNE_Working_Object.rds"
)

raw <- read.csv(raw_file, check.names = FALSE)
active_z <- read.csv(active_z_file, check.names = FALSE)
frozen_pc <- read.csv(active_pc_file, check.names = FALSE)
frozen_hc <- read.csv(active_hc_file, check.names = FALSE)
source <- read.csv(source_file, check.names = FALSE)
frozen_class <- read.csv(class_file, check.names = FALSE)
id <- "MSN_unique_ID"
features <- setdiff(names(active_z), id)

stopifnot(nrow(raw) == 494L, nrow(active_z) == 493L, length(features) == 18L)
stopifnot(all(features %in% names(raw)), setequal(active_z[[id]], source[[id]]))

# Refit the frozen primary transformation in the original 494-cell complete-case cohort.
x0 <- as.matrix(raw[, features, drop = FALSE])
storage.mode(x0) <- "double"
mins <- apply(x0, 2, min)
shifted <- sweep(x0, 2, mins, "-")
normalized <- sweep(shifted, 2, colSums(shifted), "/") * 10000
logged <- log1p(normalized)
z494 <- scale(logged, center = TRUE, scale = TRUE)
rownames(z494) <- raw[[id]]

z_active <- z494[active_z[[id]], , drop = FALSE]
frozen_z <- as.matrix(active_z[, features, drop = FALSE])
max_z_diff <- max(abs(z_active - frozen_z))

# Refit PCA, Ward.D2 HC and SNN/Louvain GC on the same 494-cell matrix.
pca_fit <- prcomp(z494, center = FALSE, scale. = FALSE)
pcs <- pca_fit$x[, 1:3, drop = FALSE]
pc_active <- pcs[frozen_pc[[id]], , drop = FALSE]
frozen_pc_matrix <- as.matrix(frozen_pc[, c("PC1", "PC2", "PC3")])
for (j in 1:3) {
  if (cor(pc_active[, j], frozen_pc_matrix[, j]) < 0) {
    pcs[, j] <- -pcs[, j]
    pca_fit$rotation[, j] <- -pca_fit$rotation[, j]
  }
}
pc_active <- pcs[frozen_pc[[id]], , drop = FALSE]
max_pc_diff <- max(abs(pc_active - frozen_pc_matrix))

hc5 <- cutree(hclust(dist(pcs), method = "ward.D2"), k = 5L)
names(hc5) <- rownames(pcs)
hc_active <- paste0("HC", hc5[frozen_hc[[id]]])
hc_exact <- identical(unname(hc_active), as.character(frozen_hc$HC5))

snn <- Seurat::FindNeighbors(
  pcs, k.param = 20L, compute.SNN = TRUE, prune.SNN = 1 / 15,
  verbose = FALSE
)[["snn"]]
gc_fit <- Seurat::FindClusters(
  snn, algorithm = 1, resolution = 1.5, random.seed = 20260722L,
  n.start = 20, n.iter = 20, verbose = FALSE
)
raw_gc <- as.integer(as.factor(gc_fit[, ncol(gc_fit)]))
names(raw_gc) <- rownames(pcs)
stopifnot(identical(sort(unique(raw_gc)), 1:11))

merge_map <- c(
  "1" = "S-1", "3" = "S-1", "9" = "S-1",
  "2" = "S-2", "5" = "S-2", "6" = "S-2", "11" = "S-2",
  "4" = "S-3", "8" = "S-4", "10" = "S-4", "7" = "S-5"
)
s_to_e <- c("S-1" = "E1", "S-2" = "E2", "S-3" = "E3", "S-4" = "E4", "S-5" = "E5")
hc_to_e <- c("1" = "E3", "2" = "E2", "3" = "E4", "4" = "E1", "5" = "E5")
merged_gc <- unname(merge_map[as.character(raw_gc)])
gc_e <- unname(s_to_e[merged_gc])
hc_e <- unname(hc_to_e[as.character(hc5)])
names(gc_e) <- names(hc_e) <- rownames(pcs)

gc_active <- gc_e[source[[id]]]
hc_e_active <- hc_e[source[[id]]]
gc_exact <- identical(unname(gc_active), as.character(source$Seurat_E))
hc_e_exact <- identical(unname(hc_e_active), as.character(source$HC_E))
consensus <- gc_active == hc_e_active
frozen_consensus <- tolower(as.character(frozen_class$HC_GC_consensus)) %in% c("true", "t", "1")
idx <- match(source[[id]], frozen_class[[id]])
consensus_exact <- identical(unname(consensus), frozen_consensus[idx])

# Refit the final fixed t-SNE from the archived exact PCA score matrix.  This
# avoids chaotic amplification of otherwise negligible (~1e-14) PCA rounding
# differences while still rerunning the t-SNE algorithm in the locked runtime.
working <- readRDS(working_rds)
archived_pcs <- as.matrix(working$pcs)
stopifnot(identical(rownames(archived_pcs), rownames(pcs)))
max_archived_pc_diff <- max(abs(archived_pcs - pcs))
set.seed(777L)
tsne <- Rtsne::Rtsne(
  archived_pcs, dims = 2, perplexity = 45, theta = 0.5,
  check_duplicates = FALSE, pca = FALSE, normalize = TRUE,
  max_iter = 1500L, eta = 200, exaggeration_factor = 4,
  verbose = FALSE
)$Y
rownames(tsne) <- rownames(pcs)
stored_tsne <- as.matrix(source[, c("tSNE_1", "tSNE_2")])
rerun_tsne <- tsne[source[[id]], , drop = FALSE]
max_tsne_diff <- max(abs(rerun_tsne - stored_tsne))

# Orthogonal Procrustes diagnostic.
a <- scale(rerun_tsne, center = TRUE, scale = FALSE)
b <- scale(stored_tsne, center = TRUE, scale = FALSE)
sv <- svd(t(a) %*% b)
rotation <- sv$u %*% t(sv$v)
aligned <- a %*% rotation
scale_factor <- sum(aligned * b) / sum(aligned * aligned)
procrustes_rmse <- sqrt(mean((aligned * scale_factor - b)^2))

rerun_assignments <- data.frame(
  MSN_unique_ID = rownames(pcs), PC1 = pcs[, 1], PC2 = pcs[, 2], PC3 = pcs[, 3],
  HC5 = paste0("HC", hc5), GC_raw = paste0("S", raw_gc), GC_merged = merged_gc,
  HC_E = hc_e, GC_E = gc_e, HC_GC_consensus = hc_e == gc_e,
  tSNE_1 = tsne[, 1], tSNE_2 = tsne[, 2], check.names = FALSE
)
write.csv(rerun_assignments, file.path(out, "Mouse_E_core_taxonomy_rerun_494.csv"), row.names = FALSE)

audit <- data.frame(
  check = c(
    "Active z-score matrix maximum absolute difference", "First three PCA scores maximum absolute difference",
    "Ward.D2 HC K5 assignments exact", "HC E labels exact", "Merged GC E labels exact",
    "GC-HC consensus mask exact", "Active consensus count", "Active E1 count", "Active E2 count",
    "Active E3 count", "Active E4 count", "Active E5 count",
    "Archived versus recomputed PCA maximum absolute difference",
    "Fixed t-SNE maximum absolute coordinate difference", "Fixed t-SNE Procrustes RMSE"
  ),
  observed = c(
    max_z_diff, max_pc_diff, hc_exact, hc_e_exact, gc_exact, consensus_exact,
    sum(consensus), sum(gc_active[consensus] == "E1"), sum(gc_active[consensus] == "E2"),
    sum(gc_active[consensus] == "E3"), sum(gc_active[consensus] == "E4"),
    sum(gc_active[consensus] == "E5"), max_archived_pc_diff, max_tsne_diff, procrustes_rmse
  ),
  expected = c(0, 0, TRUE, TRUE, TRUE, TRUE, 450, 142, 160, 45, 61, 42, 0, 0, 0),
  pass = c(
    max_z_diff < 1e-10, max_pc_diff < 1e-10, hc_exact, hc_e_exact, gc_exact,
    consensus_exact, sum(consensus) == 450,
    sum(gc_active[consensus] == "E1") == 142, sum(gc_active[consensus] == "E2") == 160,
    sum(gc_active[consensus] == "E3") == 45, sum(gc_active[consensus] == "E4") == 61,
    sum(gc_active[consensus] == "E5") == 42, max_archived_pc_diff < 1e-10,
    max_tsne_diff < 1e-8, procrustes_rmse < 1e-8
  ),
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(out, "core_taxonomy_concordance.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
writeLines(capture.output(installed.packages()[, c("Package", "Version")]), file.path(out, "R_installed_packages.txt"))

if (!all(audit$pass)) {
  print(audit)
  stop("At least one frozen core-taxonomy reproduction check failed.")
}
print(audit)
