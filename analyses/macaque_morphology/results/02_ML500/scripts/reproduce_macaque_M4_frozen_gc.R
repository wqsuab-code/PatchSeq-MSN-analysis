#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(Seurat);library(mclust);library(jsonlite)})
root <- "C:/Users/53461/OneDrive/Documentos/Patch-seq_Mouse_Acb_MSN_T-type_Visualization"
out <- file.path(root, "outputs/Macaque_M4_full_ML500/tables")
pc <- fread(file.path(out, "03_recomputed_PCA_scores_126.csv")); setorder(pc, cell_label)
frozen <- fread(file.path(root, "macaque_m/m18_tempfreeze_NPC5_HCK4_res2.3/01_temp_frozen_assignments_126.csv")); setorder(frozen, cell_label)
stopifnot(identical(pc$cell_label, frozen$cell_label), nrow(pc) == 126L)
x <- as.matrix(pc[, paste0("PC", 1:5), with=FALSE]); rownames(x) <- pc$cell_label
snn <- FindNeighbors(x, k.param=20L, compute.SNN=TRUE, prune.SNN=1/15, verbose=FALSE)[["snn"]]
set.seed(777L)
fit <- FindClusters(snn, algorithm=1, resolution=2.3, random.seed=777L, n.start=30, n.iter=30, verbose=FALSE)
gc <- as.integer(as.factor(fit[, ncol(fit)]))
raw_ari <- adjustedRandIndex(gc, frozen$GC_raw_K13)
tab <- table(HC=factor(frozen$HC_K4, 1:4), GC=factor(gc, 1:13))
map <- apply(tab, 2, which.max); merged <- unname(map[gc])
merged_ari <- adjustedRandIndex(merged, frozen$GC_merged_K4)
agreement <- mean(merged == frozen$GC_merged_K4)
cell <- data.table(cell_label=pc$cell_label, recomputed_GC_raw=gc,
                   recomputed_GC_merged=merged, frozen_GC_raw=frozen$GC_raw_K13,
                   frozen_GC_merged=frozen$GC_merged_K4,
                   raw_same=gc==frozen$GC_raw_K13, merged_same=merged==frozen$GC_merged_K4)
fwrite(cell, file.path(out, "07_frozen_GC_reproduction_cells.csv"), bom=TRUE)
result <- list(raw_GC_K=uniqueN(gc), raw_GC_ARI_label_invariant=raw_ari,
               merged_GC_ARI=merged_ari, merged_GC_exact_agreement=agreement,
               merged_mismatch_n=sum(!cell$merged_same), parameters=list(NPC=5,kNN=20,
               prune_SNN=1/15,algorithm="Louvain 1",resolution=2.3,seed=777,n_start=30,n_iter=30))
write_json(result, file.path(out, "07_frozen_GC_reproduction_summary.json"), pretty=TRUE, auto_unbox=TRUE)
print(result)
