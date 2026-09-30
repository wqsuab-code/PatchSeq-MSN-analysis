#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)
set.seed(777)

zscore_file <- file.path(
  "outputs", "e_type_qc", "pca_shift_sum_final18",
  "ShiftSumNorm10000_Log1p_Zscore_Final18_PCA_input.csv"
)
raw_file <- file.path(
  "outputs", "e_type_qc", "parallel_preprocessing",
  "Ephys_QCpass_494_raw_25features.csv"
)
output_dir <- file.path(
  "outputs", "e_type_qc", "hc5_tree_distillation"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

zscore_dt <- fread(zscore_file)
raw_dt <- fread(raw_file)
feature_names <- setdiff(names(zscore_dt), "MSN_unique_ID")
stopifnot(
  nrow(zscore_dt) == 494L,
  length(feature_names) == 18L,
  uniqueN(zscore_dt$MSN_unique_ID) == 494L,
  uniqueN(raw_dt$MSN_unique_ID) == 494L
)

z_matrix <- as.matrix(zscore_dt[, ..feature_names])
storage.mode(z_matrix) <- "double"
rownames(z_matrix) <- zscore_dt$MSN_unique_ID
stopifnot(!anyNA(z_matrix), all(is.finite(z_matrix)))

pca_fit <- prcomp(z_matrix, center = FALSE, scale. = FALSE)
pcs <- pca_fit$x[, 1:3, drop = FALSE]
hc <- hclust(dist(pcs, method = "euclidean"), method = "ward.D2")
hc5 <- cutree(hc, k = 5L)
n_cells <- nrow(pcs)

node_leaves <- function(node) {
  if (node < 0L) {
    return(-node)
  }
  c(node_leaves(hc$merge[node, 1L]), node_leaves(hc$merge[node, 2L]))
}

node_classes <- function(node) {
  sort(unique(unname(hc5[node_leaves(node)])))
}

node_size <- function(node) {
  length(node_leaves(node))
}

internal_rows <- list()
collect_compressed_nodes <- function(node, parent = NA_integer_, depth = 0L) {
  classes <- node_classes(node)
  if (length(classes) <= 1L) {
    return(invisible(NULL))
  }
  left <- hc$merge[node, 1L]
  right <- hc$merge[node, 2L]
  left_classes <- node_classes(left)
  right_classes <- node_classes(right)
  row_id <- length(internal_rows) + 1L
  internal_rows[[row_id]] <<- data.table(
    Tree_node = node,
    Parent_tree_node = parent,
    Depth = depth,
    Classes = paste(classes, collapse = "|"),
    N = node_size(node),
    Left_node = left,
    Left_classes = paste(left_classes, collapse = "|"),
    Left_N = node_size(left),
    Right_node = right,
    Right_classes = paste(right_classes, collapse = "|"),
    Right_N = node_size(right),
    Ward_height = hc$height[node]
  )
  if (length(left_classes) > 1L) {
    collect_compressed_nodes(left, parent = node, depth = depth + 1L)
  }
  if (length(right_classes) > 1L) {
    collect_compressed_nodes(right, parent = node, depth = depth + 1L)
  }
  invisible(NULL)
}

root_node <- n_cells - 1L
collect_compressed_nodes(root_node)
topology <- rbindlist(internal_rows)
setorder(topology, Depth, -N)
stopifnot(nrow(topology) == 4L)

pc_export <- data.table(
  MSN_unique_ID = rownames(pcs),
  PC1 = pcs[, 1L],
  PC2 = pcs[, 2L],
  PC3 = pcs[, 3L],
  HC5 = paste0("HC", unname(hc5))
)
fwrite(
  pc_export,
  file.path(output_dir, "NPC3_Frozen_PC_Scores_HC5.csv")
)
fwrite(
  topology,
  file.path(output_dir, "NPC3_WardD2_HC5_Compressed_Topology.csv")
)

raw_aligned <- raw_dt[
  match(pc_export$MSN_unique_ID, MSN_unique_ID),
  c("MSN_unique_ID", feature_names),
  with = FALSE
]
stopifnot(identical(raw_aligned$MSN_unique_ID, pc_export$MSN_unique_ID))
fwrite(
  raw_aligned,
  file.path(output_dir, "NPC3_HC5_Raw_Final18_Features.csv")
)

loading_export <- as.data.table(
  pca_fit$rotation[, 1:3, drop = FALSE],
  keep.rownames = "Feature"
)
fwrite(
  loading_export,
  file.path(output_dir, "NPC3_PCA_Loadings_Final18.csv")
)

topology_text <- vapply(
  seq_len(nrow(topology)),
  function(i) {
    sprintf(
      "Node %d [n=%d; HC%s] -> left HC%s [n=%d] | right HC%s [n=%d]",
      topology$Tree_node[i],
      topology$N[i],
      gsub("\\|", ",HC", topology$Classes[i]),
      gsub("\\|", ",HC", topology$Left_classes[i]),
      topology$Left_N[i],
      gsub("\\|", ",HC", topology$Right_classes[i]),
      topology$Right_N[i]
    )
  },
  character(1)
)

log_lines <- c(
  "NPC3 Ward.D2 HC5 tree-distillation input export",
  sprintf("Run time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Cells: %d", nrow(pc_export)),
  sprintf("Final electrophysiological metrics: %d", length(feature_names)),
  "Preprocessing: per-feature shift -> sum normalization x10,000 -> log1p -> z score",
  "PCA: prcomp(center=FALSE, scale.=FALSE)",
  "HC: Euclidean distance on PC1-PC3; Ward.D2; cutree(k=5)",
  "",
  "Compressed five-leaf topology:",
  topology_text
)
writeLines(
  log_lines,
  file.path(output_dir, "NPC3_WardD2_HC5_Topology_run_log.txt")
)
cat(paste(log_lines, collapse = "\n"), "\n")
