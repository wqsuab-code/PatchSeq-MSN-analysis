# Shared permanent cell exclusions for all E-type analyses.

load_e_type_exclusions <- function(root = ".") {
  path <- file.path(root, "config", "e_type_global_exclusions.txt")
  if (!file.exists(path)) {
    stop("Global E-type exclusion file is missing: ", path)
  }
  ids <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  ids <- ids[nzchar(ids) & !grepl("^#", ids)]
  if (!length(ids)) {
    stop("Global E-type exclusion file is empty: ", path)
  }
  unique(ids)
}

filter_e_type_exclusions <- function(data, id_column = "MSN_unique_ID", root = ".") {
  if (!id_column %in% names(data)) {
    stop("ID column is absent from the input table: ", id_column)
  }
  excluded <- load_e_type_exclusions(root)
  removed <- sort(intersect(as.character(data[[id_column]]), excluded))
  result <- data[!as.character(data[[id_column]]) %in% excluded, , drop = FALSE]
  attr(result, "globally_excluded_cells_removed") <- removed
  result
}
