#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

project_dir <- normalizePath(
  "C:/Users/53461/OneDrive/Documentos/Patch-seq_Mouse_Acb_MSN_T-type_Visualization",
  winslash = "/",
  mustWork = TRUE
)

raw_file <- file.path(
  project_dir,
  "outputs", "e_type_qc",
  "fixed_npc3_res1.5_merged_k5_tsne",
  "consensus_k5_raw_feature_statistics",
  "ConsensusK5_Raw18_Consensus451_Values.csv"
)
label_file <- file.path(
  project_dir,
  "outputs", "e_type_qc",
  "fixed_npc3_res1.5_merged_k5_tsne",
  "p45_i1500_seed777_four_label_views",
  "Ttype_on_Econsensus",
  "D1_D2_subtype_tSNE_cell_data.csv"
)
output_dir <- file.path(
  project_dir,
  "outputs", "e_type_qc",
  "fixed_npc3_res1.5_merged_k5_tsne",
  "D1_vs_D2_within_Eclass"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

e_levels <- paste0("S-", 1:5)
major_levels <- c("D1", "D2")

feature_map <- data.table(
  Feature_order = seq_len(18L),
  Feature = c(
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
  ),
  Feature_label = c(
    "Holding MP",
    "Input resistance",
    "Membrane time constant",
    "Rheobase",
    "Sag ratio",
    "Sag time",
    "AP threshold",
    "AP amplitude",
    "AP width",
    "Upstroke/downstroke ratio",
    "Afterhyperpolarization",
    "Maximum number of APs",
    "Latency",
    "Latency at -20 pA",
    "ISI adaptation index",
    "ISI coefficient of variation",
    "AP amplitude adaptation index",
    "AP coefficient of variation"
  )
)

raw_dt <- fread(raw_file, check.names = FALSE)
label_dt <- fread(label_file, check.names = FALSE)

required_label_columns <- c(
  "MSN_unique_ID",
  "Consensus",
  "Final_Consensus_K5",
  "RPCA_major_class",
  "Final_Subtype"
)
stopifnot(
  nrow(raw_dt) == 451L,
  uniqueN(raw_dt$MSN_unique_ID) == 451L,
  all(feature_map$Feature %in% names(raw_dt)),
  all(required_label_columns %in% names(label_dt)),
  nrow(label_dt) == 494L,
  uniqueN(label_dt$MSN_unique_ID) == 494L
)

analysis_dt <- merge(
  raw_dt[
    ,
    c(
      "MSN_unique_ID",
      feature_map$Feature
    ),
    with = FALSE
  ],
  label_dt[
    Consensus == TRUE,
    ..required_label_columns
  ],
  by = "MSN_unique_ID",
  all = FALSE,
  sort = FALSE
)
setcolorder(
  analysis_dt,
  c(
    "MSN_unique_ID",
    "Final_Consensus_K5",
    "RPCA_major_class",
    "Final_Subtype",
    feature_map$Feature
  )
)

stopifnot(
  nrow(analysis_dt) == 451L,
  uniqueN(analysis_dt$MSN_unique_ID) == 451L,
  all(analysis_dt$Final_Consensus_K5 %in% e_levels),
  all(analysis_dt$RPCA_major_class %in% major_levels),
  !anyNA(analysis_dt[, ..feature_map$Feature])
)

analysis_dt[
  ,
  Final_Consensus_K5 := factor(
    Final_Consensus_K5,
    levels = e_levels
  )
]
analysis_dt[
  ,
  RPCA_major_class := factor(
    RPCA_major_class,
    levels = major_levels
  )
]

cliffs_delta <- function(x, y) {
  comparisons <- outer(x, y, FUN = "-")
  (
    sum(comparisons > 0) -
      sum(comparisons < 0)
  ) / length(comparisons)
}

summarize_vector <- function(values) {
  c(
    N = length(values),
    Mean = mean(values),
    SD = sd(values),
    Median = median(values),
    Q1 = unname(quantile(values, 0.25, type = 7)),
    Q3 = unname(quantile(values, 0.75, type = 7))
  )
}

result_records <- vector(
  "list",
  length(e_levels) * nrow(feature_map)
)
record_index <- 0L

for (e_class in e_levels) {
  e_dt <- analysis_dt[
    Final_Consensus_K5 == e_class
  ]
  for (feature_index in seq_len(nrow(feature_map))) {
    feature_name <- feature_map$Feature[feature_index]
    d1_values <- e_dt[
      RPCA_major_class == "D1",
      get(feature_name)
    ]
    d2_values <- e_dt[
      RPCA_major_class == "D2",
      get(feature_name)
    ]
    d1_summary <- summarize_vector(d1_values)
    d2_summary <- summarize_vector(d2_values)
    test_result <- wilcox.test(
      d1_values,
      d2_values,
      alternative = "two.sided",
      exact = FALSE,
      correct = FALSE
    )
    delta <- cliffs_delta(d1_values, d2_values)
    record_index <- record_index + 1L
    result_records[[record_index]] <- data.table(
      E_class = e_class,
      Feature_order = feature_map$Feature_order[feature_index],
      Feature = feature_name,
      Feature_label = feature_map$Feature_label[feature_index],
      N_D1 = as.integer(d1_summary["N"]),
      N_D2 = as.integer(d2_summary["N"]),
      Mean_D1 = d1_summary["Mean"],
      SD_D1 = d1_summary["SD"],
      Median_D1 = d1_summary["Median"],
      Q1_D1 = d1_summary["Q1"],
      Q3_D1 = d1_summary["Q3"],
      Mean_D2 = d2_summary["Mean"],
      SD_D2 = d2_summary["SD"],
      Median_D2 = d2_summary["Median"],
      Q1_D2 = d2_summary["Q1"],
      Q3_D2 = d2_summary["Q3"],
      Median_difference_D1_minus_D2 =
        d1_summary["Median"] - d2_summary["Median"],
      Cliff_delta_D1_vs_D2 = delta,
      Effect_direction = fifelse(
        delta > 0,
        "D1 higher",
        fifelse(delta < 0, "D2 higher", "No direction")
      ),
      Wilcoxon_W = unname(test_result$statistic),
      P_value = test_result$p.value
    )
  }
}

result_dt <- rbindlist(result_records)
result_dt[
  ,
  BH_q_global := p.adjust(P_value, method = "BH")
]
result_dt[
  ,
  BH_q_within_E_class := p.adjust(P_value, method = "BH"),
  by = E_class
]
result_dt[
  ,
  Significance_global := fifelse(
    BH_q_global < 0.001,
    "***",
    fifelse(
      BH_q_global < 0.01,
      "**",
      fifelse(BH_q_global < 0.05, "*", "ns")
    )
  )
]
result_dt[
  ,
  Significant_global_FDR_0.05 := BH_q_global < 0.05
]
result_dt[
  ,
  Significant_within_E_FDR_0.05 := BH_q_within_E_class < 0.05
]
setorder(result_dt, Feature_order, E_class)

significant_dt <- copy(
  result_dt[Significant_global_FDR_0.05 == TRUE]
)
significant_dt[
  ,
  Abs_Cliff_delta := abs(Cliff_delta_D1_vs_D2)
]
setorder(
  significant_dt,
  BH_q_global,
  -Abs_Cliff_delta
)

sample_size_dt <- dcast(
  analysis_dt[
    ,
    .N,
    by = .(
      E_class = Final_Consensus_K5,
      T_major = RPCA_major_class
    )
  ],
  E_class ~ T_major,
  value.var = "N",
  fill = 0
)
sample_size_dt[
  ,
  Total := D1 + D2
]

feature_summary_dt <- result_dt[
  ,
  {
    strongest_index <- which.max(
      abs(Cliff_delta_D1_vs_D2)
    )
    list(
      Significant_E_classes_global_FDR =
        sum(Significant_global_FDR_0.05),
      Significant_E_classes_within_E_FDR =
        sum(Significant_within_E_FDR_0.05),
      Strongest_E_class = E_class[strongest_index],
      Strongest_Cliff_delta =
        Cliff_delta_D1_vs_D2[strongest_index],
      Strongest_direction =
        Effect_direction[strongest_index],
      Minimum_global_q = min(BH_q_global)
    )
  },
  by = .(
    Feature_order,
    Feature,
    Feature_label
  )
]
setorder(
  feature_summary_dt,
  -Significant_E_classes_global_FDR,
  Minimum_global_q,
  Feature_order
)

effect_matrix_dt <- dcast(
  result_dt,
  Feature_order + Feature + Feature_label ~ E_class,
  value.var = "Cliff_delta_D1_vs_D2"
)
q_matrix_dt <- dcast(
  result_dt,
  Feature_order + Feature + Feature_label ~ E_class,
  value.var = "BH_q_global"
)
setorder(effect_matrix_dt, Feature_order)
setorder(q_matrix_dt, Feature_order)

fwrite(
  result_dt,
  file.path(
    output_dir,
    "D1_vs_D2_Within_Eclass_All18_Wilcoxon.csv"
  )
)
fwrite(
  significant_dt,
  file.path(
    output_dir,
    "D1_vs_D2_Within_Eclass_GlobalFDR_Significant.csv"
  )
)
fwrite(
  feature_summary_dt,
  file.path(
    output_dir,
    "D1_vs_D2_Within_Eclass_Feature_Summary.csv"
  )
)
fwrite(
  sample_size_dt,
  file.path(
    output_dir,
    "D1_vs_D2_Within_Eclass_Sample_Sizes.csv"
  )
)
fwrite(
  effect_matrix_dt,
  file.path(
    output_dir,
    "D1_vs_D2_Within_Eclass_CliffDelta_Matrix.csv"
  )
)
fwrite(
  q_matrix_dt,
  file.path(
    output_dir,
    "D1_vs_D2_Within_Eclass_GlobalQ_Matrix.csv"
  )
)
fwrite(
  analysis_dt,
  file.path(
    output_dir,
    "D1_vs_D2_Within_Eclass_Analysis_Data_451.csv"
  )
)

sample_labels <- setNames(
  sprintf(
    "%s\n%d/%d",
    sample_size_dt$E_class,
    sample_size_dt$D1,
    sample_size_dt$D2
  ),
  sample_size_dt$E_class
)

heatmap_dt <- copy(result_dt)
heatmap_dt[
  ,
  E_class := factor(E_class, levels = rev(e_levels))
]
heatmap_dt[
  ,
  Feature_label := factor(
    Feature_label,
    levels = feature_map$Feature_label
  )
]
heatmap_dt[
  ,
  Heatmap_label := sprintf("%.2f", Cliff_delta_D1_vs_D2)
]
heatmap_dt[
  ,
  Label_colour := ifelse(
    abs(Cliff_delta_D1_vs_D2) >= 0.55,
    "white",
    "black"
  )
]

heatmap_plot <- ggplot(
  heatmap_dt,
  aes(
    x = Feature_label,
    y = E_class,
    fill = Cliff_delta_D1_vs_D2
  )
) +
  geom_tile(
    colour = "#D9D9D9",
    linewidth = 0.35 / ggplot2::.pt
  ) +
  geom_tile(
    data = heatmap_dt[Significant_global_FDR_0.05 == TRUE],
    fill = NA,
    colour = "black",
    linewidth = 0.8 / ggplot2::.pt
  ) +
  geom_text(
    aes(
      label = Heatmap_label,
      colour = Label_colour
    ),
    angle = 0,
    size = 4.7 / ggplot2::.pt,
    family = "Arial"
  ) +
  scale_colour_identity() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Cliff's delta\nD1 vs D2"
  ) +
  scale_x_discrete(
    position = "bottom",
    drop = FALSE
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "D1 vs D2 within consensus E classes",
    subtitle = "Raw E-feature values; black borders indicate global BH q<0.05",
    x = NULL,
    y = NULL,
    caption = paste0(
      "Positive Cliff's delta: D1 higher; negative Cliff's delta: D2 higher."
    )
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal(
    base_size = 7,
    base_family = "Arial"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 5.2,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 6,
      face = "bold"
    ),
    plot.title = element_text(
      colour = "black",
      size = 8,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      colour = "black",
      size = 6,
      hjust = 0
    ),
    plot.caption = element_text(
      colour = "black",
      size = 5.5,
      hjust = 0
    ),
    legend.title = element_text(
      colour = "black",
      size = 6
    ),
    legend.text = element_text(
      colour = "black",
      size = 5.5
    ),
    legend.key.height = grid::unit(14, "mm"),
    plot.margin = margin(4, 4, 4, 4, unit = "pt")
  )

png_file <- file.path(
  output_dir,
  "D1_vs_D2_Within_Eclass_CliffDelta_Heatmap.png"
)
pdf_file <- file.path(
  output_dir,
  "D1_vs_D2_Within_Eclass_CliffDelta_Heatmap.pdf"
)
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(
    png_file,
    heatmap_plot,
    width = 6.0,
    height = 4.2,
    units = "in",
    dpi = 600,
    device = ragg::agg_png,
    bg = "white"
  )
} else {
  ggsave(
    png_file,
    heatmap_plot,
    width = 6.0,
    height = 4.2,
    units = "in",
    dpi = 600,
    bg = "white"
  )
}
ggsave(
  pdf_file,
  heatmap_plot,
  width = 6.0,
  height = 4.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

significant_counts <- result_dt[
  ,
  .(
    Significant_global_FDR_0.05 =
      sum(Significant_global_FDR_0.05),
    Significant_within_E_FDR_0.05 =
      sum(Significant_within_E_FDR_0.05)
  ),
  by = E_class
]

writeLines(
  c(
    "D1 versus D2 within each consensus E class",
    sprintf("Analysis cells: %d", nrow(analysis_dt)),
    sprintf("Features: %d", nrow(feature_map)),
    sprintf(
      "Comparisons: %d",
      nrow(result_dt)
    ),
    paste0(
      "Sample sizes: ",
      paste(
        sprintf(
          "%s D1=%d D2=%d",
          sample_size_dt$E_class,
          sample_size_dt$D1,
          sample_size_dt$D2
        ),
        collapse = "; "
      )
    ),
    "Input values: original raw electrophysiological measurements",
    "Primary test: two-sided Wilcoxon rank-sum, exact=FALSE, no continuity correction",
    "Effect size: Cliff's delta; positive values indicate D1 higher than D2",
    "Primary multiplicity control: Benjamini-Hochberg across all 90 tests",
    "Secondary multiplicity control: Benjamini-Hochberg within each E class (18 tests)",
    sprintf(
      "Global-FDR significant comparisons: %d/%d",
      nrow(significant_dt),
      nrow(result_dt)
    ),
    paste0(
      "Global-FDR significant by E class: ",
      paste(
        sprintf(
          "%s=%d",
          significant_counts$E_class,
          significant_counts$Significant_global_FDR_0.05
        ),
        collapse = "; "
      )
    ),
    "Caution: S-5 contains only 6 D2 cells; effect estimates are reported but power is limited.",
    sprintf(
      "Output directory: %s",
      normalizePath(output_dir, winslash = "/")
    )
  ),
  file.path(
    output_dir,
    "D1_vs_D2_Within_Eclass_run_log.txt"
  )
)

cat(normalizePath(output_dir, winslash = "/"), "\n")
cat(sprintf(
  "Global-FDR significant comparisons: %d/%d\n",
  nrow(significant_dt),
  nrow(result_dt)
))
print(significant_counts)
