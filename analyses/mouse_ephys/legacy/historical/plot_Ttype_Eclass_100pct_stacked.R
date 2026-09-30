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
input_file <- file.path(
  project_dir,
  "outputs", "e_type_qc",
  "fixed_npc3_res1.5_merged_k5_tsne",
  "p45_i1500_seed777_four_label_views",
  "Ttype_on_Econsensus",
  "D1_D2_subtype_tSNE_cell_data.csv"
)
subtype_color_file <- file.path(
  "C:/Users/53461/Documents/Codex/2026-07-13/zh/outputs",
  "05_branch_specific_subtype_mapping",
  "60_MSN_PC16_tSNE_complete_package",
  "labels",
  "MSN_fixed_subtype_color_mapping.csv"
)
output_dir <- file.path(
  dirname(input_file),
  "eclass_100pct_stacked"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

e_levels <- paste0("S-", 1:5)
major_levels <- c("D1", "D2")
subtype_levels <- c(
  paste0("D1_", 1:8),
  paste0("D2_", 1:8)
)
major_colors <- c(
  "D1" = "#E31A1C",
  "D2" = "#377EB8"
)

cell_dt <- fread(input_file)
stopifnot(
  nrow(cell_dt) == 494L,
  uniqueN(cell_dt$MSN_unique_ID) == 494L,
  sum(cell_dt$Consensus) == 451L,
  all(c(
    "Final_Consensus_K5",
    "RPCA_major_class",
    "Final_Subtype"
  ) %in% names(cell_dt))
)
cell_dt <- cell_dt[Consensus == TRUE]
stopifnot(
  nrow(cell_dt) == 451L,
  all(cell_dt$Final_Consensus_K5 %in% e_levels),
  all(cell_dt$RPCA_major_class %in% major_levels),
  all(cell_dt$Final_Subtype %in% subtype_levels)
)

subtype_color_dt <- fread(subtype_color_file)
stopifnot(
  all(c("Subtype", "color") %in% names(subtype_color_dt)),
  all(subtype_levels %in% subtype_color_dt$Subtype)
)
subtype_colors <- setNames(
  subtype_color_dt$color,
  subtype_color_dt$Subtype
)[subtype_levels]

make_composition_table <- function(
    input_dt,
    label_column,
    label_levels,
    label_level_name) {
  complete_grid <- CJ(
    E_class = e_levels,
    T_label = label_levels,
    unique = TRUE
  )
  observed <- input_dt[
    ,
    .N,
    by = .(
      E_class = Final_Consensus_K5,
      T_label = get(label_column)
    )
  ]
  output <- merge(
    complete_grid,
    observed,
    by = c("E_class", "T_label"),
    all.x = TRUE,
    sort = FALSE
  )
  output[is.na(N), N := 0L]
  output[
    ,
    E_total := sum(N),
    by = E_class
  ]
  output[
    ,
    Percent_within_E_class := 100 * N / E_total
  ]
  output[
    ,
    `:=`(
      T_level = label_level_name,
      E_class = factor(E_class, levels = e_levels),
      T_label = factor(T_label, levels = label_levels)
    )
  ]
  setorder(output, E_class, T_label)
  output
}

major_composition <- make_composition_table(
  cell_dt,
  "RPCA_major_class",
  major_levels,
  "Major"
)
subtype_composition <- make_composition_table(
  cell_dt,
  "Final_Subtype",
  subtype_levels,
  "Subtype"
)

stopifnot(
  all(
    abs(
      major_composition[
        ,
        sum(Percent_within_E_class),
        by = E_class
      ]$V1 - 100
    ) < 1e-10
  ),
  all(
    abs(
      subtype_composition[
        ,
        sum(Percent_within_E_class),
        by = E_class
      ]$V1 - 100
    ) < 1e-10
  )
)

e_totals <- cell_dt[
  ,
  .N,
  by = Final_Consensus_K5
][
  match(e_levels, Final_Consensus_K5)
]
e_axis_labels <- setNames(
  sprintf(
    "%s\nn=%d",
    e_totals$Final_Consensus_K5,
    e_totals$N
  ),
  e_totals$Final_Consensus_K5
)

major_composition[
  ,
  Label := sprintf(
    "%d\n%.1f%%",
    N,
    Percent_within_E_class
  )
]

major_plot <- ggplot(
  major_composition,
  aes(
    x = E_class,
    y = Percent_within_E_class,
    fill = T_label
  )
) +
  geom_col(
    width = 0.72,
    position = position_stack(reverse = TRUE),
    colour = "white",
    linewidth = 0.35 / ggplot2::.pt
  ) +
  geom_text(
    aes(label = Label),
    position = position_stack(
      vjust = 0.5,
      reverse = TRUE
    ),
    colour = "white",
    family = "Arial",
    size = 5.4 / ggplot2::.pt,
    lineheight = 0.88
  ) +
  scale_fill_manual(
    values = major_colors,
    limits = major_levels,
    drop = FALSE,
    name = NULL
  ) +
  scale_x_discrete(
    labels = e_axis_labels,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 50, 100),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "D1/D2 composition within E classes",
    x = NULL,
    y = "Cells (%)"
  ) +
  theme_classic(
    base_size = 6,
    base_family = "Arial"
  ) +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5 / ggplot2::.pt
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5 / ggplot2::.pt
    ),
    axis.ticks.length = grid::unit(0.8, "mm"),
    axis.text = element_text(
      colour = "black",
      size = 5.2
    ),
    axis.title.y = element_text(
      colour = "black",
      size = 5.6
    ),
    plot.title = element_text(
      colour = "black",
      face = "bold",
      size = 6,
      hjust = 0.5,
      margin = margin(b = 2, unit = "pt")
    ),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(
      colour = "black",
      size = 5
    ),
    legend.key.width = grid::unit(3.2, "mm"),
    legend.key.height = grid::unit(2.2, "mm"),
    legend.spacing.x = grid::unit(0.7, "mm"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    plot.margin = margin(2, 3, 2, 3, unit = "pt")
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

hex_luminance <- function(hex_color) {
  rgb_values <- col2rgb(hex_color) / 255
  0.2126 * rgb_values[1, ] +
    0.7152 * rgb_values[2, ] +
    0.0722 * rgb_values[3, ]
}
subtype_text_colors <- ifelse(
  hex_luminance(subtype_colors) < 0.52,
  "white",
  "black"
)
names(subtype_text_colors) <- names(subtype_colors)

subtype_composition[
  ,
  Label := fifelse(
    Percent_within_E_class >= 5,
    sprintf(
      "%s\n%d (%.1f%%)",
      as.character(T_label),
      N,
      Percent_within_E_class
    ),
    ""
  )
]
subtype_composition[
  ,
  Label_y := cumsum(Percent_within_E_class) -
    Percent_within_E_class / 2,
  by = E_class
]
subtype_composition[
  ,
  Label_color := subtype_text_colors[
    as.character(T_label)
  ]
]

subtype_plot <- ggplot(
  subtype_composition,
  aes(
    x = E_class,
    y = Percent_within_E_class,
    fill = T_label
  )
) +
  geom_col(
    width = 0.72,
    position = position_stack(reverse = TRUE),
    colour = "white",
    linewidth = 0.25 / ggplot2::.pt
  ) +
  geom_text(
    aes(
      y = Label_y,
      label = Label,
      colour = Label_color,
      group = T_label
    ),
    family = "Arial",
    size = 4.4 / ggplot2::.pt,
    lineheight = 0.82
  ) +
  scale_colour_identity() +
  scale_fill_manual(
    values = subtype_colors,
    limits = subtype_levels,
    drop = FALSE,
    name = NULL
  ) +
  scale_x_discrete(
    labels = e_axis_labels,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "D1/D2 subtype composition within E classes",
    subtitle = "Labels shown for segments >=5%",
    x = NULL,
    y = "Cells (%)"
  ) +
  theme_classic(
    base_size = 6,
    base_family = "Arial"
  ) +
  theme(
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5 / ggplot2::.pt
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5 / ggplot2::.pt
    ),
    axis.ticks.length = grid::unit(0.8, "mm"),
    axis.text = element_text(
      colour = "black",
      size = 5.2
    ),
    axis.title.y = element_text(
      colour = "black",
      size = 5.6
    ),
    plot.title = element_text(
      colour = "black",
      face = "bold",
      size = 7,
      hjust = 0.5,
      margin = margin(b = 1, unit = "pt")
    ),
    plot.subtitle = element_text(
      colour = "black",
      size = 5.3,
      hjust = 0.5,
      margin = margin(b = 2, unit = "pt")
    ),
    legend.position = "right",
    legend.direction = "vertical",
    legend.text = element_text(
      colour = "black",
      size = 5
    ),
    legend.key.width = grid::unit(2.6, "mm"),
    legend.key.height = grid::unit(2.6, "mm"),
    legend.spacing.y = grid::unit(0.25, "mm"),
    legend.margin = margin(0, 0, 0, 1, unit = "pt"),
    plot.margin = margin(3, 3, 2, 3, unit = "pt")
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = FALSE
    )
  )

save_plot <- function(
    filename,
    plot_object,
    width,
    height) {
  if (grepl("\\.png$", filename, ignore.case = TRUE)) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggsave(
        filename,
        plot_object,
        width = width,
        height = height,
        units = "in",
        dpi = 600,
        device = ragg::agg_png,
        bg = "white"
      )
    } else {
      ggsave(
        filename,
        plot_object,
        width = width,
        height = height,
        units = "in",
        dpi = 600,
        bg = "white"
      )
    }
  } else {
    ggsave(
      filename,
      plot_object,
      width = width,
      height = height,
      units = "in",
      device = cairo_pdf,
      bg = "white"
    )
  }
}

major_stem <- file.path(
  output_dir,
  "Pct100_D1_D2_within_E_W2p4_H1p4"
)
subtype_stem <- file.path(
  output_dir,
  "Pct100_16subtypes_within_E_W4p8_H3p2"
)
save_plot(
  paste0(major_stem, ".png"),
  major_plot,
  2.4,
  1.4
)
save_plot(
  paste0(major_stem, ".pdf"),
  major_plot,
  2.4,
  1.4
)
save_plot(
  paste0(subtype_stem, ".png"),
  subtype_plot,
  4.8,
  3.2
)
save_plot(
  paste0(subtype_stem, ".pdf"),
  subtype_plot,
  4.8,
  3.2
)

remove_text_layers <- function(plot_object) {
  is_text_layer <- vapply(
    plot_object$layers,
    function(layer) {
      inherits(layer$geom, "GeomText")
    },
    logical(1)
  )
  plot_object$layers <- plot_object$layers[!is_text_layer]
  plot_object
}

major_nolabel_plot <- remove_text_layers(major_plot)
subtype_nolabel_plot <- remove_text_layers(subtype_plot) +
  labs(subtitle = NULL)

major_nolabel_png <- file.path(
  output_dir,
  "Pct100_D1_D2_within_E_NoSegmentLabels_W2p4_H1p4.png"
)
subtype_nolabel_png <- file.path(
  output_dir,
  "Pct100_Subtype_NoLabels_W4p8_H3p2.png"
)
save_plot(
  major_nolabel_png,
  major_nolabel_plot,
  2.4,
  1.4
)
save_plot(
  subtype_nolabel_png,
  subtype_nolabel_plot,
  4.8,
  3.2
)

composition_output <- rbindlist(
  list(
    major_composition[
      ,
      .(
        T_level,
        E_class = as.character(E_class),
        T_label = as.character(T_label),
        N,
        E_total,
        Percent_within_E_class
      )
    ],
    subtype_composition[
      ,
      .(
        T_level,
        E_class = as.character(E_class),
        T_label = as.character(T_label),
        N,
        E_total,
        Percent_within_E_class
      )
    ]
  ),
  use.names = TRUE
)
fwrite(
  composition_output,
  file.path(
    output_dir,
    "Eclass_D1_D2_Subtype_100pct_Composition.csv"
  )
)

percentage_check <- composition_output[
  ,
  .(
    Percent_sum = sum(Percent_within_E_class)
  ),
  by = .(
    T_level,
    E_class
  )
]
fwrite(
  percentage_check,
  file.path(
    output_dir,
    "Eclass_100pct_Composition_Check.csv"
  )
)

writeLines(
  c(
    "D1/D2 and subtype 100% stacked composition within each E class",
    "Cohort: 451 E-consensus QC-pass cells",
    "Denominator: all consensus cells within each E class",
    "Each S-1 to S-5 bar sums to 100%",
    paste0(
      "E-class totals: ",
      paste(
        sprintf(
          "%s=%d",
          e_totals$Final_Consensus_K5,
          e_totals$N
        ),
        collapse = "; "
      )
    ),
    "Major colours: D1=#E31A1C; D2=#377EB8",
    "Subtype colours: established MSN fixed subtype palette",
    "Major labels: count and within-E percentage for every segment",
    "Subtype labels: subtype, count and within-E percentage for segments >=5%",
    "No-segment-label PNG versions retain axes, titles and legends",
    sprintf(
      "Output directory: %s",
      normalizePath(output_dir, winslash = "/")
    )
  ),
  file.path(
    output_dir,
    "Eclass_100pct_Stacked_run_log.txt"
  )
)

cat(normalizePath(output_dir, winslash = "/"), "\n")
print(percentage_check)
