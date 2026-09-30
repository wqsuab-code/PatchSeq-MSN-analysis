# Analysis-status register

## Final primary analysis

- 549-cell stage-1 input and 494-cell E-QC complete-case cohort.
- 493-cell, 18-feature shift-sum-log-Z taxonomy using three PCs.
- HC K = 5 and merged GC K = 5; 450-cell GC-HC consensus.
- Raw-value 18-feature statistics in the 450 consensus cells.
- Nested stratified linear-SVM validation in the 450 consensus cells.
- Strict 441-cell D1/D2 T-E RRR using merged GC K = 5 and HC K = 5.

## Sensitivity and robustness analyses

- Recording-date-grouped seven-model machine-learning validation.
- Yeo-Johnson preprocessing branch.
- GC-unmerged K = 11 versus HC K = 13 RRR.
- Alternative heat-map display limits of ±2 and ±4.

## Historical analyses

- The 551-cell stage-1 workbook.
- A legacy D1/D2 composition assigning all 450 consensus cells without the
  strict stability requirement.
- Outputs with 439, 451 or other superseded cohort sizes.

## Display-only processing

- Fixed t-SNE coordinates and class boundaries.
- Heat-map truncation.
- Radar radial mapping and clipping.
- Interactive-editor layout JSON and label positions.
- Representative traces retained as descriptive provenance only.

Reference modules in this repository were used for organization and visual
conventions only. No macaque or morphology data, labels, predictors, PCA scores
or preprocessing parameters enter the mouse E taxonomy.
