# Data dictionary

## Cell identifiers

`MSN_unique_ID` is the frozen cell key used across electrophysiology matrices,
PCA scores, GC-HC assignments, statistics, machine-learning outputs and T-E RRR.
The integrity test requires uniqueness within every cell-level table and exact
set relationships between analysis stages.

## Feature matrices

- `Ephys_QCpass_494_raw_25features.csv`: raw complete-case measurements for the
  25 candidate features after E-QC.
- `Active_Final18_PCA_input.csv`: 493-cell transformed and standardized matrix
  used for the primary taxonomy. The first column is the cell ID and the
  remaining 18 columns follow `config/frozen_feature_order.txt` exactly.
- `NPC3_PCA_Loadings_Final18.csv`: feature loadings for the frozen PCA.
- `NPC3_Frozen_PC_Scores_HC5.csv`: frozen PC scores and HC assignments.

## Labels

- `HC_GC_only_cell_assignments.csv`: 493 active cells with GC, HC, consensus and
  fixed t-SNE fields. Legacy filenames retain `HC_GC`; prose uses GC-HC.
- `e_stability_all493.csv`: transcriptomic identity stability for every active
  E-analysis cell.
- `e_strict_t_identity.csv`: the 450 GC-HC consensus cells with strict identity
  flag and D1/D2 label where stable.
- `RRR_cell_scores_n441.csv`: T-E RRR scores for the 441 strict cells.

## Missing data

Primary matrices are complete case and contain no missing or nonfinite feature
values. Controlled raw sources use their original missing-data conventions.
No statistical outlier removal was applied to the 450-cell raw-value comparisons.

## Observation and grouping units

The observation unit is a cell. The recording-date-grouped robustness analysis
kept all cells recorded on the same date in the same train/test partition.
Sequencing-submission batch was retained as metadata. A unique animal identifier
is not present in the archived E-analysis tables, so animal count was not
reconstructed or substituted by recording date.
