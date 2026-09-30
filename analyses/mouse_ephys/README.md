# Frozen mouse nucleus accumbens electrophysiology E-type analysis

This module archives the audited mouse nucleus accumbens medium spiny neuron
electrophysiology taxonomy and its transcriptome-to-electrophysiology (T-E)
integration. The primary taxonomy was defined independently by graph clustering
(GC) and hierarchical clustering (HC); expert E-type labels were not used to fit
or define the frozen classes.

## Frozen primary analysis

- Raw electrophysiology archive: 585 records.
- Main eligible cohort before stage-1 exclusions: 556 cells.
- Stage-1 cohort: 549 cells.
- E-QC complete-case cohort: 494 cells with 25 candidate features.
- Active taxonomy cohort: 493 cells with 18 frozen features.
- Classification space: the first three principal components.
- HC: Euclidean distance, Ward.D2 linkage, K = 5.
- GC: SNN/Louvain, k = 20, pruning = 1/15, resolution = 1.5; 11 raw
  communities merged to five prespecified groups.
- GC-HC consensus: 450 of 493 cells (91.28%); E1 = 142, E2 = 160,
  E3 = 45, E4 = 61 and E5 = 42.
- Strict transcriptomic identity: 441 cells; D1 = 209 and D2 = 232. Nine
  additional consensus cells were transcriptomically ambiguous.
- T-E RRR: 1,000 variable genes, 20 transcriptomic PCs, 18 E-features,
  rank = 3 and seed = 777.

## Primary numerical transformation

For each E-feature, the cohort minimum was subtracted, shifted values were
divided by their column sum and multiplied by 10,000, and values were transformed
with `log1p` and feature-wise Z-scoring. Yeo-Johnson preprocessing is retained
only as a sensitivity branch and was not used for the frozen taxonomy.

## Repository contents

- `config/`: frozen parameters, feature order and exclusion list.
- `code/`: analysis scripts organized by stage.
- `data/`: source manifest, cohort/QC tables, frozen inputs, transformed matrices
  and cell-level labels.
- `results/`: statistics, machine-learning validation and strict-441 RRR tables.
- `figures/`: the current main SVG and selected vector/raster extended figures.
- `interactive/`: offline figure editors and frozen display payloads.
- `docs/`: Methods, figure legend, data dictionary and release notes.
- `environment/`: recorded R and Python package snapshots.
- `audit/`: conflict resolution, portability scans, file indices and checksums.
- `legacy/`: historical, sensitivity or display-only materials excluded from the
  primary inferential chain.

## Validation designs

The nested linear-SVM analysis used five outer and four inner stratified folds,
with all transformations fitted within training folds. Its out-of-fold accuracy
was 0.9200. A separate robustness analysis used 200 recording-date-grouped
75:25 train-test splits and seven classifier families. Sequencing-submission
batches were retained as metadata, not used as the grouping variable. These two
performance estimates have different estimands and must not be conflated.

## Data policy

Analysis-ready derived data and cell-level labels required to audit conclusions
are versioned here. Original Excel workbooks, ABF recordings and the 143 MB RNA
count object are not committed pending final redistribution and accession review.
Their filenames, sizes and SHA-256 hashes are recorded in
`data/00_source_manifest/`. The RNA object and frozen review ZIP are staged in a
local `release-assets` directory outside Git.

## Reproducibility status

The 24 archived concordance checks passed. `tests/test_archive_integrity.py`
verifies the source manifest, frozen feature order, all principal cohort sizes,
class counts, cell-ID relationships and critical missingness using only files in
a clean clone. Full raw-data reruns still require the externally staged source
files and path configuration described in `audit/PORTABILITY_AUDIT.md`.

This is a collaborator-review module, not a public release. License, author list,
permanent accessions and archival DOI remain pending.
