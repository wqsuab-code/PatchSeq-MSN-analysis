# Frozen macaque MSN electrophysiology E1-E4 and T-E analysis

This module contains the audited macaque medium spiny neuron (MSN)
electrophysiology analysis from the public HMBA macaque Patch-seq release through
the transcriptome-to-electrophysiology (T-E) reduced-rank regression display.
The scope is restricted to caudate (Ca), putamen (Pu) and nucleus accumbens
(NAc). Macaque morphology and all mouse data are excluded from E-class
construction.

## Frozen analysis

- Stage-1 Ca/Pu/NAc cohort: 455 cells.
- Strict 19-feature complete-case cohort: 390 cells from 55 donors; no
  imputation.
- Transformation: positive ratio features were log2 transformed; the other
  features were minimum-shifted, column-sum normalized to 10,000 and log1p
  transformed; all features were Z-scored with `ddof=0`.
- Classification space: PC1-PC3 (54.6211% cumulative variance).
- HC: Euclidean Ward.D2, K=4.
- GC: Seurat SNN (`k.param=20`, `prune.SNN=1/15`) and Louvain
  (`algorithm=1`, resolution 3.0, seed 777), producing 15 raw communities.
- GC merge: maximum-overlap surjective mapping of GC15 to the four HC classes.
- Consensus cohort: 368/390 cells (94.36%) from 54 donors; E1=59, E2=133,
  E3=57 and E4=119.
- Consensus transcriptomic annotations: D1=163, D2=183 and Hybrid=22.
- T-E RRR: 346 stable D1/D2 consensus cells, 1,000 variable genes, 20
  transcriptomic PCs, 19 frozen E features and rank 3. Rank 3 is a fixed
  descriptive display choice, not a donor-held-out predictive rank.

Archived internal labels C1-C4 map one-to-one to manuscript labels E1-E4.

## Organization

- `config/`: frozen parameters, feature names and colors.
- `code/`: analysis code organized by pipeline stage.
- `data/`: source manifests and analysis-ready frozen data.
- `results/`: E-class statistics, ML500 tables and T-E RRR source data.
- `figures/`: final main figure and supporting vector/raster panels.
- `interactive/`: offline editors, embedded source data and final layouts.
- `docs/`: Methods, workflow records, figure legends and reporting materials.
- `environment/`: software versions and missing R-session disclosure.
- `tests/`: clean-clone integrity tests.
- `audit/`: conflicts, portability findings, file index and checksums.
- `legacy/`: exclusions and non-primary historical branches.

## Source-data policy

The third-party `Macaque-PatchSeq-BG.zip` archive and the 63.8-MB expression
H5AD are not redistributed here. Their official URLs, release version, byte
sizes and SHA-256 values are recorded in `data/00_source_manifest`. All
analysis-ready derived inputs required to audit the frozen labels are included.

## Important limitations

The exact historical script that generated the final19 Seurat GC scan was not
located, and the original R/Seurat `sessionInfo()` was not recorded. Frozen GC
assignments, merge maps, parameters, seed scans and downstream results are
included. Python preprocessing, PCA and Ward clustering were independently
reconstructed; Ward labels reproduce the frozen partition with ARI=1. These are
provenance limitations, not conflicting scientific results.

Some archived scripts preserve historical Windows paths. See
`audit/PORTABILITY_AUDIT.md`. The integrity test itself is repository-relative
and runs in a clean clone.

## Integrity

From the repository root run:

```bash
python analyses/macaque_ephys/tests/test_archive_integrity.py
```

## Review status

This is a collaborator-review branch, not a formal release. Do not create a
release, tag or archival DOI until licenses, authorship and repository-wide
publication policy have been finalized.
