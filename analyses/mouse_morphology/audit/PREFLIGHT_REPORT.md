# Submission preflight report

Date: 2026-09-30

## Frozen scientific scope

- Starting morphology/QC table: 228 cells.
- Post-QC morphology cohort: 193 cells.
- Frozen taxonomy: 187 cells and 10 ordered features.
- Frozen classes: M1=67, M2=23, M3=62, M4=35.
- HC-GC consensus: 181/187 (96.79%).
- Strict T-M RRR: 168 cells; D1=74, D2=94; M1=63, M2=22,
  M3=55, M4=28.

## Validation performed

- Source-manifest sizes and SHA-256 values: pass.
- 228-to-193-to-187 cohort transitions: pass.
- Ten-feature count and order: pass.
- Frozen input, PCA, labels and heat-map source IDs: pass.
- Frozen class counts and HC-GC consensus count: pass.
- 187-cell ML OOF IDs and 168-cell strict RRR IDs: pass.
- Key result tables contain no unexpected missing tokens: pass.
- ML predictor boundary (`M_` features only; no E predictors): pass.
- Common secret/token/private-key/password/e-mail scan: no matches.
- Files larger than 25 MB in the Git module: none.
- Files larger than 50 MB in the Git module: none.
- Changes outside `analyses/mouse_morphology/`: only root `README.md` and
  `analyses/README.md` module-index updates.
- Changes under frozen `analyses/macaque_morphology/`: none.

## Reproducibility level

PCA/HC and strict RRR were previously rerun from frozen inputs and reproduced to
floating-point precision. This archive run performs a fast consistency audit of
saved 500-repeat ML outputs, parameters, split logs, row counts and hashes; it
does not recompute all 500 repeats.

## Remaining portability findings

The absolute-path scan found 790 path strings across 105 files. Most are export
metadata embedded in matched PNG/PDF/SVG triplets; the remainder are archived
provenance scripts, source-table path columns and interactive builder defaults.
They are enumerated in `absolute_path_scan.txt`. Repository-relative frozen data
and the integrity test do not depend on them. Exact R/Seurat versions, original
ASC files, the strict identity classifier training generator and the future
public mouse-data accession remain outstanding.
