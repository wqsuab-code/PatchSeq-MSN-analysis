# Mouse M morphology taxonomy and T-M analysis

Audited frozen module for mouse nucleus accumbens medium spiny neuron morphology.
The primary taxonomy contains 187 cells (M1=67, M2=23, M3=62, M4=35), of
which 181 show HC-GC agreement. Cross-modal T-M reduced-rank regression (RRR)
uses only 168 HC-GC-consensus cells with strict RPCA-Pearson-concordant D1/D2
identity (D1=74, D2=94).

## Frozen primary analysis

Ten de-redundant morphology features were transformed by feature-wise minimum
shift, shifted column-sum normalization to 10,000, `log1p`, and Z-scoring.
PC1-PC3 were used for Euclidean Ward.D2 hierarchical clustering at K=4.
Independent graph clustering used a Seurat SNN graph (k=20, prune=1/15),
Louvain algorithm 1, and resolution 0.50. M1-M4 are the frozen HC-derived
labels; HC-GC agreement defines the 181-cell consensus subset.

## Scope and interpretation

- `n=187`: complete frozen morphology taxonomy and internal ML recovery.
- `n=181`: HC-GC morphology-consensus subset.
- `n=180`: broad descriptive D1/D2 identity among consensus cells; not an RRR input.
- `n=168`: strict RPCA-Pearson D1/D2 and HC-GC consensus; authoritative T-M RRR cohort.
- Recording day is used as a mouse proxy because true animal IDs are unavailable
  and one mouse was used per day.
- ML evaluates recoverability and stability of a taxonomy derived from the same
  morphology features; it is not an independent biological validation cohort.

## Directory guide

- `config/`: frozen feature order and analytical parameters.
- `code/`: ordered QC, preprocessing, classification, statistics, ML, RRR,
  interactive and packaging code.
- `data/`: source manifests, cohort audit, frozen inputs, PCA products and labels.
- `results/`: clustering scan, bootstrap, ML, ASC audit and strict RRR results.
- `figures/`: final composite and source publication figures.
- `interactive/`: local feature, RRR and ML editors with embedded frozen data.
- `docs/`: Methods and reusable Methods sections.
- `environment/`: recorded Python environment and missing R-lock disclosure.
- `tests/`: fast frozen-integrity checks runnable in a clean clone.
- `audit/`: version, conflict, path, data-source and release-asset audits.
- `legacy/`: explicitly non-primary historical configurations only.

## Reproduction status

The frozen PCA and Ward.D2 labels were reproduced exactly (ARI=1.0; score
differences at floating-point precision). Strict T-M RRR outputs were reproduced
to floating-point precision in the recorded Python environment. Saved ML run
tables, split logs, parameters and summaries are verified here; the expensive
500-repeat suite is not rerun by the fast archive test. The standalone generator
for the strict transcriptomic identity classifier has not been recovered, so its
frozen predictions are auditable but cannot yet be refitted end-to-end.

Original ASC files are not redistributed or embedded. See
`audit/RAW_ASC_acquisition_status_for_final187.csv` and
`data/00_source_manifest/source_manifest.csv`.

Run `python tests/test_archive_integrity.py` from this directory for the minimum
clean-clone verification, and consult `RUN_ORDER.md` for the complete sequence.
