# Portability and provenance audit

## Portable components

- Frozen analysis-ready inputs, labels, statistical tables, ML run tables and
  RRR scores are stored with relative paths under this module.
- `tests/test_archive_integrity.py` runs without the original project tree.
- `code/02_preprocessing/reproduce_frozen_preprocessing_pca_hc.py` independently
  reproduces the archived transform, PCA and Ward K=4 partition.
- Third-party raw data are represented by accessions, byte sizes and SHA-256
  values rather than redistributed files.

## Known nonportable historical components

- Several archived historical scripts retain absolute Windows paths. They are
  preserved as provenance and require path configuration before reuse.
- The exact historical R script and `sessionInfo()` that generated the final
  Seurat GC15 result were not recovered. Frozen assignments, merge maps, seed
  scans and confusion tables are archived, but this one step is not claimed as
  clean-clone executable.
- Interactive builders may refer to the original project tree; the built static
  site and layout JSON files are archived separately.

## Scientific boundary

No Mouse data, macaque morphology variables, alternate E labels, or Yeo–Johnson
primary preprocessing are used in the frozen Macaque E classification. Any
Yeo–Johnson branch is a sensitivity analysis and is not the primary result.
