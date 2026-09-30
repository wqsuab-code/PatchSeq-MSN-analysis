# Frozen macaque MSN morphology M1-M4 and M-T analysis

This review repository contains the frozen macaque medium spiny neuron (MSN)
morphology analysis from source data through the transcriptome-to-morphology
(M-T) reduced-rank regression analysis. The scope is restricted to macaque MSN
cells from caudate, putamen and nucleus accumbens. Mouse morphology results and
macaque E-class features are not used as M-class predictors.

## Frozen analysis

- Source cohort: 486 macaque MSN profiles in Ca, Pu and NAc.
- Complete morphology cohort: 126 cells from 42 donors.
- Frozen features: 18 nonredundant dendritic/somatic measurements.
- Transformation: feature-wise Yeo-Johnson when absolute adjusted skewness was
  at least 0.5, otherwise direct use, followed by cohort-wide Z-scoring.
- Classification space: PC1-PC5.
- HC: Euclidean distance, Ward.D2, K=4.
- GC: SNN/Louvain, resolution 2.3, 13 raw communities merged to four by maximum
  overlap with HC.
- Consensus cohort: 117 cells from 41 donors; M1=43, M2=42, M3=20, M4=12.
- M-T RRR: 1,000 variable genes, 20 transcriptomic PCs, 18 morphology features,
  rank 3; D1, D2 and hybrid cells retained.

## Repository organization

- `data/`: source archive, QC tables, frozen inputs, transformed matrices, PCA
  scores and frozen labels.
- `code/`: exact analysis scripts organized by pipeline stage.
- `results/`: feature statistics, full 500-repeat ML validation and M-T RRR
  source data.
- `figures/`: main SVG, figure panels and vector/raster validation figures.
- `interactive/`: offline interactive morphology, radar and RRR figure editors.
- `docs/`: English and Chinese Methods and figure-legend materials.
- `audit/`: frozen-result consistency report, portability report and checksums.

## Important review status

This is a review package, not yet a public release. Exact historical scripts are
preserved, including some Windows absolute paths. See
`audit/PORTABILITY_AUDIT.md` before attempting a clean-machine rerun. Display-only
t-SNE contraction and boundary fitting do not affect clustering or statistical
testing.

## Large files

The source ZIP and large ML raw tables are included. `.gitattributes` marks them
for Git LFS. Install Git LFS before pushing this repository.

## Integrity

Run `python tests/test_archive_integrity.py` from the repository root. File-level
SHA-256 hashes are stored in `audit/SHA256SUMS.txt` and `audit/FILE_INDEX.csv`.

## License and citation

The scientific-data license and code license have not been selected. Resolve the
placeholders in `LICENSE_PENDING.md` and `CITATION.cff` before public release.
