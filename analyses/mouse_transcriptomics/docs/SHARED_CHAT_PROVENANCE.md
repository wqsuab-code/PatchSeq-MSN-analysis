# Shared-chat provenance: Ref IN MSN Mapping

## Source

- Shared conversation: `https://chatgpt.com/share/6abc9d75-4fb0-83ea-b358-c04be90ffb24`
- Page title: `Ref IN MSN Mapping`
- Reviewed: 2026-09-30
- Scope visible in the shared copy: approximately 246 user prompts covering reference reconstruction, broad-class QC, MSN/IN mapping exploration, plotting and code revisions.

The shared conversation is retained as historical provenance. It is not itself an analysis result and does not supersede the frozen tables, scripts or checksums in this repository.

## Source paths recovered from the history

The conversation records the following server-side inputs used in an earlier broad-class QC/bootstrap workflow:

- Reference: `/home/qw4/Project_NAc_Mouse_patch-seq/ForThePaper/Transcriptomic mapping_20260626/00_input/YZ_all_final_fixed.rds`
- Query: `/home/qw4/Project_NAc_Mouse_patch-seq/ForThePaper/Transcriptomic mapping_20260626/00_input/query_BroadClass_FrozenHVG_gene90_ref90_500iter_query_with_QCmapping.rds`
- Earlier local QC table candidate: `C:/Users/53461/OneDrive/Desktop/Patch-seq/keyfile/Nature_supplement_QC/STAR_alignment_QC_cleaned.csv`

These locations are historical pointers. Availability and file identity must be verified before reuse. The authoritative local inputs for the frozen package are listed in `RAW_SOURCE_REGISTRY.md` with file sizes and SHA-256 hashes.

## Historical broad-class QC workflow

The shared conversation contains a proposed broad-class robustness analysis with these settings:

- Reference label column: `CellType`.
- Reference-only variable-gene ranking; first 2,000 ranked genes present in Query.
- 100 iterations.
- 90% stratified Reference-cell sampling within each broad class.
- Seurat anchor transfer and SingleR run in parallel as two classifiers.
- 20 PCs and `k.anchor = 5` in the recorded script.
- Per-cell outputs included consensus label frequency, entropy-derived stability, top-to-second margin and Seurat-SingleR agreement.

This workflow concerns broad-class QC and sensitivity to Reference-cell resampling. It is distinct from the frozen MSN subtype procedure.

## Separation from the frozen MSN subtype result

The formal MSN subtype analysis in this package uses:

- 588 globally stable MSN Query cells.
- A Top750 Reference-only dispersion-z ranking, yielding 721 usable shared features.
- RPCA with 16 PCs, `k.anchor = 10`, `k.filter = 200` and `k.weight = 30`.
- Pearson centroid assignment with 500 gene-bootstrap replicates.
- In each Pearson replicate, 721 genes are sampled with replacement from the complete 721-gene set.

Therefore, the historical 90% Reference-cell bootstrap must not be described as the 500-replicate MSN subtype bootstrap, and its Seurat-SingleR agreement must not be substituted for the frozen Pearson-RPCA agreement statistics.

## Additional design decisions documented in the history

- Reference subtype colors should remain fixed across t-SNE, UMAP, violin and agreement figures.
- Broad classes use D1 red, D2 blue and IN green.
- Reference cells are colored and Query cells are black when displaying Query-over-Reference embeddings for broad-class stability.
- Query and Reference distributions should use matched transparency in expression-distribution panels.
- IN parameter selection was recognized as branch-specific and should not copy MSN HVG or PC settings without independent validation.

## Verification status

- The shared page and the passages summarized above were inspected directly.
- Historical server paths have not been asserted to exist on this Windows workstation.
- No numerical result from the shared conversation is treated as frozen unless it is independently represented in the package data and validation checks.
- IN subtype mapping remains pending, as recorded in `KNOWN_GAPS.md`.
