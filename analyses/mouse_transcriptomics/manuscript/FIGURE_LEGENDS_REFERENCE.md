# Figure legends: working reference

These legends reflect the supplied draft and are not yet the final manuscript text. Numerical statements must be checked against `config/frozen_parameters.yml` and `docs/FREEZE_DECISIONS.md`.

## Fig. 1 | Reference-guided transcriptomic classification of Patch-seq medium spiny neurons

**a**, t-SNE representation of 20,029 reference medium spiny neurons, colored by eight D1 and eight D2 transcriptomic subtypes. Subtype labels were determined independently of the t-SNE visualization. **b**, Projection of 588 Patch-seq cells with stable MSN identities onto the reference embedding. Reference cells are shown in light gray and Patch-seq cells are colored by their modal Pearson-correlation subtype across 500 gene-bootstrap replicates. **c**, Chord diagram comparing Pearson- and RPCA-based subtype assignments. Pearson and RPCA classifications agreed for 422 of 588 cells at the exact subtype level and for 572 of 588 cells after collapse to D1 and D2. **d**, Expression distributions of MSN identity and subtype-associated genes across the 16 reference and Patch-seq populations. Reference and query counts were independently library-size normalized to 10,000 and transformed as natural-log `log1p`. Each gene was scaled to its maximum across the combined display data; the unscaled maximum is shown at right. Bottom bars show the subtype percentages and exact cell counts within each dataset.

## Extended Data Fig. 1 | Patch-seq workflow and transcriptomic quality control

The supplied layout contains the experimental workflow, detected-gene and MT/Ribo summaries for 643 sequenced libraries, two independent reference-correlation QC diagnostics for 641 analyzed cells, sample attrition, reads/genes and exon/intron relationships, and neural/non-neural marker expression. The final legend must retain the distinction between method-specific risk tiers and the joint 622 Good/19 Risk annotation.

## Extended Data Fig. 2 | Stability of major cell-type and MSN subtype assignments

The final legend must state the actual global RPCA scan values (10, 12, ..., 30 PCs), identify PC10 as exploratory and define final stability across 12, 14, ..., 30 PCs. The 500-bootstrap heatmap contains 339 High, 175 Moderate and 74 Low/Inconsistent stable MSN cells. Bootstrap support measures sensitivity to feature resampling and is not a calibrated classification probability.
