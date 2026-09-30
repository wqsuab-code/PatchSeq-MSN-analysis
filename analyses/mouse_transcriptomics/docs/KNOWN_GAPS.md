# Known gaps and manuscript checks

## Must resolve before final submission

1. **ED Fig. 1a provenance**: register the original electrophysiology image, trace, SMART-Seq workflow schematic, fluorescence image and morphology reconstruction. The supplied composite alone is not enough for a reproducible source package.
2. **ED Fig. 1b-c assembly**: raw data sources are identified, but a clean script and exported plot-data table for the 643-cell panels still need to be created.
3. **Composite assembly scripts**: exact placement, lettering and cropping for all three supplied layouts are not currently scripted. Individual scientific panels are available.
4. **Sequencing processing fields**: adapter trimmer, aligner version/parameters, mouse genome build and annotation release remain manuscript metadata fields until confirmed from the primary sequencing pipeline.
5. **IN subtype analysis**: branch-specific subtype prediction remains pending and must not be represented as frozen.

## Already reconciled

- The query detected-gene median in the supplied QC layout is 9,394, not 3,394.
- The 500-replicate MSN bootstrap uses 721 draws with replacement per replicate, not a 90% draw.
- The Fig. 1 query subtype distribution is the modal Pearson distribution; it is not the RPCA distribution.
- The 422 exact agreements refer to Pearson versus RPCA subtype labels. A separate RPCA-centroid consensus analysis uses a different definition and must not be substituted.
- ED Fig. 2a-c has now been rebuilt with the actual `10,12,...,30` nPC axis. The supplied draft layout remains preserved only as a historical visual reference.
