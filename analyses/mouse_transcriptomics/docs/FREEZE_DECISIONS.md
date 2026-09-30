# Frozen definitions and non-negotiable distinctions

## Cohort definitions

- `sequenced`: all 643 transcriptomic libraries.
- `mapped`: 641 libraries after removal of O9268_Batch2 and O9273_Batch2 for extremely low read depth.
- `Good`: retained by both the Pearson-centroid and five-PC QC diagnostics.
- `Risk`: not retained by one or both QC diagnostics. Risk is an annotation, not an exclusion from mapping.
- `stable MSN`: `Final_Root_stability_from_CellType == "Stable_MSN"` (n = 588).
- `stable IN`: `Final_Root_stability_from_CellType == "Stable_IN"` (n = 34).

## Global RPCA stability

The source scan contains 11 actual nPC settings: 10, 12, 14, ..., 30. PC10 is exploratory. The final CellType and root trajectories contain ten values and are explicitly named `PC12_PC30`; they therefore use 12, 14, ..., 30.

Any figure axis that labels these settings as consecutive `PC10` to `PC20`, or any text that describes the retained settings as `11-20 PCs`, must be corrected or explicitly relabeled as scan indices. The numerical nPC values must not be silently changed.

## MSN subtype reporting

- RPCA is the frozen transfer method: Top750 reference-only dispersion ranking, 721 effective genes, PC16, `k.anchor=10`, `k.filter=200`, `k.weight=30`.
- Pearson is an independent centroid-based support analysis using the same 721-gene feature pool.
- The formal 500-replicate run samples 721 genes with replacement in every replicate. It is a conventional full-size bootstrap, not a 90% subsample.
- The modal Pearson subtype is used in Fig. 1b and Fig. 1d.
- RPCA-Pearson exact agreement is 422/588. D1/D2 agreement is 572/588.
- Bootstrap support is sensitivity to feature resampling, not a calibrated posterior probability.

## Confidence tiers

- High: modal support at least 0.80 and top1-top2 vote margin at least 0.50.
- Moderate: modal support at least 0.50 and margin at least 0.20, excluding High.
- Low/Inconsistent: all remaining cells.

The manuscript may rename High/Moderate/Low as Highly consistent/Moderately consistent/Inconsistent, but the numerical thresholds and membership must remain unchanged.

## Visualization

The locked manuscript t-SNE uses reference PC1-PC16, perplexity 15 and learning rate 100. This differs from the earlier visual-quality scan recommendation of perplexity 50 and learning rate 500. Because the supplied submission figure uses the 15/100 layout, the present framework records 15/100 as the display freeze while keeping subtype labels independent of this choice.

## Color mapping

The 16 subtype colors are fixed in `data/msn/MSN_fixed_subtype_color_mapping.csv`. Broad D1/D2/IN colors remain red/blue/green. Subtype plots must not use gray or black for subtype identities; gray and black are reserved for reference context, outlines or QC status.

## IN subtype branch

IN subtype mapping is not frozen. Its HVG number, feature set, PC count and mapping method must be optimized independently because the available reference and query anchors are much fewer than in the MSN branch.
