Macaque Ca+Pu+NAC four-panel transcriptomic/electrophysiological RRR display

Cells
- HC-GC consensus D1/D2 cells: n=346
- D1=163; D2=183
- E classes among these cells: C1=44, C2=130, C3=57, C4=115
- Hybrid cells are excluded because D1/D2 domains define the shared RRR model and ellipses.

Inputs
- Transcriptomic domain: log2 expression H5AD; technical/ribosomal/mitochondrial genes excluded;
  top 1,000 variable unique gene symbols; column Z-score; PCA up to 20 PCs.
- Electrophysiological domain: the final 19 transformed and Z-scored E features used for PCA/clustering.

Model and display
- Pooled T-to-E reduced-rank regression; rank=3 fitted for the shared latent basis.
- Panel order from left to right:
  1) T space, Component 2 on X / Component 1 on Y, points colored by D1/D2;
  2) T space, Component 3 on X / Component 1 on Y, points colored by D1/D2;
  3) E space, Component 2 on X / Component 1 on Y, points colored by C1-C4;
  4) E space, Component 3 on X / Component 1 on Y, points colored by C1-C4.
- All panels: D1/D2 covariance ellipses covering 90% under a bivariate-normal approximation.
- Arrows are correlation loadings: top 9 genes or top 10 E features by 2D vector magnitude.
- Cells are scaled by the 99th percentile radial distance for display inside the unit circle.
- Canvas: 5.6 x 1.0 inches; PNG 1200 dpi plus vector PDF.

Interpretation boundary
- Component 1 is the primary supported shared direction.
- Higher components are exploratory display directions and should not be described as independently validated axes.
