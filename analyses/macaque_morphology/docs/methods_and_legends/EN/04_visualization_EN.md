# Low-dimensional displays, dendrogram, radar plots and heat map

For visualization, exact t-SNE was applied to PC1–PC5. Perplexities of 10, 15, 20, 25 and 30 and random seeds 777–781 were compared using 10-nearest-neighbour preservation and hierarchical-class silhouette width. The selected embedding used perplexity 25, seed 777 and 2,000 iterations. For the displayed compact version only, within-class deviations from each hierarchical-cluster centroid were multiplied by 0.85; clustering, labels and inferential results were unchanged. Dashed t-SNE boundaries are 80% bivariate-normal regions fitted to consensus cells.

The circular dendrogram displays the complete Ward.D2 tree for all 126 cells. Leaf order was taken directly from the tree. The inner ring shows hierarchical labels and the outer ring shows merged graph labels.

Radar plots use ten prespecified features drawn from the frozen 18-feature standardized matrix. The final feature order and labels are stored in the accompanying radar-layout JSON. All classes share the same cohort-wide Z-score scale from −3 to 3; values outside this range are clipped only for display. Pale curves show individual cells and black closed curves show class medians. No within-class normalization or interior fill was used.

The heat map shows all 18 frozen feature Z-scores for the 117 consensus cells, clipped to −2.5 to 2.5 for display. Cells were grouped by M class and shuffled within class using a fixed seed. Features were assigned to the class with the largest mean Z-score and ordered within each module by Ward.D2 clustering. Annotation tracks show M class, T class and region.
