# Complete methods: frozen Macaque M1–M4 morphology analysis

# Macaque MSN cohort and morphology quality control

Analyses were restricted to macaque medium spiny neurons annotated as STR D1 MSN, STR D2 MSN or STR Hybrid MSN and sampled from the caudate nucleus, putamen or nucleus accumbens. Cell metadata and morphology measurements were read directly from `Macaque-PatchSeq-BG.zip`, whose SHA-256 checksum was verified before analysis. Mouse morphology results and macaque electrophysiology features were not used.

The source metadata contained 486 eligible MSN profiles. Morphology classification required complete observations across 18 nonredundant dendritic and somatic measurements. Axonal variables, axonal Sholl measurements, the circular variable `axon_exit_theta_coronal`, and the derived `3_Sholl_PC1` score were excluded. A total of 126 cells from 42 donors met the complete-case criterion. No imputation was performed and no cell was removed solely because of an extreme PCA score.


# Morphology transformation and principal component analysis

The 18 retained morphology features were transformed independently. Adjusted sample skewness was calculated on the 126-cell cohort. Features with an absolute adjusted skewness of at least 0.5 underwent maximum-likelihood Yeo–Johnson transformation, with the transformation parameter optimized over −5 to 5; the remaining variables were left on their original scale. Each resulting feature was centered and divided by its sample standard deviation. Principal component analysis was then performed by singular-value decomposition without further centering or scaling. The frozen classification used PC1–PC5.

The transformation and PCA were fitted to all 126 morphology-complete cells because they were part of the unsupervised discovery procedure. In donor-held-out machine-learning validation, transformation, standardization and PCA parameters were instead estimated within each training split and applied unchanged to the corresponding test cells.


# Hierarchical and graph-based consensus classification

Hierarchical clustering was performed in PC1–PC5 space using Euclidean distance and Ward's minimum-variance linkage (Ward.D2), and the dendrogram was cut at four groups. In parallel, a shared-nearest-neighbour graph was constructed with 20 nearest neighbours and an SNN pruning threshold of 1/15. Louvain community detection used algorithm 1, resolution 2.3 and random seed 777, yielding 13 graph communities. Each graph community was mapped to the hierarchical cluster containing the largest number of its cells, producing four merged graph classes.

Cells with identical hierarchical and merged graph assignments were defined as the frozen consensus set. Agreement was observed for 117 of 126 cells (92.86%; adjusted Rand index 0.825), yielding M1=43, M2=42, M3=20 and M4=12. Nine discordant cells were retained in displays where indicated but excluded from consensus-class summaries and supervised label-recovery analyses. T class, anatomical region and donor identity were not used to construct the M labels.


# Low-dimensional displays, dendrogram, radar plots and heat map

For visualization, exact t-SNE was applied to PC1–PC5. Perplexities of 10, 15, 20, 25 and 30 and random seeds 777–781 were compared using 10-nearest-neighbour preservation and hierarchical-class silhouette width. The selected embedding used perplexity 25, seed 777 and 2,000 iterations. For the displayed compact version only, within-class deviations from each hierarchical-cluster centroid were multiplied by 0.85; clustering, labels and inferential results were unchanged. Dashed t-SNE boundaries are 80% bivariate-normal regions fitted to consensus cells.

The circular dendrogram displays the complete Ward.D2 tree for all 126 cells. Leaf order was taken directly from the tree. The inner ring shows hierarchical labels and the outer ring shows merged graph labels.

Radar plots use ten prespecified features drawn from the frozen 18-feature standardized matrix. The final feature order and labels are stored in the accompanying radar-layout JSON. All classes share the same cohort-wide Z-score scale from −3 to 3; values outside this range are clipped only for display. Pale curves show individual cells and black closed curves show class medians. No within-class normalization or interior fill was used.

The heat map shows all 18 frozen feature Z-scores for the 117 consensus cells, clipped to −2.5 to 2.5 for display. Cells were grouped by M class and shuffled within class using a fixed seed. Features were assigned to the class with the largest mean Z-score and ordered within each module by Ward.D2 clustering. Annotation tracks show M class, T class and region.


# T-class composition and transcriptome-to-morphology reduced-rank regression

Transcriptomic D1, D2 and hybrid identities were treated as post hoc annotations and were not used to define morphology classes. The 117-cell consensus cohort contained 47 D1, 62 D2 and 8 hybrid MSNs. T-class composition was summarized as counts and within-M-class percentages.

Transcriptome-to-morphology coupling was examined by pooled reduced-rank regression (RRR) in the 117 consensus cells. Transcriptomic predictors were the first 20 principal components derived from 1,000 eligible variable genes. Responses were the 18 frozen transformed and standardized morphology features; `3_Sholl_PC1` was excluded. A rank-three model was used for display. Panels show components 2 versus 1 and 3 versus 1 in transcriptomic and morphological score space. D1 and D2 cells were assigned distinct colors and 90% bivariate-normal ellipses; hybrid cells were retained as grey points. The ten genes or morphology features with the largest correlation-loading magnitudes were shown in each component plane, with vectors uniformly rescaled for display.


# Donor-grouped machine-learning validation

The frozen 117-cell M1–M4 consensus set was evaluated using 500 repeated donor-grouped analyses. Donors were kept intact so that cells from the same donor could not occur in both training and test sets. In every split or resample, feature transformation, standardization, PCA, feature selection and model fitting were estimated using training data only. Complete-case analysis was maintained without imputation.

Unsupervised discovery stability was assessed by repeatedly retaining approximately 80% of donors, refitting preprocessing and PCA, and comparing de novo Ward, k-means, Gaussian-mixture and spectral solutions with the frozen labels after optimal label matching. The scan covered 2–10 PCs and K=2–8. Supervised recoverability was evaluated on repeated 75%/25% donor-grouped partitions that contained all four classes in both subsets. Balanced accuracy and macro-F1 were the primary metrics. Label-permutation controls preserved donor structure. Cell-level stability was summarized from pairwise co-clustering probabilities, and feature robustness was evaluated using training-fold permutation importance, Extra Trees impurity importance, multinomial-logistic coefficients and progressive feature ablation. These analyses measure internal recoverability and resampling stability, not independent biological validation.


# Statistical analysis and interpretation

The main clustering figure is descriptive. HC–GC agreement was summarized by the proportion of identical assignments and the adjusted Rand index. T-class and regional compositions are reported as cell counts and percentages. Radar plots show cell-level profiles and class medians; heat-map values are standardized feature measurements. No null-hypothesis test was used to define M1–M4.

For the separate 18-feature class-comparison analysis, each feature was compared across M1–M4 using a two-sided Kruskal–Wallis test. The 18 omnibus P values were adjusted using the Benjamini–Hochberg procedure. Post hoc pairwise comparisons used Dunn tests with Holm correction. Cells were the displayed observational units; donor-grouped resampling was used for machine-learning validation to address within-donor dependence.
