Mouse HC-GC consensus E-class machine-learning robustness audit
================================================================

Cohort
------
450 HC-GC consensus cells with the frozen 18-feature matrix:
E1 = 142, E2 = 160, E3 = 45, E4 = 61, E5 = 42.

Grouping and leakage control
----------------------------
The experimental grouping variable was the recording date extracted from the
structured cell key. All cells recorded on the same date were kept together in
either training or test data. The analysis contains 160 experimental dates,
with 1-7 cells per date. The Batch2-Batch8 suffix was retained only as the
sequencing-submission batch metadata and was not used as the train-test grouping
variable.

Two hundred group-aware 75:25 train-test splits were generated. Every accepted
split contained E1-E5 in both train and test sets. Within each split, feature
minima and shifted column sums were estimated from training cells only. Values
were shifted to the training-fold floor, column-sum normalized to 10,000,
log1p transformed and standardized with training-fold means and sample standard
deviations. PCA (3 components) was then fitted in the training fold, and all
fitted operations were applied unchanged to the held-out cells.

Supervised model comparison
---------------------------
Seven classifiers were evaluated: Extra Trees, distance-weighted k-nearest
neighbors, RBF SVM, linear SVM, class-weighted logistic regression, random
forest and histogram gradient boosting. The RBF SVM had the highest median
held-out balanced accuracy.

RBF SVM results across 200 date-grouped splits:
- balanced accuracy: mean 0.8678, s.d. 0.0636, median 0.8705
- macro-F1: mean 0.8555, s.d. 0.0681, median 0.8607
- repeated held-out probability aggregation: 411/450 correct (91.33%)

Aggregated confusion counts (rows: frozen HC-GC consensus class; columns:
predicted class):
           E1   E2   E3   E4   E5
E1        136    4    0    1    1
E2         10  146    0    4    0
E3          0    2   41    2    0
E4          1    6    4   47    3
E5          1    0    0    0   41

Additional robustness analyses
------------------------------
- Label-permutation control: 100 experimental-date-block permutations, each
  evaluated on five fixed date-grouped splits; empirical P = 0.0099.
- Unsupervised discovery: Ward, k-means, GMM and spectral clustering were
  assessed across 200 date-resampled training sets.
- Full-data agreement with the frozen E labels (ARI): Ward 0.625, k-means
  0.290, GMM 0.283 and spectral 0.306.
- Metadata association: T class, bias-corrected Cramer's V = 0.341,
  experimental-date-block permutation P = 0.0010; sequencing-submission batch,
  V = 0.109, P = 0.0100.
- Feature-rank agreement between Extra Trees MDI and absolute multinomial
  logistic coefficients: Spearman rho = 0.767.

Interpretation
--------------
The frozen HC-GC E1-E5 partition is recoverable from held-out electrophysiology
when same-day cells are excluded from the corresponding training set. The
drop from the previously reported cell-random cross-validation estimate to the
date-grouped estimate is expected and provides a stricter, lower-leakage test.
The 91.33% aggregated cell-level value is based on averaging repeated held-out
probabilities and must not be reported as the median split-level balanced
accuracy (87.05%).

Reproducibility
---------------
Random seed: 20260916.
Analysis script: scripts/validate_mouse_E5_macaqueM_style.py
Interactive data: sites/ephys-core18-explorer/dist/data/mouse_ml_validation.json
Interactive editor: sites/ephys-core18-explorer/dist/ml-validation-editor.html
All panel-level source tables are the numbered CSV files in this directory.
