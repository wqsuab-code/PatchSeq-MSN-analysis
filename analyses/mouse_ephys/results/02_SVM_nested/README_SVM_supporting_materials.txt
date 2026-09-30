SVM MAIN-PANEL SUPPORTING MATERIALS
===================================

Cohort
------
450 HC × GC-consensus cells from the active 493-cell cohort.
Class support: E1=142, E2=160, E3=45, E4=61, E5=42.
Predictors: 18 frozen electrophysiological features.

Main panel
----------
00_main_SVM_confusion_panel: standalone manuscript panel showing the row-normalized
nested out-of-fold confusion matrix.

Necessary supporting figures
-----------------------------
01a_confusion_raw_counts: raw nested OOF cell counts.
01b_confusion_row_normalized: percentages normalized within each true class.
02_cross_validation_workflow: leakage-controlled nested-CV design.
03_regularization_C_scan: repeated held-out performance across C values.
04_classwise_performance: precision, recall, specificity, F1 and class support.
05_one_vs_rest_ROC_PR: one-vs-rest ROC and precision-recall curves.
06_nested_OOF_probability_matrix: held-out Platt-scaled class probabilities.
07_misclassified_cells_fixed_tSNE: 36 nested-OOF errors on unchanged GC t-SNE coordinates.
08_repeated_CV_stability: 50 held-out folds from 10 repeated 5-fold CV runs.
09_label_permutation_test: 999-permutation null distribution.
10_signed_linear_SVM_feature_weights: signed one-vs-rest surrogate coefficients.

Associated data
---------------
Every figure has a PDF and 900-dpi PNG. CSV/JSON files contain the plotted values,
cell-level error list, fold-level scores, permutation null, AUC values, class metrics,
and coefficient means/SDs. The Python script reproducing the analysis is included.

Interpretation boundary
-----------------------
The target labels and predictors were derived from the same electrophysiological
feature space. These analyses quantify learnability and internal reproducibility;
they are not an independent biological validation of the taxonomy.

