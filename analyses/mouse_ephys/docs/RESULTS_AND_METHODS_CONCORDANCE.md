# Results and Methods concordance

The final R 4.5.1 and Python 3.14.3 rerun reproduced the frozen preprocessing,
18-feature matrix, PCA, HC and GC assignments, GC-HC consensus mask, E1-E5
counts, fixed t-SNE coordinates, univariate statistics, nested-SVM outputs,
recording-date-grouped machine-learning outputs and strict-441 RRR tables. All
24 prespecified concordance checks passed.

Numbers approved for Methods, Results and legends:

- 549 cells entered E-QC; 494 passed complete-case QC; 493 entered taxonomy.
- 450/493 cells were GC-HC consensus (91.28%).
- E1-E5 counts: 142, 160, 45, 61 and 42.
- Strict identity: 441 cells, comprising 209 D1 and 232 D2 cells; nine of the
  450 consensus cells were ambiguous.
- Strict E1-E5 counts: 140, 155, 45, 59 and 42.
- Nested linear-SVM OOF accuracy 0.9200, balanced accuracy 0.9320 and macro-F1
  0.9205; selected C = 0.1.
- Recording-date-grouped robustness: 200 splits; mean balanced accuracy 0.8678,
  median 0.8705; aggregated cell-level accuracy 411/450 (91.33%).
- T-E RRR: 1,000 genes, 20 transcriptomic PCs, 18 E-features, rank 3,
  in-sample multivariate R² = 0.2958.

The 92.0% and 91.33% values derive from different validation designs. Neither is
a measure of unsupervised GC-HC agreement. The RRR R² is descriptive and must
not be described as held-out prediction.
