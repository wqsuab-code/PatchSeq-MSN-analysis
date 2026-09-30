# Environment record

The frozen Python reporting environment is recorded in
`requirements-frozen.txt` and the RRR reproduction audit. It used Python 3.14.3,
NumPy 2.4.6, pandas 2.3.3, SciPy 1.17.1, scikit-learn 1.8.0, joblib 1.5.3,
Matplotlib 3.10.9, seaborn 0.13.2 and NeuroM 4.0.5.

An exact R/Seurat lock file and `sessionInfo()` were not present in the frozen
archive. The HC implementation is reproducible in Python/SciPy; exact
re-execution of the Seurat graph-clustering branch requires recovery of the
original R/Seurat environment. This is a documented portability gap, not a
change to the frozen labels.
