# Full ML stability and reproducibility audit of M1–M4

## Design

- Frozen cohort: 187 cells; 105 recording-day groups used as the animal proxy.
- Frozen predictors: 10 de-redundant morphology features.
- Primary target: HC-derived M1–M4 labels (M1=67, M2=23, M3=62, M4=35).
- Leakage control: every transformation is fitted on the training fold only; recording days never cross train/test folds.

## Results

- Nested repeated grouped OOF: accuracy 0.941, balanced accuracy 0.930, macro-F1 0.935, ARI 0.854.
- Nested leave-one-recording-day-out: accuracy 0.941, balanced accuracy 0.930, macro-F1 0.935, ARI 0.853.
- Thirty partition seeds: mean balanced accuracy 0.929 ± 0.008; mean pairwise prediction ARI 0.919.
- Label permutation: observed balanced accuracy 0.941; null mean 0.254; empirical P=0.0020 (500 permutations).
- Probability quality: macro OVR ROC-AUC 0.996, macro average precision 0.988, ECE 0.101.
- Integrated high-confidence reproducible cells: 157/187 under the prespecified HC–GC agreement, nested-OOF correctness, LODO correctness, OOF confidence ≥0.70 and 30-seed reference agreement ≥0.80 rule.

## Interpretation

These results test whether the four HC-derived morphology labels are recoverable in unseen recording-day groups and robust to model family, partition seed, feature removal and training-set size. They are internal reproducibility evidence. They do not independently establish that M1–M4 are discrete biological cell types because the labels and predictors originate from the same morphology measurements; external reconstruction cohorts, blinded manual annotation, spatial association or orthogonal transcriptomic/electrophysiological evidence remain independent validation layers.
