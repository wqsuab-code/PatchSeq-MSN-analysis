Frozen Macaque MSN morphology M1-M4: 500-repeat stability and reproducibility validation
======================================================================================
Scope: Macaque MSN only, Ca+Pu+NAC. Of 126 morphology-complete cells, 117 HC-GC consensus cells entered machine learning (M1=43, M2=42, M3=20, M4=12; 41 donors).

Leakage control: transformation, scaling, and PCA were re-estimated within every training or resampled donor set. Donors never crossed training and test sets. Frozen workflow reproduction: ARI=1.000, NMI=1.000.

Interpretation boundary: supervised performance measures recoverability of frozen labels from the same morphology domain; donor-resampled clustering measures de novo discovery stability; the consensus margin measures cell-level assignment stability; E/T/ROI associations are orthogonal annotations or confounding checks, not classifier performance. None alone establishes naturally discrete biological types.

Availability: consensus E4 annotation was available for 91/117 cells. The ZIP contains no independent reconstruction-batch or recording-date field; these were marked untestable rather than inferred. Deep learning was not run because n=117 and M4=12 are inadequate for primary neural-network evidence.
