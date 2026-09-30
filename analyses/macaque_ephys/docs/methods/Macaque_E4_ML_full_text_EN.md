# Machine-learning validation of the frozen macaque E1–E4 electrophysiological classification

## Proposed figure title

Donor-aware machine-learning validation supports the recoverability of four frozen macaque electrophysiological classes while revealing moderate de novo clustering stability.

## Figure legend

**Figure X | Donor-aware validation, stability and feature dependence of the frozen macaque E1–E4 electrophysiological classification.** The analysis included medium spiny neurons from caudate (Ca), putamen (Pu) and nucleus accumbens (NAc). Unsupervised analyses used all 390 cells with complete measurements for 19 electrophysiological features from 55 donors; supervised analyses used the 368 cells from 54 donors for which Ward hierarchical clustering and merged graph-based clustering agreed (E1, n=59; E2, n=133; E3, n=57; E4, n=119). Unless otherwise stated, resampling was performed 500 times at the donor level, and preprocessing parameters were estimated without access to held-out donors.

**a,b,** Balanced accuracy (**a**) and macro-averaged F1 score (**b**) for seven classifiers evaluated in 500 donor-held-out splits. Each split assigned 25% of donors to the test set and was accepted only when all four classes occurred in both training and test partitions. Boxes show the median and interquartile range; whiskers extend to 1.5 times the interquartile range. The dotted line marks chance-level balanced accuracy for four classes (0.25). ET, Extra Trees; kNN, k-nearest neighbours; RBF, radial-basis-function support-vector machine; Linear, linear support-vector machine; Logit, multinomial logistic regression; RF, random forest; HGB, histogram-based gradient boosting.

**c,** Row-normalized confusion matrix for Extra Trees predictions aggregated across the 500 donor-held-out splits. Values denote the percentage of predictions within each frozen class.

**d,** Matched permutation control for Extra Trees. The complete training-fold preprocessing, three-component PCA and classifier pipeline was repeated after independently shuffling E-class labels within each donor. Grey bars show the distribution of mean fivefold donor-grouped balanced accuracy across 500 permutations. The orange line shows the observed mean balanced accuracy (0.982); the null median was 0.304 and the empirical permutation P value was 0.0020.

**e,f,** One-versus-rest receiver-operating-characteristic (**e**) and precision–recall (**f**) curves derived from cell-level Extra Trees probabilities. For each cell, probabilities were averaged across every donor-held-out split in which that cell appeared in the test partition. The legends report area under the ROC curve (AUC) and average precision (AP), respectively. The dashed diagonal in **e** denotes random discrimination.

**g,** Sensitivity of Ward clustering to the number of retained principal components (NPC, 2–10) and requested clusters (K, 2–8). For each parameter combination, the heat map shows the median adjusted Rand index (ARI) relative to the frozen four-class labels over 500 resamples containing 80% of donors. The white box marks the frozen NPC=3, K=4 setting. The scan is a sensitivity analysis rather than a new optimization of the frozen solution.

**h,** Recovery of the frozen four-class partition by Ward, K-means, full-covariance Gaussian mixture modelling (GMM) and nearest-neighbour spectral clustering after 80% donor resampling. Violin distributions represent 500 resamples; internal boxes summarize their central distributions.

**i,** Full-data agreement between each unsupervised algorithm and the frozen Ward partition, quantified by ARI. Ward reproduces the frozen partition by construction and serves as the value of 1.0.

**j,** Cell-level co-clustering margin calculated from donor-resampled Ward solutions. For each cell, the margin is the mean probability of co-clustering with cells in its frozen class minus the largest mean probability of co-clustering with any other frozen class. Positive values support the frozen assignment; the dashed line marks zero.

**k,** Co-clustering margins stratified by frozen class. Violin plots show the distributions and embedded boxes show their central tendency.

**l,** Association of the frozen E class with donor, transcriptomic class and anatomical region, measured using bias-corrected Cramér's V. Significance was assessed with 5,000 structured permutations: a global E-label shuffle for donor and within-donor E-label shuffling for transcriptomic class and region.

**m,** Leakage-safe permutation importance calculated in 25 donor-grouped test folds (five repeats of fivefold StratifiedGroupKFold). A random forest containing 500 trees was fitted after estimating all transformations in the training fold; each feature was permuted 20 times in the held-out fold. Bars show the mean reduction in balanced accuracy and error bars show the standard deviation across grouped folds and permutations.

**n,** Stability of Extra Trees mean-decrease-in-impurity (MDI) importance across 500 donor-held-out splits. Bars show median importance and error bars show the empirical 2.5th–97.5th percentile interval.

**o,** Agreement between feature ranks obtained from Extra Trees MDI and absolute multinomial logistic-regression coefficients across donor-held-out splits. The dashed line denotes identical ranks; Spearman rho=0.832, P=1.02×10^-5.

**p,** Nested progressive feature ablation. In each of 25 outer donor-grouped folds, features were re-ranked using threefold grouped cross-validation restricted to the outer training data. A 500-tree random forest was then evaluated on the untouched outer test donors using all 19 features, all except the two highest-ranked features, all except the five highest-ranked features, or only the three lowest-ranked features. The dotted line marks chance-level balanced accuracy (0.25).

## Results

### Frozen-class definition and exact computational recovery

The electrophysiological analysis was restricted to 390 macaque medium spiny neurons from Ca, Pu and NAc that had complete measurements for 19 nonredundant electrophysiological features. The frozen classification had been defined by Ward hierarchical clustering of PC1–PC3 together with an independently derived graph-clustering solution. Agreement between the four Ward clusters and the graph clusters after maximum-overlap merging retained 368 consensus cells (94.36%), comprising E1 (n=59), E2 (n=133), E3 (n=57) and E4 (n=119). Reapplication of the frozen transformations, PCA and Ward K=4 procedure recovered the stored labels exactly (ARI=1.000).

### Frozen E classes are recoverable in unseen donors

Across 500 donor-held-out splits, all seven supervised classifiers recovered the frozen labels with high accuracy. Extra Trees performed best, with a median balanced accuracy of 0.973 (95% empirical interval, 0.904–1.000) and median macro-F1 of 0.975 (0.914–1.000). kNN and the RBF support-vector machine showed similarly high median balanced accuracies of 0.971 and 0.968, respectively, indicating that recoverability was not confined to a single model family. Aggregate Extra Trees predictions were strongly diagonal: row-normalized recalls were 98%, 99%, 94% and 97% for E1–E4, respectively. Across individual splits, median class recall was 1.000 for E1, 1.000 for E2, 0.952 for E3 and 0.972 for E4.

The matched permutation test confirmed that performance could not be explained by donor-specific class composition alone. When labels were shuffled within donors and the full Extra Trees pipeline was rerun, the null median balanced accuracy was 0.304 (95% interval, 0.250–0.360), compared with an observed mean of 0.982 (empirical P=0.001996). Cell-averaged held-out probabilities yielded near-perfect one-versus-rest discrimination: ROC AUC values were 1.000, 1.000, 0.999 and 0.999 for E1–E4, and the corresponding AP values were 1.000, 1.000, 0.997 and 0.999. These estimates quantify internal donor-grouped recoverability and should not be interpreted as validation in an independent external cohort.

### De novo discovery is reproducible but less invariant than supervised recovery

The unsupervised analyses produced a more qualified result. Under 80% donor resampling, median ARI relative to the frozen E4 partition was 0.476 for Ward (95% interval, 0.253–0.776), 0.570 for K-means (0.430–0.760), 0.320 for GMM (0.160–0.410) and 0.719 for spectral clustering (0.520–0.820). On the complete dataset, ARI relative to frozen Ward clustering was 0.593 for K-means, 0.401 for GMM and 0.742 for spectral clustering. Thus, the four labels are highly learnable once defined, but the identical partition is not an algorithm-invariant outcome of every de novo clustering analysis.

The NPC-by-K scan further showed that the frozen NPC=3, K=4 setting was defensible but was not the numerical maximum of the sensitivity surface: its median resampled ARI was 0.48, whereas NPC=3, K=5 reached 0.54. The K=5 full-data solution contained a one-cell cluster, however, and therefore did not provide a preferable stable biological partition. The frozen K=4 solution should consequently be described as a prespecified consensus classification supported by complementary clustering and cluster-size constraints, not as a universally optimal value selected from the heat map.

Cell-level co-clustering margins were positive for most cells. Median margins were 0.44, 0.56, 0.61 and 0.44 for E1–E4, respectively; 16 of 390 cells had negative margins (E1, 4; E2, 10; E3, 0; E4, 2). These cells identify assignment-sensitive observations and explain why de novo resampling agreement is lower than supervised label recovery.

### Classification depends on distributed adaptation and waveform features

Leakage-safe permutation importance identified AHP delay at the fifth spike as the strongest individual contributor to held-out performance (mean decrease in balanced accuracy, 0.075). The next four features were upstroke adaptation ratio (0.054), width adaptation ratio (0.053), fifth-spike AHP-delay ratio (0.046) and the upstroke/downstroke ratio at rheobase (0.030). Extra Trees MDI rankings were stable across 500 held-donor splits and were led by the same broad group of adaptation-related variables. Feature rankings from Extra Trees and multinomial logistic regression were strongly correlated (Spearman rho=0.832), supporting a signal that was not specific to a single importance definition.

Nested ablation indicated that classification information was distributed rather than restricted to one or two variables. Median balanced accuracy was 0.891 (95% empirical interval, 0.788–0.954) with all 19 features and remained 0.895 (0.756–0.935) after removal of the two highest-ranked training-fold features. Performance decreased to 0.776 (0.710–0.863) after removing the five highest-ranked features and to 0.535 (0.436–0.638) when only the three lowest-ranked features were retained. The absence of a decline after removing two features, together with the clear loss after removing five, is consistent with partially redundant information distributed across a compact electrophysiological feature set.

### Metadata structure is present but does not account for donor-held-out recoverability

Bias-corrected Cramér's V indicated associations between E class and transcriptomic class (V=0.329, within-donor permutation P=0.00020), donor (V=0.247, global permutation P=0.00020) and ROI (V=0.164, within-donor permutation P=0.00160). These effects demonstrate that the E classes are biologically and sampling-structure aware rather than metadata independent. However, strong performance in held-out donors and the within-donor label-permutation result show that donor identity alone does not explain the recoverability of the frozen labels. The metadata associations remain observational and do not establish causal or independent biological validation.

### Overall interpretation

Together, these analyses support E1–E4 as a reproducible frozen operational classification of macaque MSN electrophysiology. The principal evidence is the high recovery of class labels in unseen donors, the failure of matched within-donor permutations to reproduce that accuracy, cross-model convergence on similar predictive features and the graded loss of performance during nested feature ablation. The more moderate stability of de novo clustering is an important boundary condition: the data support a useful and reproducible four-class reference system, but not the stronger claim that four sharply separated natural groups must emerge under every donor sample, clustering algorithm or parameter choice.

## Methods

### Analysis population and frozen labels

The analysis included macaque Patch-seq medium spiny neurons with anatomical region labels Ca, Pu or NAc. Cells from other regions were excluded. Strict complete-case filtering across the frozen 19-feature panel retained 390 cells from 55 donors; no missing values were imputed. The frozen reference labels were derived using PC1–PC3, Euclidean Ward clustering cut at K=4, and a Seurat shared-nearest-neighbour graph generated with `k.param=20`, `prune.SNN=1/15`, Louvain `algorithm=1`, `resolution=3.0` and `random.seed=777`. Fifteen raw graph clusters were mapped surjectively onto the four Ward clusters by maximizing overlap. Cells with identical Ward and merged graph labels were retained for supervised validation (n=368).

### Electrophysiological features and preprocessing

The 19 frozen features were: width at rheobase; fast-trough voltage at rheobase; peak delta voltage at rheobase; peak voltage at rheobase; post-action-potential slope at rheobase; threshold voltage at rheobase; trough time at rheobase; trough voltage at rheobase; upstroke/downstroke ratio at rheobase; fifth-spike AHP delay; fifth-spike AHP-delay ratio; post-action-potential slope at the hero sweep; trough time at the hero sweep; downstroke adaptation ratio; peak-voltage adaptation ratio; threshold-voltage adaptation ratio; upstroke adaptation ratio; width adaptation ratio; and threshold voltage from the short-square protocol.

For each analysis fold or donor resample, preprocessing was fitted anew. Strictly positive ratio variables were log2-transformed. For every other feature, the training-set minimum was subtracted, the shifted training column was scaled to a sum of 10,000, and values were transformed as log(1+x). Each transformed feature was standardized using the training-set mean and population standard deviation (`ddof=0`). Held-out data were transformed using the corresponding training parameters; values below a training-set minimum were mapped to zero before log transformation. No donor, ROI or transcriptomic labels were regressed from the electrophysiological measurements.

### Principal-component and clustering sensitivity analyses

PCA used full singular-value decomposition. The frozen clustering used PC1–PC3, which explained 54.62% of total variance. For the NPC-by-K sensitivity analysis, 80% of donors were sampled without replacement in each of 500 repetitions. Preprocessing and a ten-component PCA were fitted independently within each resample. Ward clustering was evaluated for NPC=2–10 and K=2–8, and agreement with the corresponding frozen labels for included cells was quantified by ARI.

For cross-algorithm discovery, four K=4 algorithms were evaluated on the first three PCs in each 80% donor resample: Ward agglomerative clustering; K-means with 20 initializations; a full-covariance GMM with five initializations and `reg_covar=1×10^-5`; and spectral clustering using a 15-nearest-neighbour affinity graph, k-means label assignment and ten initializations. Cell-level co-clustering probabilities were computed as the number of resamples in which a cell pair received the same Ward label divided by the number of resamples in which both cells were observed.

### Donor-held-out supervised validation

Supervised validation was restricted to the 368 frozen consensus cells. Five hundred `GroupShuffleSplit` partitions were generated with donor as the grouping variable and `test_size=0.25`; a split was retained only if all four classes were represented in both partitions. All feature transformations and the three-component PCA were fitted using training cells only and then applied to test cells.

Seven classifiers were tested: multinomial logistic regression (`max_iter=3000`, `class_weight='balanced'`); linear and RBF support-vector machines (`C=1`, `class_weight='balanced'`, default scaled gamma for the RBF model); random forest (200 trees, balanced class weights); Extra Trees (200 trees, balanced class weights); distance-weighted kNN (`n_neighbors=7`); and histogram-based gradient boosting (`max_iter=180`, `learning_rate=0.06`). Balanced accuracy, macro-F1, Matthews correlation coefficient and per-class recall were calculated in each held-out split. Reported intervals are empirical 2.5th and 97.5th percentiles across the 500 splits.

### ROC and precision–recall analysis

Extra Trees class probabilities were stored for every held-out prediction. For each cell, probabilities were averaged over all splits in which its donor was held out. One-versus-rest ROC and precision–recall curves were then calculated from these cell-level mean probabilities. ROC AUC and AP were calculated separately for each class. Because these predictions arise from repeated internal donor-grouped validation, they are not estimates from an external validation cohort.

### Matched label-permutation control

A fivefold `StratifiedGroupKFold` partition was selected such that all four classes were present in every training and test fold. The observed pipeline used training-fold preprocessing, three-component PCA and Extra Trees with 200 trees and balanced class weights. For each of 500 null repetitions, labels were independently permuted within each donor, and the identical folds, transformations, PCA and classifier family were rerun. The test statistic was mean balanced accuracy across the five folds. The empirical one-sided P value was calculated as `(1 + number of null values greater than or equal to observed)/(1 + number of permutations)`.

### Feature importance and nested ablation

Leakage-safe permutation importance used five repeats of fivefold `StratifiedGroupKFold`, giving 25 donor-grouped outer test folds. Transformations were fitted within each training fold. A random forest with 500 trees, `max_features='sqrt'`, `min_samples_leaf=1` and `class_weight='balanced_subsample'` was fitted to the transformed 19-feature matrix without PCA. Each feature was permuted 20 times in the corresponding held-out fold, and importance was defined as the reduction in balanced accuracy.

MDI importance was obtained from 200-tree Extra Trees classifiers fitted to the transformed 19-feature training matrix in each of the 500 donor-held-out splits. Absolute coefficients from class-balanced multinomial logistic regression were averaged across classes, and median feature ranks across splits were compared with Extra Trees ranks using Spearman correlation.

For nested ablation, the same 25 outer donor-grouped folds were used. Within each outer training set, features were ranked by random-forest permutation importance using threefold grouped cross-validation and ten permutations per feature. The ranking was therefore independent of the outer test donors. Separate 500-tree random forests were evaluated using all 19 variables, all except the two highest-ranked variables, all except the five highest-ranked variables, or the three lowest-ranked variables. Balanced accuracy and macro-F1 were evaluated only in the untouched outer test fold.

### Metadata association tests

Associations between frozen E class and donor, transcriptomic class or ROI were quantified using bias-corrected Cramér's V. Significance was assessed using 5,000 permutations. For the donor association, E labels were shuffled globally because within-donor shuffling would leave the donor-by-class table unchanged. For transcriptomic class and ROI, E labels were shuffled independently within donor, preserving donor composition while disrupting the association of interest. Empirical one-sided P values used a plus-one correction.

### Software and reproducibility

Analyses were performed in Python 3.14.3 using NumPy 2.4.6, pandas 2.3.3, SciPy 1.17.1, scikit-learn 1.8.0, joblib 1.5.3, Matplotlib 3.10.9 and seaborn 0.13.2 on Windows 11. The primary validation seed was 20260910; submission-completion analyses used seed 20260912. Numerical source data, fold-level results, cell-level probabilities, transformation parameters and analysis scripts are retained with the project outputs.

## Reporting and interpretation boundaries

1. The supervised tests evaluate recovery of frozen consensus labels in unseen donors; they do not independently rediscover biological classes.
2. ROC AUC and AP are based on repeated internal donor-grouped validation, not an external cohort.
3. The K=4 solution is a frozen consensus reference, not the numerical maximum of every parameter scan.
4. Metadata associations are observational and should not be described as causal.
5. E1–E4 should be called reproducible electrophysiological classes or a frozen operational classification, not universally discrete natural cell types.

## Data-availability wording

Source data underlying all panels include fold-level performance, confusion matrices, donor-resampled clustering results, cell-level co-clustering margins, held-out probabilities, feature-importance estimates, nested-ablation results and metadata contingency tables. The frozen 390-cell transformed matrix, PC scores, graph-cluster mapping and 368-cell consensus assignments should accompany the source-data deposit. Analysis scripts and software versions should be released with a fixed version identifier or archive checksum.
