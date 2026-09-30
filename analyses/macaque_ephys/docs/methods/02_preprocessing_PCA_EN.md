# Feature transformation and principal-component analysis

All transformations were fitted to the 390 complete cases. The seven strictly positive ratio variables were log2-transformed. For each of the remaining 12 variables, the column minimum was subtracted, the shifted column was normalized to a total of 10,000, and `log1p` was applied. Every transformed variable was then standardized to zero mean and unit variance using the population standard deviation (`ddof=0`). No row normalization, imputation, batch correction, covariate regression, or Yeo-Johnson transformation was used in the frozen classification.

Principal-component analysis was performed on the 390 x 19 standardized matrix by full singular-value decomposition. The first three components were used for clustering and explained 23.2751%, 17.5400%, and 13.8060% of the variance, respectively (cumulative, 54.6211%). PC4 was examined only in sensitivity analyses.
