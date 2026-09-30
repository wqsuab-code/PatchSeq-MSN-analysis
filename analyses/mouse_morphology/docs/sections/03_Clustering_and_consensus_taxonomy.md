# Hierarchical and graph clustering

Hierarchical clustering was performed in the PC1-PC3 space using Euclidean distances and Ward.D2 linkage. The dendrogram was cut at K=4, yielding the HC-derived M1-M4 labels used as the complete 187-cell taxonomy (M1, n=67; M2, n=23; M3, n=62; M4, n=35).

Graph clustering was performed independently in the same PC1-PC3 space. A shared-nearest-neighbour graph was constructed with k=20 and pruning threshold 1/15. Louvain clustering (Seurat algorithm 1; 30 starts; 30 iterations; random seed 20260825) was run across resolutions 0.50-3.00. Resolution 0.50 produced four graph clusters and was selected for the frozen comparison. Graph labels were matched to HC labels by maximum overlap. HC and graph clustering agreed for 181 of 187 cells (96.79%). The term HC-GC consensus therefore refers only to these 181 cells; the 187-cell M taxonomy retains the HC-derived label for all included cells.
