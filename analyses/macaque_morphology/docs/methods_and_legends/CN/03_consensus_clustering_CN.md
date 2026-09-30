# HC–GC共识分类

HC在PC1–PC5空间使用欧氏距离和Ward.D2，并切为K=4。GC使用k=20、prune.SNN=1/15的SNN图及Louvain algorithm 1（resolution=2.3，seed=777），得到13个原始社区；每个社区按与HC的最大细胞重叠映射至4类。HC与合并GC标签一致的117/126个细胞构成冻结共识集（92.86%，ARI=0.825）：M1=43、M2=42、M3=20、M4=12。9个不一致细胞不进入共识类汇总和监督恢复性分析。T class、ROI和供体均未参与M标签构建。
