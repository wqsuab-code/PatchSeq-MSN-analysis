# 冻结 Macaque MSN 形态学 M1-M4：结论与解释边界

## 冻结输入

仅使用 Macaque MSN、Ca + Pu + NAC、通过当前完整形态病例标准的细胞。原始 ZIP 重新筛选得到 126 个完整病例；HC-GC 共识细胞为 117 个，细胞级计数为 M1=43、M2=42、M3=20、M4=12，来自 41 个供体。18 项冻结树突/胞体指标均完整，因此未进行插补。E class 仅作为外部注释；117 个细胞中 91 个具有冻结 E4 共识注释。Mouse 数据、E 特征和 E 预处理参数均未进入分析。

## 六个预设问题

1. **冻结分类能否精确重现？** 能。由原始 ZIP 重算18项变换、Z-score、PCA5和 Ward.D2 K=4 后，HC相对冻结标签的 ARI=1.000、NMI=1.000、映射准确率=1.000，无不一致细胞。重算 Seurat SNN/Louvain（k=20、prune=1/15、resolution=2.3、seed=777）得到 raw K=13，原始GC标签不变性 ARI=1.000；合并GC准确率=1.000，无不一致细胞。

2. **更换供体后能否预测M标签？** 在同一形态特征域内可以较可靠恢复，但区间较宽。500个75/25供体隔离划分中，RBF SVM的 balanced accuracy 中位数为0.909（95%经验区间0.730-1.000），macro-F1为0.899（0.713-1.000）。M4 recall中位数为1.000，但2.5%分位仅0.286，反映M4样本少且测试供体组成敏感。供体内标签置换的balanced accuracy中位数为0.376，经验P=0.001996。

3. **从头重聚类时四类能否稳定再次出现？** 不能称为稳定再发现。500次保留约80%供体、每次重估变换和PCA后，冻结位置 NPC5/K4 的 Ward ARI中位数仅0.322（95%区间0.107-0.626），NMI中位数0.437，silhouette中位数0.155。最小群大小中位数为10，未出现小于等于2个细胞群，但“无极小群”不等于四类结构稳定。K-means、GMM和spectral也未显示接近完全恢复的跨供体发现稳定性。

4. **哪些细胞和类别最不稳定？** 共有16个细胞 consensus margin <=0，27个低于预设阈值0.10。M1和M2的margin分布最低且包含最多负值；M4的margin中位数最高（0.538），但仅有12个细胞，不能据此宣称其供体外稳定性最强。完整风险名单包含 donor、ROI、T class、E class与可用QC字段。

5. **哪些特征贡献较稳定？** Extra Trees MDI最靠前的是 basal_dendrite_stem_exit_MedialLateral、basal_dendrite_total_length、basal_dendrite_extent_medial、basal_dendrite_num_branches 和 basal_dendrite_mean_diameter。Extra Trees与多项logistic排序的Spearman rho=0.674，说明中等偏强但并非完全一致。训练折内嵌套 permutation importance 的区间普遍跨0，提示单个特征的独立贡献不稳，分类依赖多特征联合结构。Extra Trees消融的balanced accuracy中位数从全部18项的0.770降至删除前5项的0.712，仅保留最不重要3项时为0.369。

6. **分类是否可能由donor、ROI、批次或QC驱动？** donor的Cramer's V=0.612，但供体结构置换P=0.188，且41个供体造成列联表稀疏，不能据V值单独判定混杂。ROI关联为V=0.288、P=0.006，值得在后续生物学解释中明确控制或分层。T class为V=0.136、P=0.513；已知E class的91个细胞中V=0.257、P=0.084。源ZIP没有独立reconstruction batch或recording date字段，recording batch proxy也无变异；完整病例的18项缺失率均为0，因此这些因素无法由当前数据检验，不能解释为“无混杂”。

## 结论边界

冻结 M1-M4 的计算流程可精确复现，标签在同一18项形态特征空间中具有较高的供体外监督可恢复性。然而，供体重抽样后从头聚类与冻结四类的一致性较弱，且部分细胞margin低或为负。因此当前证据支持“可精确复现、监督可恢复的冻结操作型形态分类”，不支持仅凭本分析宣称四个天然离散生物学类型。独立生物学验证仍需来自未用于定义M的E class、T class、转录组、解剖位置或连接特征。

# Frozen Macaque MSN morphology M1-M4: conclusions and interpretation limits

The frozen workflow was reproduced exactly from the source ZIP (HC ARI/NMI/mapped accuracy=1.000; raw and merged GC reproduction=1.000). Donor-held-out supervised recoverability was high for the best model (RBF SVM median balanced accuracy 0.909, 95% empirical interval 0.730-1.000; macro-F1 0.899, 0.713-1.000), and exceeded the within-donor permutation distribution (empirical P=0.001996). In contrast, de novo donor-resampled discovery was modest at frozen NPC5/K4 (Ward median ARI 0.322, 0.107-0.626; median NMI 0.437). Sixteen cells had consensus margin <=0 and 27 were below 0.10. ROI showed a detectable association with M class, whereas T class and available E class did not reach the prespecified permutation threshold. Batch/date and reconstruction-completeness confounding could not be tested because the source data lacked informative fields.

These findings support a precisely reproducible and supervised-recoverable operational morphology classification. They do not establish four naturally discrete biological types. Orthogonal validation must use measurements not used to define M1-M4.
