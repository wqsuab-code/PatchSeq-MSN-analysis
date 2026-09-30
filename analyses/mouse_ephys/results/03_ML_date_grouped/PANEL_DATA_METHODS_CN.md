# Mouse E1–E5 机器学习稳健性分析：逐图数据与方法要点

## 共同分析框架

- 分析对象：450 个 HC–GC 共识细胞；E1 = 142、E2 = 160、E3 = 45、E4 = 61、E5 = 42。
- 输入：冻结的 18 个电生理特征；不重新聚类、不重新定义 E 类别。
- 泄漏控制：按实验记录日期分组，同一天的所有细胞只能整体进入训练集或测试集。
- 主要监督分析：200 次 group-aware 75:25 留出拆分；每个拆分的训练集和测试集均包含 E1–E5。
- 数据变换：按训练折特征最小值平移、按训练折列和归一化至10,000、log1p和Z-score；需要PCA时，3-PC PCA也仅在训练集拟合。所有拟合参数原样应用于对应测试集。
- 随机种子：20260916。
- sequencing batch（Batch2–Batch8）仅作为元数据分析，不作为主要训练/测试分组变量。

## a | Held-out balanced accuracy

- 数据文件：`01_date_grouped_model_metrics.csv`；汇总文件：`02_model_summary.csv`。
- 数据单位：每个点是一次日期分组留出拆分的 held-out balanced accuracy；200 点/模型，7 个模型，共 1,400 行。
- 模型：Extra Trees、distance-weighted kNN、RBF SVM、linear SVM、class-weighted logistic regression、random forest、histogram gradient boosting。
- 箱线图概括 200 次拆分的分布；散点是原始拆分结果，不是单细胞。
- RBF SVM 的中位数最高：0.8705；2.5–97.5 百分位数为 0.7441–0.9699。

## b | Held-out macro-F1

- 数据文件：`01_date_grouped_model_metrics.csv`；汇总文件：`02_model_summary.csv`。
- 数据单位与 a 相同；每个点是一项日期分组测试拆分的 macro-F1。
- macro-F1 为五类 F1 的非加权平均，避免多数类别完全支配结果。
- RBF SVM：中位数 0.8607；2.5–97.5 百分位数为 0.7129–0.9651。

## c | Confusion (%)

- 原始计数：`04_confusion_counts.csv`；真实类别内百分比：`05_confusion_row_fraction.csv`。
- 行：冻结的 HC–GC 共识类别；列：RBF SVM 预测类别。
- 每个细胞在多次作为测试细胞时获得预测概率；先对该细胞的 held-out 概率求平均，再取最高概率类别。因此每个细胞只在最终矩阵中出现一次。
- 总体正确：411/450（91.33%）。该数值是重复 held-out 概率聚合后的细胞级正确率，不等于 a 图的拆分级 balanced accuracy。

## d | Permutation control

- 数据文件：`07_permutation_null.csv`。
- 方法：保持实验日区块结构，在相同细胞数的日期区块间置换标签；100 次置换，每次在固定的 5 个日期分组拆分中训练并评估 Extra Trees。
- 灰色分布：置换标签后的平均 balanced accuracy；观测竖线：未置换标签的 5 拆分平均值。
- 观测值 0.8778；置换均值 0.1984；置换范围 0.1689–0.2453；经验 P = (1 + 极端置换数)/(100 + 1) = 0.0099。
- 注意：该图检验“结果是否高于结构化随机标签”，不是模型间性能比较。

## e | One-vs-rest ROC

- 数据文件：`03_best_model_cell_predictions.csv` 和 `06_roc_pr_auc.csv`。
- 使用聚合的 RBF SVM held-out 概率，对每个 E 类进行 one-vs-rest ROC 分析。
- ROC-AUC：E1 0.9941、E2 0.9878、E3 0.9971、E4 0.9866、E5 0.9996。
- ROC 在类别不均衡时可能偏乐观，应与 f 图 PR 曲线共同解释。

## f | One-vs-rest precision–recall

- 数据文件：`03_best_model_cell_predictions.csv` 和 `06_roc_pr_auc.csv`。
- 横轴 recall，纵轴 precision；使用同一组聚合 held-out 概率。
- Average precision：E1 0.9858、E2 0.9816、E3 0.9764、E4 0.9197、E5 0.9966。
- PR 曲线对小类别与假阳性更敏感，是对 ROC 的必要补充。

## g | NPC and K sensitivity

- 数据文件：`08_NPC_K_sensitivity.csv`，63 个 NPC × K 组合。
- 全部 450 个细胞按冻结Mouse-E流程进行最小值平移、列和归一化至10,000、log1p和Z-score，再计算最多10个PC；对 NPC = 2–10、K = 2–8 扫描 Ward 层次聚类。
- 每个格子为该无监督划分与冻结 E1–E5 标签之间的 adjusted Rand index（ARI）。
- 最高 ARI = 0.6250（NPC = 3，K = 5），与预设的3-PC、5类结构一致。
- 这是参数敏感性分析，不是独立 held-out 预测结果。

## h | Experimental-day-resampled discovery

- 数据文件：`09_date_resampled_discovery.csv`，200 次日期分组重采样 × 4 种无监督算法 = 800 行。
- 每次仅在对应训练细胞中拟合最小值平移、列和归一化、log1p、Z-score、3-PC PCA与K = 5聚类。
- 算法：Ward、K-means、GMM、spectral clustering；评价为聚类结果相对冻结标签的 ARI。
- 该图检验在改变纳入实验日期后，无监督结构能否重复出现。

## i | Full-data algorithm agreement

- 数据文件：`10_full_data_algorithm_agreement.csv`。
- 全部 450 个细胞采用相同的 3-PC 表示，分别运行 K = 5 的 Ward、K-means、GMM 和 spectral clustering。
- 每个算法只有一个全数据 ARI，因此每根柱上仅有一个真实点，不能解释为重复实验分布。
- ARI：Ward 0.6250、K-means 0.2903、GMM 0.2831、spectral 0.3062。

## j | Cell-level consensus margin

- 数据文件：`11_cell_consensus_margin.csv`，每个细胞一行，共 450 行。
- 基于 h 图 200 次日期重采样中的 Ward 聚类，计算每对细胞在共同被抽到时的共聚类概率。
- 每个细胞的 margin = 同一冻结类别内平均共聚类概率 − 与其他类别中最高的平均共聚类概率。
- margin > 0 表示该细胞更稳定地与本类聚集；margin < 0 表示其更常与某一其他类聚集。
- 总体中位数 0.3423；56/450 个细胞 margin < 0；89/450 个细胞 margin < 0.1。

## k | Assignment stability by class

- 数据文件：`11_cell_consensus_margin.csv`。
- 将 j 图的细胞级 margin 按冻结 E1–E5 分组显示；每个点是一个细胞。
- 各类中位数：E1 0.45、E2 0.28、E3 0.37、E4 −0.03、E5 0.51。
- 该图用于识别哪一类的无监督重发现稳定性较弱；不是分类准确率。

## l | Association with metadata

- 数据文件：`12_metadata_association.csv`。
- 效应量：bias-corrected Cramér's V。
- 显著性：在保持实验日区块结构的条件下置换 E 标签 1,000 次，经验 P = (1 + 极端置换数)/1,001。
- T class：N = 438，V = 0.3414，P = 0.0010；sequencing batch：N = 450，V = 0.1089，P = 0.0100。
- V 表示关联强度，不表示方向或因果关系。

## m | Held-out permutation importance

- 数据文件：`13_permutation_importance.csv`，50 个日期分组拆分 × 18 个特征 = 900 行。
- 在每个拆分中，只用训练集拟合最小值平移、列和归一化、log1p、Z-score及class-weighted Extra Trees；随后在测试集逐个打乱某一特征8次。
- 指标为 held-out balanced accuracy 的平均下降量：`原始测试分数 − 特征打乱后的测试分数`。
- 正值越大，说明模型在未见日期的测试细胞上越依赖该特征；接近 0 表示影响弱；负值表示打乱后偶然改善，不能解释为“负生物学贡献”。
- 图中按 50 个拆分的平均值排序，并显示前 12 项。最高的两个特征为 input resistance 与 AP amplitude（平均下降均约 0.10），随后为 rheobase（约 0.06）。
- 这是模型预测依赖性，而不是单变量组间效应，也不证明因果作用。

## n | MDI stability

- 数据文件：`14_MDI_stability.csv`，50 拆分 × 18 特征 = 900 行。
- 每个拆分仅在训练集拟合 250-tree Extra Trees；MDI 为树中由该特征产生的加权不纯度下降。
- 图中使用各特征跨 50 拆分的 MDI 中位数并显示前 12 项。
- MDI 可能偏好可产生较多切分的连续特征，并会在相关特征之间分摊重要性；需与 m 图配合解释。

## o | Cross-model rank agreement

- 数据文件：`16_cross_model_ranks.csv`。
- 每个拆分分别计算 Extra Trees MDI 排名和 class-weighted multinomial logistic regression 绝对系数排名；再对 50 个拆分取平均排名。
- 点代表一个电生理特征；虚线表示两种模型排名完全一致。
- 两模型平均排名的 Spearman 相关为 ρ = 0.767，表明总体排序较一致，但并非所有特征一致。

## p | Progressive feature ablation

- 数据文件：`15_progressive_ablation.csv`，50 个日期分组拆分 × 4 个特征子集 = 200 行。
- 在每个拆分中依据训练集 Extra Trees MDI 排名定义子集，再重新训练 Extra Trees 并在对应未见日期的测试集计算 balanced accuracy。
- 子集：全部 18 项、去除前 2 项、去除前 5 项、仅保留排名最低 3 项。
- 中位 balanced accuracy：全部 18 项 0.86、去除前 2 项 0.80、去除前 5 项 0.75、仅保留最低 3 项 0.34。
- 该图验证高排名特征对分类性能的联合贡献；不能把性能下降完全归因于某一个特征。

## 报告时应避免混用的三个数值

1. a 图 RBF SVM 拆分级 balanced accuracy 中位数：0.8705。
2. b 图 RBF SVM 拆分级 macro-F1 中位数：0.8607。
3. c 图重复 held-out 概率聚合后的细胞级正确率：411/450 = 91.33%。

三者的统计单位和计算方式不同，不应统一写成一个“交叉验证准确率”。
