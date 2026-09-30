# Macaque E1–E4机器学习验证：方法、结果与可重复性披露

## 分析范围与冻结标签

分析仅包括Macaque MSN及Ca、Pu和NAc三个ROI。冻结分类的输入为390个在19项非冗余E-feature上均无缺失的细胞；HC与合并GC一致的368个细胞进入监督机器学习分析（E1=59、E2=133、E3=57、E4=119），涉及55个供体。冻结标签由PC1–PC3上的Ward.D2层次聚类（K=4）与Seurat SNN/Louvain图聚类（`k.param=20`、`prune.SNN=1/15`、`resolution=3.0`、`algorithm=1`、`seed=777`）的合并四类共识得到。标签是在机器学习验证之前冻结的，机器学习没有重新定义E1–E4。

## 与冻结分类一致的数据变换

主分析没有使用Yeo–Johnson。冻结及机器学习主流程均采用特征特异的变换：

1. 名称含`ratio`的7项严格正值指标先进行`log2(x)`；
2. 其余12项指标依次执行训练集最小值平移、按训练集平移后列和归一化至10,000、`log1p`；
3. 所有指标使用训练集均值和总体标准差（`ddof=0`）进行Z-score；
4. 需要PCA的模块仅在训练集拟合full-SVD PCA，并使用PC1–PC3；
5. 测试集完全沿用训练集估计的最小值、列和、均值、标准差和PCA载荷。测试值低于训练最小值时截到训练下限。

在全部390个完整病例上重建该流程，与冻结矩阵的最大绝对差为`1.78×10⁻15`；符号对齐后的PCA score最大绝对差为`5.71×10⁻13`，解释方差比例最大绝对差为`1.67×10⁻16`。重建Ward K=4与冻结HC标签的ARI=1.000、NMI=1.000。因此，当前a–p主结果与冻结Macaque E1–E4预处理一致。任何Yeo–Johnson版本只能单独标为敏感性分析，不得替代主结果。

## 供体感知验证设计

监督分析以冻结共识细胞为对象，生成500个有效的供体分组训练/测试拆分；约75%的供体用于训练、25%用于测试，并要求训练与测试均包含四类。每次拆分均从原始19项指标开始重新拟合变换和PCA，供体不跨训练集和测试集。比较了Extra Trees、k-nearest neighbours、RBF SVM、linear SVM、multinomial logistic regression、random forest和histogram gradient boosting。主指标为balanced accuracy和macro-F1；随机基线为0.25。

Extra Trees的累计供体外预测用于混淆矩阵、类别召回率、one-vs-rest ROC和precision–recall分析。ROC/PR基于每个细胞在其所有供体外出现中的平均预测概率，保证368个细胞均至少被留出一次。

标签置换采用固定的严格供体分组五折结构，在每个供体内部置换E标签500次；每次从头拟合训练折变换、PCA和Extra Trees。经验P值按`(1 + null ≥ observed)/(1 + permutations)`计算。

## 无监督发现稳定性

每次重复随机抽取80%的供体，重新拟合变换和PCA，并在NPC=2–10、K=2–8下评估Ward聚类。另在PC1–PC3、K=4条件下比较Ward、K-means、full-covariance Gaussian mixture和nearest-neighbour spectral clustering。所得标签与冻结E4标签使用ARI比较。该部分衡量从扰动数据中重新发现四类的能力，不等同于监督标签可恢复性。

## 细胞稳定性、元数据与特征解释

细胞共聚类margin定义为同一冻结类内平均共聚类概率减去最高其他类平均共聚类概率；margin≤0表示在重抽样中其他类的共聚类支持不低于本类。元数据关联使用偏差校正Cramér's V；T class和ROI使用供体内标签置换，donor使用全局标签置换，各5,000次。

特征重要性包括两条互补证据链：随机森林的严格测试折permutation importance，以及500个供体留出拆分中的Extra Trees mean decrease in impurity。跨模型排序一致性使用Spearman相关。渐进消融在外层供体分组测试之前，仅用训练数据中的内层供体分组验证对特征排序；随后比较全部19项、去除前2项、去除前5项和仅保留末3项。

## 主要结果

- Extra Trees供体外balanced accuracy中位数为0.973，95%经验区间为0.904–1.000。
- 严格供体分组五折Extra Trees平均balanced accuracy为0.982；匹配的供体内置换零分布中位数为0.304，经验`P=0.001996`。
- E1–E4的one-vs-rest ROC-AUC为1.000、1.000、0.999和0.999；average precision为1.000、1.000、0.997和0.999。
- 80%供体重抽样下Ward的ARI中位数为0.476（95%经验区间0.253–0.776），spectral clustering中位数为0.719，提示从头发现四群对抽样和算法较敏感。
- 390个细胞的共聚类margin中，16个≤0，34个<0.1，提示少数边界细胞的无监督归属不稳定。
- 偏差校正Cramér's V：donor=0.247、T class=0.329、ROI=0.164；结构化置换P分别为0.00020、0.00020和0.00160。这些是关联/混杂审计，不是分类器性能。
- 训练折限定的前5项重要指标为AHP delay（5-spike）、upstroke adaptation ratio、AP-width adaptation ratio、AHP-delay ratio（5-spike）和rheobase upstroke/downstroke ratio。
- 嵌套消融balanced accuracy中位数：全部19项0.891、去除前2项0.895、去除前5项0.776、仅末3项0.535。前两项存在冗余，但前五项整体携带不可替代的信息。

## 解释边界

结果支持：冻结E1–E4标签可以从相同电生理特征域中在未见供体上高精度恢复，并显著高于匹配置换基线。结果同时显示：从供体扰动数据中无监督地重新发现完全相同四群只有中等稳定性且具有算法依赖性。因此，E1–E4应表述为“稳定、可操作且供体外可恢复的冻结电生理分类”，而不能仅凭这些结果宣称为算法无关、天然离散的生物类型。

## 软件与代码入口

- 主500次验证：`scripts/validate_macaque_E4_full_ml500.py`
- 基础变换和跨算法重抽样：`scripts/validate_macaque_EM_stability_500.py`
- ROC/PR：`scripts/add_macaque_E4_roc_pr.py`
- 严格置换、训练折特征重要性、嵌套消融和元数据置换：`scripts/complete_macaque_E4_submission_ml_audit.py`
- 冻结流程数值审计：`scripts/audit_macaque_E4_ml_preprocessing.py`
- 网页数据构建：`outputs/Macaque_E4_feature_explorer_site/build-ml-data.mjs`
- 交互页面：`outputs/Macaque_E4_feature_explorer_site/dist/ml-editor.html`

固定随机种子分别记录在各脚本与输出参数文件中。主要验证使用Python、NumPy、pandas、SciPy、scikit-learn、joblib、Matplotlib和seaborn；精确版本应以运行时生成的参数/版本记录为准。
