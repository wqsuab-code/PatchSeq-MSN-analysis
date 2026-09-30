Macaque E4 投稿补充分析（2026-09-12）
=====================================

范围：Macaque MSN，Ca+Pu+NAC；冻结HC-GC共识细胞n=368，19项E特征，55个donor。

本目录补齐原机器学习证据链中四项投稿风险：

1. 匹配的Extra Trees置换检验
   - 真实标签与null均使用Extra Trees 200树、训练折内E变换和3-PC PCA。
   - 5折donor-grouped外层验证；每折训练和测试均含C1-C4。
   - 在每个donor内部置换E标签，保留donor的类别构成。
   - 500次置换；真实平均balanced accuracy=0.98174，null中位数约0.304，经验P=0.001996。

2. 严格训练折内permutation importance
   - 5次重复×5折StratifiedGroupKFold，共25个donor-grouped测试折。
   - 每折均从原始19项数据重新拟合训练集变换和Z-score。
   - Random forest 500树；每特征置换20次；scoring=balanced accuracy。
   - 前五项：AHP delay 5-spike、upstroke adaptation ratio、width adaptation ratio、AHP delay ratio 5-spike、upstroke/downstroke ratio at rheobase。

3. 嵌套feature ablation
   - 外层5×5 donor-grouped folds；每个外层训练集内以3折grouped CV重新估计特征排名。
   - 外层测试donor不参与特征选择。
   - Balanced accuracy中位数（95%经验区间）：
       全19项 0.891（0.788–0.954）
       删除top 2 0.895（0.756–0.935）
       删除top 5 0.776（0.710–0.863）
       bottom 3 0.535（0.436–0.638）

4. Metadata association修正
   - 采用bias-corrected Cramer's V和5000次结构化置换。
   - donor：V=0.2468，P=0.00020（全局E标签置换）。
   - T class：V=0.3286，P=0.00020（donor内E标签置换）。
   - ROI：V=0.1636，P=0.00160（donor内E标签置换）。

解释边界：这些结果支持冻结E4标签在未见donor中的可恢复性、特征扰动稳健性和超出donor结构的统计关联；它们仍不是独立生物学验证。Ward从头恢复同一四分类的donor-resampled稳定性仅为中等水平。

文件：
- 00_submission_completion_summary.json：关键结果
- tables/01-02：匹配置换检验
- tables/03-04：严格permutation importance
- tables/05-07：嵌套ablation
- tables/08-09：metadata检验及列联表
- 10_software_and_parameters.json：软件环境和参数
- Supplementary_ML4_submission_corrections.png/.pdf：四项补充结果图

运行脚本：scripts/complete_macaque_E4_submission_ml_audit.py
