const fmt=(x,d=3)=>Number(x).toFixed(d);
async function json(path){const r=await fetch(path);if(!r.ok)throw Error(`${path}: ${r.status}`);return r.json()}
function card(value,label){return `<article><strong>${value}</strong><span>${label}</span></article>`}
async function main(){
 const [audit,manifest]=await Promise.all([json('ml-data/preprocessing_audit.json'),json('ml-data/manifest.json')]);
 const s=manifest.summary,c=manifest.corrected_summary;
 document.querySelector('#heroMetrics').innerHTML=[card('390 / 368','完整病例 / 冻结共识细胞'),card('55','供体'),card('19','E-features'),card('500×','供体感知重复')].join('');
 document.querySelector('#auditCards').innerHTML=[
  card(audit.max_abs_difference_rebuilt_vs_frozen_z.toExponential(2),'重建Z-score与冻结矩阵的最大绝对差'),
  card(audit.max_abs_difference_rebuilt_vs_frozen_pc_after_sign_alignment.toExponential(2),'符号对齐后PCA score最大绝对差'),
  card(fmt(audit.ward_k4_ARI_vs_frozen),'Ward K=4对冻结标签的ARI'),
  card(audit.yeo_johnson_used?'YES':'NO','主分析是否使用Yeo–Johnson')].join('');
 const meta=Object.fromEntries(c.metadata.map(x=>[x.metadata,x]));
 document.querySelector('#resultGrid').innerHTML=[
  card(fmt(s.best_supervised_BA_median),'Extra Trees供体外BA中位数'),
  card(fmt(c.observed_ET_grouped5fold_mean_BA),'严格分组五折平均BA'),
  card(fmt(c.matched_within_donor_null_median),'供体内置换零分布中位数'),
  card(`P=${c.matched_permutation_empirical_p.toPrecision(3)}`,'500次匹配置换检验'),
  card(fmt(s.Ward_resampling_ARI_median),'Ward供体重抽样ARI中位数'),
  card(fmt(s.Spectral_resampling_ARI_median),'Spectral供体重抽样ARI中位数'),
  card(fmt(meta.donor_label.cramers_v_bias_corrected),'Donor偏差校正Cramér’s V'),
  card(fmt(meta.T_class.cramers_v_bias_corrected),'T class偏差校正Cramér’s V'),
  card(fmt(meta.Lib_region_of_interest_label.cramers_v_bias_corrected),'ROI偏差校正Cramér’s V'),
  card('16 / 34','margin≤0 / margin<0.1的细胞数')].join('');
}
main().catch(e=>{document.body.insertAdjacentHTML('afterbegin',`<div style="padding:10px;background:#ffd9d9;color:#7a1111">Data load failed: ${e.message}</div>`)});
