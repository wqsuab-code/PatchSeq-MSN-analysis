# 当前绘图参数

权威参数文件为`layout/Macaque_E4_ML_layout.json`。本页是便于人工查阅的摘要；若与JSON冲突，以JSON为准。

## 全局参数

- 画布：1200×900 SVG units；4列；水平/垂直gap均20。
- 背景：`#FFFFFF`；字体：Arial。
- 默认标题10、panel字母11、坐标刻度6、线宽1、点大小1.5、点透明度0.7。
- E1 `#F8766D`；E2 `#7CAE00`；E3 `#00BFC4`；E4 `#C77CFF`。
- model `#149174`；bar `#4C78A8`；neutral `#BDBDBD`；accent `#E76F00`；heat low `#203B74`；heat high `#F6E928`。
- 每个panel外框尺寸均为275×205；位置按四列为x=20、315、610、905，四行为y=20、245、470、695。
- 默认plot margins：left=48、right=14、top=31、bottom=36；每个panel可在JSON中独立覆盖。

## 逐panel轴与标记参数摘要

|Panel|ID|Y范围|X范围|主色/次色|点大小|点透明度|分类标签角度|
|---|---|---|---|---|---:|---:|---:|
|a|ba|0.5–1|auto|#FFFFFF / #333333|0.3|1|默认|
|b|f1|0.5–1|auto|#FFFFFF / #333333|0.3|1|默认|
|c|confusion|auto|auto|#FF6B6B / #FFF5F5|1.5|0.7|默认|
|d|permutation|0–100|0.05–1|#1F2B82 / #FF0000|3|0.7|默认|
|e|roc|0–1|0–1|#008B74 / #999999|1.5|0.7|默认|
|f|pr|0–1|0–1|#008B74 / #999999|1.5|0.7|默认|
|g|sensitivity|auto|auto|#FFF4F2 / #D96A6A|1.5|0.7|默认|
|h|discovery|0–1|auto|#FFFFFF / #333333|0.3|1|默认|
|i|agreement|0–1|categorical|#5A5AA5 / #333333|1.5|0|0°|
|j|margin|0–40|−1–1|#5A5AA5 / #EAE6E1|1.5|0.7|默认|
|k|classmargin|−1–1|categorical|#FFFFFF / #333333|1|0.9|默认|
|l|metadata|auto|auto|#D9D9DF / #333333|1.5|0.7|默认|
|m|importance|auto|0–auto|#FF5B55 / #333333|1.5|0.7|默认|
|n|mdi|auto|0–auto|#FF6F6F / #333333|1.5|0.7|默认|
|o|rank|auto|auto|#D7DCF7 / #999999|3|1|默认|
|p|ablation|0.15–1|categorical|#FFFFFF / #333333|1|0.9|0°|

## 箱线图规范

- 空白填充箱体；中线为median；须为1.5×IQR。
- 透明灰色点显示各次重抽样/外层fold。
- 实心点表示mean。
- 不叠加mean 95% CI误差条，避免与1.5×IQR须重复并避免对相关重抽样给出误导性独立样本CI。
