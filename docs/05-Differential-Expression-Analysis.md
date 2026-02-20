# Bulk RNA-seq 数据分析最佳实践 (Part 5)

## 第五章：差异表达分析 (Differential Expression Analysis)

本章按固定格式组织：先讲原理，再给代码，再给代码输出结果。

---

## 5.1 原理：DEG 的核心是“带离散度建模的统计检验”

- RNA-seq 是计数型数据，通常存在过离散
- DESeq2 用负二项模型估计 size factor、dispersion，再做检验
- 结果解释要同时看显著性（`padj`）和效应量（`log2FoldChange`）

---

## 5.2 代码：运行 DESeq2 并导出结果

```r
library(DESeq2)

# 假设 dds 已由 tximport + colData 构建
dds <- DESeq(dds)
resultsNames(dds)

# WT: AvrRpm1 vs mock（LFC shrink）
res_wt <- lfcShrink(dds, coef = "condition_AvrRpm1_vs_mock", type = "apeglm")
res_wt_df <- as.data.frame(res_wt)
res_wt_df <- res_wt_df[order(res_wt_df$padj), ]
write.csv(res_wt_df, "artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_shrunk.csv")

# 阈值筛选
deg_wt <- subset(res_wt_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg_wt, "artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv")

# 交互项（clf vs WT 的处理响应差异）
res_inter <- results(dds, name = "genotypeclf.conditionAvrRpm1")
res_inter_df <- as.data.frame(res_inter)
res_inter_df <- res_inter_df[order(res_inter_df$padj), ]
write.csv(res_inter_df, "artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT.csv")

deg_inter <- subset(res_inter_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg_inter, "artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT_sig.csv")
```

---

## 5.3 代码：统计 DEG 数量并查看结果头部

```bash
# 显著 DEG 数量（WT）
tail -n +2 artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv | wc -l

# 交互项显著 DEG 数量
tail -n +2 artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT_sig.csv | wc -l

# 查看 WT 显著 DEG 前 5 行
head -n 6 artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv

# 查看交互项显著 DEG 前 5 行
head -n 6 artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT_sig.csv
```

### 输出结果

```text
1541
56
"","baseMean","log2FoldChange","lfcSE","pvalue","padj"
"AT5G15950",1058.46429885084,-1.45376652095534,0.132731770663756,5.16733272722286e-29,4.50074680541111e-25
"AT1G09935",24.1197681031312,3.74698979819632,0.364417164046629,2.44079011892465e-26,1.41728546238891e-22
"AT3G17420",88.8579253411927,1.63508829634183,0.162226020945527,4.51659669259174e-25,1.57358228769896e-21
"AT1G09310",2457.14852891811,-1.39998771718618,0.142523808393732,6.84344059187218e-24,1.7030390730059e-20
"AT1G64060",266.611658384728,1.18063546805895,0.122747998195324,7.40067456565989e-23,1.61149688667244e-19
"","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
"AT1G67070",49.4803711208944,-1.5755619595937,0.313552297126098,-5.02487774458903,5.03753339073799e-07,0.00168686843122252
"AT3G61060",127.599719325857,-1.22285834313158,0.253060089891281,-4.83228447305515,1.34975200454455e-06,0.00251098864578772
"AT3G55560",39.2719988293596,-1.19171445089373,0.253284246997521,-4.70504764911572,2.53806668243077e-06,0.00265592815399615
"AT4G38470",206.993199756586,-1.0452717288759,0.221874001988044,-4.71110503939182,2.46377201395389e-06,0.00265592815399615
"AT5G13400",140.757962415624,1.19134842011623,0.250851013492529,4.7492270552525,2.04195583734993e-06,0.00265592815399615
```

---

## 5.4 代码：结果图（MA + Volcano）

```bash
Rscript scripts/05_generate_case_figures.R
ls docs/assets/validated_case/ch5_*.png
```

### 输出结果

```text
docs/assets/validated_case/ch5_ma_wt.png
docs/assets/validated_case/ch5_volcano_wt.png
```

火山图：

![Volcano WT](assets/validated_case/ch5_volcano_wt.png)

MA 图：

![MA WT](assets/validated_case/ch5_ma_wt.png)

---

## 5.5 本章产物（网页可见）

- `artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_shrunk.csv`
- `artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv`
- `artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT.csv`
- `artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT_sig.csv`
- `docs/assets/validated_case/ch5_volcano_wt.png`
- `docs/assets/validated_case/ch5_ma_wt.png`
