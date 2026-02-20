# Bulk RNA-seq 数据分析最佳实践 (Part 5)

## 第五章：差异表达分析 (Differential Expression Analysis)

> 核心原则：差异分析不是“找大 Fold Change”，而是在噪声中做可靠推断。

---

## 5.1 统计模型基础：负二项分布与离散度

RNA-seq 计数通常存在过离散（方差 > 均值），DESeq2 用负二项模型处理这种特征。

### 为什么不是 t-test

- t-test 假设近似正态且方差结构简单。
- 计数数据离散且异方差明显，直接 t-test 容易误判。

DESeq2 的关键优势是：

1. 估计 size factor（样本深度校正）
2. 估计 dispersion（基因层噪声建模）
3. 进行稳健的广义线性模型检验

---

## 5.2 设计公式先于检验

示例：

```r
colData$condition <- relevel(factor(colData$condition), ref = "Control")
colData$batch <- factor(colData$batch)

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = colData,
  design = ~ batch + condition
)
```

> Best Practice：先把 `reference level` 定好，再跑 `DESeq()`，避免后面反复重算。

---

## 5.3 标准 DESeq2 流程

```r
library(DESeq2)

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast = c("condition", "Treated", "Control"))
res <- lfcShrink(dds, contrast = c("condition", "Treated", "Control"), type = "apeglm")

res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj), ]
write.csv(res_df, file = "results/DEG_results_shrunk.csv")
```

解释：

- `padj` 是 BH 校正后的 FDR。
- `lfcShrink` 让低计数基因的 fold change 更稳定，更适合排序和展示。

---

## 5.4 阈值设置与解释

常用阈值（建议作为初始标准）：

- `padj < 0.05`
- `|log2FoldChange| > 1`

```r
deg <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg, file = "results/DEG_sig.csv", row.names = TRUE)
```

> 避坑指南：阈值不应机械固定。对低效应但高一致性的通路，GSEA 往往更合适。

---

## 5.5 可视化：MA 图与火山图

```r
# MA plot
plotMA(res, ylim = c(-5, 5))
```

```r
library(EnhancedVolcano)

EnhancedVolcano(
  res_df,
  lab = rownames(res_df),
  x = 'log2FoldChange',
  y = 'padj',
  pCutoff = 0.05,
  FCcutoff = 1
)
```

图形解读应回答两件事：

1. 差异信号是否集中在合理范围。
2. 是否存在由少数异常样本驱动的“假阳性高 FC”。

---

## 5.6 高级设置：effect-size 检验

当你更关心“是否达到生物学意义的效应量”时，可用 `lfcThreshold`：

```r
res_lfc <- results(
  dds,
  contrast = c("condition", "Treated", "Control"),
  lfcThreshold = 1,
  altHypothesis = "greaterAbs"
)
```

这比事后硬筛 `|log2FC|` 更统计一致。

---

## 5.7 常见翻车点

1. 技术重复当生物学重复，显著性虚高。
2. 忘记设置参考组，导致方向解释反了。
3. 只看 p-value 不看效应量，得到“统计显著但生物意义弱”的结果。
4. 直接使用未收缩 LFC 排名做下游叙事。

---

## 5.8 本章检查清单

- 设计公式与实验问题一致。
- 对比方向（Treated vs Control）已核实。
- `padj` 与 `log2FC` 阈值已记录。
- DEG 结果表与图形输出已保存。

下一章进入功能富集与生物学解释。
