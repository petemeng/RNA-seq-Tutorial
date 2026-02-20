# 📘 Bulk RNA-seq 进阶专题 (Part 8)

## 第八章：时间序列分析 —— 捕捉动态转录组

时间序列数据的关键价值是“轨迹信息”。
如果你只做 `T1 vs T0`、`T2 vs T0` 这种两两比较，会丢失动态过程中的整体模式。

---

## 8.1 什么时候必须用时间序列模型？

满足以下任意条件，就不建议只做 pairwise：

1. 时间点 `>= 3`。
2. 你关心“先升后降”或“短暂脉冲”这类非单调变化。
3. 你需要比较不同处理在时间维度上的响应差异（interaction）。

> 💡 **核心概念：LRT（似然比检验）**
> 用 Full model 和 Reduced model 比较，直接回答：
> “这个基因的表达是否随时间（或时间交互）显著变化？”

---

## 8.2 设计公式怎么选？

| 生物学问题 | DESeq2 设计 | LRT reduced |
| :--- | :--- | :--- |
| 单条件时间变化 | `~ time` | `~ 1` |
| 两条件共同时间效应 | `~ treatment + time` | `~ treatment` |
| 关注条件-时间交互 | `~ treatment + time + treatment:time` | `~ treatment + time` |

> ⚠️ 避坑指南：`time` 默认建议当作 **factor**（分类变量），除非你明确假设线性趋势。

---

## 8.3 实战演练：DESeq2 做时间序列 LRT

```r
library(DESeq2)

# PRJDB11848 示例：WT 子集，time 为 0h/0.5h/3h，condition 为 mock/AvrRpm1
colData$time <- factor(colData$time, levels = c("0h", "0.5h", "3h"))
colData$condition <- relevel(factor(colData$condition), ref = "mock")

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = colData,
  design = ~ condition + time + condition:time
)

# 检验“是否存在时间相关响应（含交互）”
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ condition + time)
res_lrt <- results(dds_lrt)

sig_dynamic <- subset(as.data.frame(res_lrt), padj < 0.05)
nrow(sig_dynamic)
```

解释：

- `padj < 0.05`：该基因存在显著时间交互动态。
- 这一步是“筛动态基因”；具体在哪个时间点变，需要下一步 Wald/contrast。

---

## 8.4 提取具体时间点效应（Wald + contrast）

```r
# 同一设计下再跑 Wald（默认）
dds_wald <- DESeq(dds)
resultsNames(dds_wald)

# 在参考时间点（0h）下，处理主效应
res_treat_0h <- results(dds_wald, name = "condition_AvrRpm1_vs_mock")

# 处理在 3h 的总效应 = 主效应 + interaction(3h)
res_treat_3h <- results(
  dds_wald,
  list(c("condition_AvrRpm1_vs_mock", "conditionAvrRpm1.time3h"))
)
```

> ⚠️ 避坑指南：交互模型中，“主效应”只对应参考水平。
> 一定先 `resultsNames(dds)` 再提取，避免拿错系数。

---

## 8.5 动态模式可视化（聚类热图）

```r
library(pheatmap)

vsd <- vst(dds_wald, blind = FALSE)
sel_genes <- rownames(sig_dynamic)

mat <- assay(vsd)[sel_genes, ]
mat_z <- t(scale(t(mat)))

anno <- as.data.frame(colData(dds_wald)[, c("time", "treatment")])
pheatmap(
  mat_z,
  show_rownames = FALSE,
  annotation_col = anno,
  cluster_cols = FALSE
)
```

教学建议：

- 先按时间顺序排列样本，再做聚类，更容易解释“早期响应/晚期响应”。
- 可对基因做 k-means 分组，得到“瞬时诱导型”“持续抑制型”等动态模块。

---

## 8.6 进阶工具：ImpulseDE2 何时值得上？

当你预期是“脉冲型/非线性”表达轨迹（例如 1h 急升、3h 回落），`ImpulseDE2` 往往优于线性模型。

```r
# 伪代码示意：具体输入格式请按 ImpulseDE2 文档整理
# library(ImpulseDE2)
# res_imp <- runImpulseDE2(
#   matCountData = counts(dds),
#   dfAnnotation = annotation_df,
#   boolCaseCtrl = TRUE,
#   scaNProc = 8
# )
```

---

## 8.7 常见翻车点清单

1. 每个时间点只有 1-2 个重复，导致动态结论不稳。
2. 没有 `0h` 基线，后续解释非常困难。
3. 把时间当 numeric 强行拟合线性，错过非线性轨迹。
4. 把受试者重复测量当作独立样本（需要 mixed model 思路）。

---

## 8.8 本章小结

时间序列分析的关键是先用 LRT 找“是否动态”，再用 contrast 定位“何时变化”。
如果你的实验问题是“响应过程”而不是“终点差异”，时间序列模型是必须项，不是可选项。

---

## 8.9 本教程实跑代码与结果（PRJDB11848）

### 代码：Chapter 8 实跑片段（对应 `scripts/04_downstream_ch4_to_ch9.R`）

```r
library(DESeq2)
library(tximport)

# WT 子集
wt_samples <- samples[samples$genotype == "WT", ]
wt_files <- file.path(base_dir, "data", "quant", wt_samples$sample_id, "quant.sf")
names(wt_files) <- wt_samples$sample_id

txi_wt <- tximport(wt_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

dds_wt <- DESeqDataSetFromTximport(
  txi_wt,
  colData = as.data.frame(wt_samples),
  design = ~ condition + time + condition:time
)

keep_wt <- rowSums(counts(dds_wt) >= 10) >= 3
dds_wt <- dds_wt[keep_wt, ]
dds_wt_lrt <- DESeq(dds_wt, test = "LRT", reduced = ~ condition + time)

res_lrt <- results(dds_wt_lrt)
res_lrt_df <- as.data.frame(res_lrt)
res_lrt_df <- res_lrt_df[order(res_lrt_df$padj), ]
write.csv(res_lrt_df, file.path(base_dir, "results/ch8/time_series_LRT_WT.csv"))

sig_dynamic <- subset(res_lrt_df, !is.na(padj) & padj < 0.05)
write.csv(sig_dynamic, file.path(base_dir, "results/ch8/time_series_LRT_WT_sig.csv"))
```

### 代码：padj 分布图（对应 `scripts/05_generate_case_figures.R`）

```r
library(ggplot2)
library(readr)

lrt <- read_csv(file.path(base_dir, "results/ch8/time_series_LRT_WT.csv"), show_col_types = FALSE)

p_lrt <- ggplot(lrt, aes(padj)) +
  geom_histogram(bins = 50, fill = "#33a02c", color = "white") +
  theme_bw(base_size = 12) +
  labs(title = "LRT adjusted p-value distribution (WT)", x = "padj", y = "Gene count")

grDevices::png("docs/assets/validated_case/ch8_lrt_padj_hist.png", width = 6.5, height = 4.8, units = "in", res = 150)
print(p_lrt)
grDevices::dev.off()
```

### 代码：验收命令

```bash
tail -n +2 artifacts/prjdb11848/results/ch8/time_series_LRT_WT_sig.csv | wc -l
head -n 3 artifacts/prjdb11848/results/ch8/time_series_LRT_WT_sig.csv
```

### 输出结果

```text
2927
"","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
"AT5G04140",4003.82180526792,-1.21881482066899,0.154193085409381,158.375766513059,4.06573686942051e-35,6.76945188758516e-31
"AT2G37710",4010.27684144751,1.10955082167625,0.179225289536766,144.702484110701,3.78664696027715e-32,3.15238359443073e-28
```

padj 分布图：

![LRT padj](assets/validated_case/ch8_lrt_padj_hist.png)
