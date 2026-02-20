# Bulk RNA-seq 数据分析最佳实践 (Part 4)

## 第四章：标准化与探索性分析 (Normalization & EDA)

> 核心原则：在统计检验前先理解数据结构，否则显著性没有意义。

---

## 4.1 为什么 DESeq2 不接受 TPM 做差异分析

### 关键点

1. DESeq2 建模对象是原始计数（负二项分布）。
2. 它内部会估计 `size factor` 来校正测序深度与组成偏差。
3. TPM/FPKM 可用于展示表达丰度，不用于 DESeq2 的假设检验输入。

> 避坑指南：把 TPM 喂给 DESeq2 是最常见错误之一，会破坏模型假设。

---

## 4.2 过滤低表达基因：减少噪声、提高统计效率

示例策略（按项目可调）：

- 保留在至少 `k` 个样本中 `count >= 10` 的基因。
- `k` 通常设为最小组重复数。

```r
library(DESeq2)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
```

> Best Practice：过滤阈值写进方法部分，保证可复现。

---

## 4.3 数据变换：VST 与 rlog

用于可视化和距离计算时，建议先做方差稳定化。

```r
vsd <- vst(dds, blind = FALSE)
# 小样本时也可考虑
# rld <- rlog(dds, blind = FALSE)
```

### 选择建议

- `vst`：速度快，适合多数场景。
- `rlog`：样本很少时可能更稳健，但更慢。

---

## 4.4 EDA 第一张图：PCA

```r
plotPCA(vsd, intgroup = c("condition", "batch"))
```

解释重点：

1. 生物重复是否聚在一起。
2. 组间分离是否符合实验预期。
3. 主成分是否被 batch 主导。

如果样本按批次分而不是按处理分，优先考虑在模型中加入 `batch`。

---

## 4.5 EDA 第二张图：样本距离热图

```r
library(pheatmap)

sample_dists <- dist(t(assay(vsd)))
sample_mat <- as.matrix(sample_dists)

pheatmap(
  sample_mat,
  annotation_col = as.data.frame(colData(dds)[, c("condition", "batch")])
)
```

热图可以快速发现离群样本和潜在标签错误。

---

## 4.6 批次效应处理策略

### 首选：在设计公式显式纳入

```r
design(dds) <- ~ batch + condition
```

### 仅用于可视化的批次去除（不要用于差异建模输入）

```r
library(limma)
mat <- assay(vsd)
mat_corrected <- removeBatchEffect(mat, batch = colData(dds)$batch)
```

> 避坑指南：不要把“批次校正后的表达矩阵”再喂回 DESeq2 重新做检验。

---

## 4.7 异常样本处理

建议流程：

1. 先通过 FastQC/MultiQC 看是否有技术原因。
2. 再结合 PCA/距离热图判断是否持续离群。
3. 如果剔除样本，必须在记录里写清楚依据。

> Best Practice：样本剔除标准提前定义，避免结果导向删样本。

---

## 4.8 常见翻车点

1. 直接在原始 count 上做 PCA。
2. 未过滤低表达基因，导致噪声主导聚类。
3. 批次效应明显却未纳入设计。
4. 看到离群点直接删除但没有可追溯理由。

---

## 4.9 本章检查清单

- 已完成低表达过滤并记录阈值。
- 已完成 VST/rlog 变换。
- 已检查 PCA 与样本距离热图。
- 已确定批次是否进入设计公式。

完成后再进入差异表达分析。

---

## 4.10 本教程实跑代码与结果（PRJDB11848）

### 代码：Chapter 4 实跑片段（对应 `scripts/04_downstream_ch4_to_ch9.R`）

```bash
Rscript - <<'RSCRIPT'
library(DESeq2)
library(pheatmap)

# dds 已在前面章节构建并完成 DESeq()
vsd <- vst(dds, blind = FALSE)

# PCA 数据导出
pca_df <- plotPCA(vsd, intgroup = c("genotype", "condition", "time"), returnData = TRUE)
write.csv(pca_df, "validation_run_downstream/results/ch4/PCA_data.csv", row.names = FALSE)

# 样本距离热图
sample_dists <- dist(t(assay(vsd)))
sample_mat <- as.matrix(sample_dists)
pdf("validation_run_downstream/results/ch4/sample_distance_heatmap.pdf", width = 8, height = 6)
pheatmap(sample_mat, annotation_col = as.data.frame(colData(dds)[, c("genotype", "condition", "time")]))
dev.off()
RSCRIPT
```

### 代码：验收命令

```bash
wc -l artifacts/prjdb11848/results/ch4/PCA_data.csv
ls docs/assets/validated_case/ch4_pca.png
```

### 输出结果

```text
37 /home/data/t060551/Codex/RNA-seq-Tutorial/artifacts/prjdb11848/results/ch4/PCA_data.csv
/home/data/t060551/Codex/RNA-seq-Tutorial/docs/assets/validated_case/ch4_pca.png
```

PCA 实跑结果图：

![PCA](assets/validated_case/ch4_pca.png)
