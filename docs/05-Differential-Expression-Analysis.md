# 📖 Bulk RNA-seq 数据分析最佳实践 (Part 5)

## 第五章：差异表达分析 —— 发现真正的生物学差异

差异分析的目标不仅仅是计算倍数变化 (Fold Change)，更重要的是评估这种变化是否具有**统计学显著性**。

### 5.1 统计学底座：负二项分布 (Negative Binomial Distribution)

为什么我们不能用简单的 T 检验？

> 📘 **核心概念：过离散 (Overdispersion)**
> 1. **泊松分布 (Poisson)**：在测序深度极高时，Reads 的计数理论上符合泊松分布（均值 = 方差）。
> 2. **现实情况**：由于个体差异，**方差往往远大于均值**。这种现象称为“过离散”。
> 3. **负二项分布**：通过引入“离散参数 (Dispersion)”，能完美建模这种高噪声数据。这是 DESeq2 和 edgeR 的核心数学基础。

---

### 5.2 离散度估算：信息共享的智慧 (Dispersion Estimation)

如果你只有 3 个重复，你几乎无法准确估计一个基因的真实方差。

#### 💡 Best Practice：离散度收缩 (Dispersion Shrinkage)
DESeq2 假设：具有相似表达水平的基因，其离散度应该是相似的。
1. 它先为每个基因计算一个单独的离散度。
2. 然后将这些点拟合到一条曲线上（全基因组平均水平）。
3. 最后，将单个基因的离散度向曲线方向“收缩 (Shrink)”。
*   **结果**：这种“信息共享”的方式极大地增强了统计效能，防止了异常样本导致的伪阳性。

---

### 5.3 LFC 收缩：让低表达基因保持安静 (LFC Shrinkage)

> ⚠️ **避坑指南 (Pitfall)：低表达基因的“虚假繁荣”**
> 在 MA 图中，低表达区域的基因往往具有极大的 Log2 Fold Change，但这通常是由于 Counts 很小时的随机波动导致的噪声。

#### 💡 Best Practice：使用 apeglm 算法进行收缩
DESeq2 会对低表达或高变异基因的 LFC 进行收缩。
*   **高表达/低噪声基因**：LFC 保持不变。
*   **低表达/高噪声基因**：LFC 会被拉回到 0。
*   **意义**：这让你的差异基因列表更加可靠，特别是在进行后续的 GSEA 富集分析时，不会被噪声干扰。

---

### 5.4 多重假设检验校正：FDR 与 p-adj

如果你测试了 20,000 个基因，按 p < 0.05 的标准，你会得到 1,000 个伪阳性。

> 📘 **核心概念：Benjamini-Hochberg (BH) 校正**
> 我们不看原始的 p-value，而看 **Adjusted p-value (p-adj)**，即 **FDR (错误发现率)**。
> **金标准：** 通常设定 **p-adj < 0.05** 且 **|log2FoldChange| > 1** 作为差异基因的标准。

---

### 🛠️ 实战演练：DESeq2 标准下游分析流程

```R
# 1. 运行差异分析主函数
dds <- DESeq(dds)

# 2. 提取结果并进行 LFC 收缩
res <- results(dds, contrast=c("condition", "Treated", "Control"))
res_shrunk <- lfcShrink(dds, coef="condition_Treated_vs_Control", type="apeglm")

# 3. 排序并保存
res_ordered <- res_shrunk[order(res_shrunk$padj),]
write.csv(as.data.frame(res_ordered), file="DEG_results.csv")
```
