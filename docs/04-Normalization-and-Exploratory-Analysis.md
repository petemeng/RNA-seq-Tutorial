# 📖 Bulk RNA-seq 数据分析最佳实践 (Part 4)

## 第四章：标准化与探索性分析 —— 消除噪声，发现信号

在拿到基因表达矩阵（Counts Matrix）后，我们不能直接进行比较。因为不同样本之间的**测序深度**不同，不同基因之间的**长度**不同。

### 4.1 标准化之争：RPKM, FPKM vs. TPM (The Big Debate)

这是生信面试和论文评审中最常被提及的话题。

> 📘 **核心概念：为什么需要标准化？**
> 1. **测序深度 (Sequencing Depth)**：样本 A 测了 10M reads，样本 B 测了 20M。
> 2. **基因长度 (Gene Length)**：基因 X 长 10kb，基因 Y 长 1kb。在同等表达水平下，X 会产生更多的 Reads。
> 3. **RNA 组成 (RNA Composition)**：极高表达基因会“占领”测序空间，导致其他基因的相对丰度下降。

#### ⚠️ 避坑指南：为什么 RPKM/FPKM 已经被学界抛弃？
*   **RPKM/FPKM** 无法保证不同样本之间标准化后的总和相等，导致样本间不可比。
*   **TPM (Transcripts Per Million)** 确保了每个样本的总和一致，是跨样本比较表达丰度的更好选择。
*   **🚨 终极警告**：**千万不要将 TPM 或 FPKM 直接输入 DESeq2！** DESeq2 需要的是原始 Count 值，它有内置的 Median-of-Ratios 标准化算法。

---

### 4.2 数据变换：为了可视化 (rlog vs. VST)

原始 Count 数据具有“异方差性”：表达量越高的基因，其绝对方差也越大。这会干扰 PCA 和聚类分析。

> 🛠️ **实战演练：在 R 中进行转换**
> *   **VST (Variance Stabilizing Transformation)**：速度快，适合样本量 n > 30。
> *   **rlog (Regularized Log)**：在样本量较小时更稳健。

```R
# VST 变换示例
vsd <- vst(dds, blind = FALSE)
```

---

### 4.3 探索性分析 (EDA)：发现潜伏的危机

在做差异分析之前，你必须先看两张图：**PCA 图** 和 **样本相关性热图**。

#### 1. PCA (主成分分析)：你的样本聚在一起了吗？
PCA 是诊断实验成败的“第一线”。
*   **理想情况**：同组样本紧密聚类，不同组样本显著分开。
*   **异常情况**：若 Treatment 组样本混入 Control 组，需检查是否存在**标签混淆**或**实验失败**。

#### 2. 样本距离热图 (Sample Distance Heatmap)
展示样本间的整体相似性，有助于发现离群点 (Outliers)。

```R
# 绘制 PCA
plotPCA(vsd, intgroup = "condition")
```

> ⚠️ **避坑指南 (Pitfall)：批次效应 (Batch Effect)**
> 若样本聚类是按“提取日期”而非“实验处理”排列，说明存在批次效应。必须在 DESeq2 模型中加入批次因子（如 `design = ~ batch + condition`）进行校正。
