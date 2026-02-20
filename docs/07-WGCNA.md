# 📖 Bulk RNA-seq 进阶专题 (Part 7)

## 第七章：WGCNA —— 发现基因共表达模块 (Co-expression Network)

当你拥有大量样本（n > 15）时，仅仅看差异基因 (DEGs) 已经不够了。你需要一种方法来发现**协同工作**的基因集合，并将其与临床/表型指标关联。

> 📘 **核心概念：加权基因共表达网络分析 (WGCNA)**
> 1. **加权 (Weighted)**：使用“软阈值 (Power)”放大强信号，保留弱信号。
> 2. **模块 (Modules)**：将表达模式相似的基因聚类为模块。
> 3. **特征基因 (Eigengene)**：模块内基因表达的“第一主成分”，代表该模块的整体水平。

### 💡 Best Practice：软阈值 (Power) 选择
必须确保网络符合 **“无尺度分布 (Scale-free Topology)”**。
*   **操作**：运行 `pickSoftThreshold`。
*   **标准**：选择使 $R^2$ 平台期达到 0.8-0.9 的最小 Power。

### ⚠️ 避坑指南：Biweight Midcorrelation (bicor)
Pearson 相关性对离群点极其敏感。在构建网络时，务必使用 **`bicor`**，它能有效抵御离群值的干扰。
