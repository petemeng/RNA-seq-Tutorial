# 📖 Bulk RNA-seq 进阶专题 (Part 9)

## 第九章：多因素实验设计 —— 交互作用 (Interaction)

这是研究“基因型 x 环境”干扰最核心的方法。

> 📘 **核心概念：交互作用 (Interaction Term)**
> `design = ~ Genotype + Treatment + Genotype:Treatment`
> *   **Interaction (交互项)**：计算的是“差异的差异”。
> *   **核心问题**：突变体对处理的反应，是否与野生型不同？

### ⚠️ 避坑指南：Interaction 的提取
在 DESeq2 中，交互项通常命名为类似 `GenotypeMutant.TreatmentTreated`。只有当这个项显著时，我们才说该基因在不同背景下的响应模式发生了**质的变化**。
