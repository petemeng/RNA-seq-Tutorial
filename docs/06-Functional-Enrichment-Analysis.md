# 📖 Bulk RNA-seq 数据分析最佳实践 (Part 6)

## 第六章：功能富集分析 —— 从基因列表到生物学意义

功能富集分析旨在将差异基因 (DEGs) 归类到已知的生物学模块中。目前最主流的工具有两类：**ORA (过表达分析)** 和 **GSEA (基因集富集分析)**。

### 6.1 过表达分析 (Over-Representation Analysis, ORA)

> 📘 **核心概念：超几何分布 (Hypergeometric Test)**
> ORA 通过 Fisher 精确检验判断你的差异基因在某个通路中的占比是否显著高于背景分布。

#### ⚠️ 避坑指南：ORA 的局限性
1.  **阈值敏感性**：严重依赖你设置的差异筛选阈值。
2.  **忽略微弱信号**：它只看“显著”的基因，会损失具有协同趋势的微弱信号。

---

### 6.2 基因集富集分析 (Gene Set Enrichment Analysis, GSEA)

GSEA 是目前最推崇的方法，它**不需要**设置人为的差异筛选阈值。

> 📘 **核心概念：排名驱动 (Rank-based)**
> GSEA 对全基因组按差异倍数排名。如果某个通路的基因集中在排名的两端，则判定为富集。
> **优势**：它能捕捉到“微弱但系统”的生物学变化。

---

### 6.3 💡 Best Practice：ID 转换与物种匹配

新手最容易翻车的地方就是 **Gene ID 转换**。大多数包使用的是 **Entrez ID**，而原始矩阵通常是 **Ensembl ID** 或 **Symbol**。

**最佳实践：** 使用 `org.Hs.eg.db` 等包进行转换。转换率低于 70% 通常意味着注释包选错了。

---

### 🛠️ 实战演练：使用 clusterProfiler 

```R
library(clusterProfiler)
library(org.Hs.eg.db)

# 1. 运行 GO 富集分析
ego <- enrichGO(gene          = gene_list_entrez,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH")

# 2. 可视化：点图 (Dotplot)
dotplot(ego, showCategory=20) + ggtitle("GO Biological Process Enrichment")
```

---

### 6.4 可视化：如何讲好生物学故事？

*   **Dotplot**：同时展示显著性与基因比例。
*   **GSEA Plot**：展示富集得分的动态过程，是证明通路激活的最有力证据。
*   **Enrichment Map**：揭示不同富集通路之间的功能重合与耦合。
