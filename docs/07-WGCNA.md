# 7. WGCNA 共表达网络分析 (Weighted Gene Co-expression Network Analysis)

## 7.1 Motivation

差异表达分析回答的是"哪些基因在两个条件之间显著变化"，但它无法告诉你这些基因之间的关系——哪些基因总是一起上调？哪些基因的表达模式高度协同？这些协同变化的基因模块是否与特定的实验条件或表型相关？

WGCNA（Weighted Gene Co-expression Network Analysis）正是为了回答这类"系统层"问题而设计的 [^langfelder2008]。它将基因表达矩阵转化为加权共表达网络，通过聚类识别出协同表达的基因模块（module），然后将模块与外部性状（如处理条件、时间点、基因型）进行关联分析。与 DEA 的"单基因显著性"视角不同，WGCNA 提供的是"基因组级别的协同关系"视角——两者结合更接近真实的生物学过程。

需要注意的是，WGCNA 对样本量有较高要求。当样本数少于 15 时，网络结构不稳定，模块的可重复性差。PRJDB11848 的 36 个样本刚好满足这个要求，使得 WGCNA 分析具有合理的统计基础。

---

## 7.2 环境与数据 (Environment and Data)

本章在 R 中完成，需要以下包：

- WGCNA（网络构建与模块检测）
- DESeq2（VST 变换）
- pheatmap（热图可视化）

输入：[第 4 章](04-Normalization-and-Exploratory-Analysis.md) 生成的 VST 变换矩阵
输出：模块分配表、模块-性状关联矩阵、hub gene 列表

---

## 7.3 输入数据准备

WGCNA 的输入应该是经过方差稳定化变换的表达矩阵，而不是原始计数。同时需要过滤低表达和低变异基因，减少噪声对网络结构的干扰：

```r
library(DESeq2)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 1) VST 变换
vsd <- vst(dds, blind = FALSE)
expr <- assay(vsd)

# 2) 过滤低表达基因（保留平均 normalized count > 10 的基因）
expr <- expr[rowMeans(counts(dds, normalized = TRUE)) > 10, ]

# 3) 转置为 sample × gene 矩阵（WGCNA 要求的格式）
datExpr <- t(expr)

# 4) 检查并移除坏样本/坏基因
gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
```

!!! warning "避坑指南"
    样本数不足（n < 15）时网络结构很不稳定，模块可重复性差。如果你的样本量有限，优先做 DEA 和富集分析，不要过度解读 WGCNA 模块。

### 样本聚类：识别离群点

在构建网络之前，先做样本聚类检查是否有明显的离群样本：

```r
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", xlab = "", sub = "")
```

如果某个样本与其他所有样本的距离都很远，应在构建网络前将其移除——离群样本会严重扭曲相关性估计。

---

## 7.4 软阈值（Power）选择

WGCNA 通过软阈值加权将相关性矩阵转化为邻接矩阵：$a_{ij} = |cor(i,j)|^\beta$。Power（$\beta$）的选择目标是让网络尽量接近无尺度拓扑（scale-free topology）：

```r
powers <- 1:20
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",     # signed 网络区分正负相关
  corFnc = "bicor",           # bicor 对离群值更稳健
  verbose = 5
)

# 选择使 scale-free R² 接近 0.8 的最小 power
softPower <- sft$powerEstimate  # 或手动选择
```

!!! warning
    不要机械套用 `power = 6`。不同数据集的分布不同，必须用 `pickSoftThreshold` 实测。如果没有任何 power 能让 R² 达到 0.8，可能说明数据中存在强烈的批次效应或样本异质性，需要先处理这些问题。

---

## 7.5 构建模块并关联性状

```r
net <- blockwiseModules(
  datExpr,
  power            = softPower,
  networkType      = "signed",
  TOMType          = "signed",
  corType          = "bicor",       # 比 Pearson 更稳健
  minModuleSize    = 30,            # 模块最小基因数
  mergeCutHeight   = 0.25,          # 相似模块合并阈值
  pamRespectsDendro = FALSE,
  numericLabels    = TRUE,
  verbose          = 3
)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs  # Module Eigengenes

# 关联处理条件
trait <- as.numeric(colData(dds)[rownames(datExpr), "condition"] == "AvrRpm1")
moduleTraitCor <- bicor(MEs, trait, use = "pairwise.complete.obs")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))
```

!!! tip "Recommendation"
    我们推荐使用 **signed 网络 + bicor 相关性**。signed 网络区分正相关和负相关基因（它们可能参与不同的生物学过程），bicor（biweight midcorrelation）比 Pearson 相关性对离群值更稳健 [^langfelder2008]。

解读模块-性状关联：

- 相关系数绝对值越大，模块与性状的关联越强。
- p < 0.05 是初筛标准，但建议结合多重检验校正和生物学一致性综合判断。
- 关联不等于因果——模块与处理条件相关，不代表模块中的基因直接参与了处理响应。

---

## 7.6 Hub Gene 提取

Hub gene 是模块内连接度最高的基因，通常被认为是该模块功能的核心代表：

```r
moduleOfInterest <- "turquoise"
inModule <- moduleColors == moduleOfInterest

kME <- signedKME(datExpr, MEs)
hub <- rownames(kME)[inModule][
  order(kME[inModule, paste0("kME", moduleOfInterest)], decreasing = TRUE)
]
head(hub, 20)
```

!!! note
    Hub gene 是网络中心节点，不等于"因果基因"。它代表的是"与模块整体表达模式最一致的基因"。建议结合文献、功能注释和实验验证进一步筛选候选基因，而不是直接将 hub gene 等同于关键调控因子。

---

## 7.7 常见翻车点

!!! warning "常见翻车点"
    1. **强批次效应未处理**：批次效应会让模块按 batch 分组而非按生物学功能分组。构建网络前应先检查 PCA，必要时使用批次校正后的矩阵（仅用于 WGCNA 输入，不用于 DESeq2）。
    2. **用 Pearson 相关性**：Pearson 对离群值极其敏感，一个异常样本就可能扭曲整个网络结构。推荐使用 bicor。
    3. **基因过滤太松或太严**：太松会引入大量噪声基因，模块结构模糊；太严会丢失重要信息。建议保留高变异基因（如 top 5,000-10,000）。
    4. **把模块-性状相关当因果结论**：相关性分析只能提示关联，不能证明因果。过度解读是 WGCNA 论文中最常见的问题。

---

## 7.8 Key Takeaways

!!! abstract "本章要点"
    1. **WGCNA 提供系统层视角**——它回答的是"哪些基因协同变化"，与 DEA 的"单基因显著性"互补。
    2. **样本量是硬性要求**——n < 15 时网络不稳定，不建议使用 WGCNA [^langfelder2008]。
    3. **signed 网络 + bicor 是推荐配置**——区分正负相关，对离群值稳健。
    4. **软阈值必须实测**——不同数据集需要不同的 power，不要套用固定值。
    5. **Hub gene ≠ 因果基因**——它是网络中心节点，需要结合功能注释和实验验证进一步筛选。
    6. **模块-性状关联是相关性，不是因果性**——谨慎解读，避免过度叙事。

---

## 7.9 References

[^langfelder2008]: Langfelder P, Horvath S. (2008) WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics* 9:559. DOI: [10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559)

