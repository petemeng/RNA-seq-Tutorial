# 📘 Bulk RNA-seq 进阶专题 (Part 7)

## 第七章：WGCNA —— 从“差异基因”到“协同模块”

当样本量足够时（建议 `n >= 15`，理想 `n >= 20`），WGCNA 可以回答 DEGs 难以回答的问题：
哪些基因总是一起上调/下调？这些“模块”是否和处理、时间或表型显著相关？

---

## 7.1 核心概念：WGCNA 在做什么？

WGCNA 的核心是把“基因表达矩阵”转为“基因网络”：

1. 先计算基因两两相关性 `cor(i, j)`。
2. 再做软阈值加权：

$$
a_{ij} = |cor(i,j)|^\beta
$$

其中 `beta`（Power）越大，强相关会被放大得越明显。

3. 基于 `TOM (Topological Overlap Matrix)` 进行聚类，切分成不同颜色模块（module）。
4. 每个模块用一个“特征基因”表示（`Module Eigengene`, 即该模块第一主成分）。

> 💡 **Why over How**
> DEGs 是“单基因显著性”，WGCNA 是“系统级协同关系”。
> 两者结合更接近真实生物过程。

---

## 7.2 输入数据与前置条件

### ✅ 输入建议

- 使用 `DESeq2::vst` 或 `rlog` 转换后的矩阵（不建议直接用原始 counts）。
- 过滤低表达基因（例如平均 normalized count > 10）。
- 去除明显离群样本（先做样本聚类检查）。

### ⚠️ 避坑指南：样本数不足

`n < 12` 时网络结构很不稳定，模块可重复性差。  
如果样本量有限，优先做 DEA/富集，不要过度解读 WGCNA 模块。

---

## 7.3 实战演练：WGCNA 标准流程（R）

```r
library(DESeq2)
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# 1) 准备表达矩阵：gene x sample -> sample x gene
vsd <- vst(dds, blind = FALSE)
expr <- assay(vsd)

# 低表达过滤（可按项目调整阈值）
expr <- expr[rowMeans(counts(dds, normalized = TRUE)) > 10, ]
datExpr <- t(expr)

# 2) 检查坏样本/坏基因
gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# 3) 样本聚类，识别离群点
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering", xlab = "", sub = "")
```

---

## 7.4 Best Practice：软阈值（Power）选择

目标是让网络尽量接近“无尺度拓扑”（scale-free topology）。

```r
powers <- 1:20
sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  networkType = "signed",
  corFnc = "bicor",
  verbose = 5
)

# 一般选择使 scale-free R^2 接近 0.8 的最小 power
softPower <- 8
```

> ⚠️ 避坑指南：不要机械套用 `power = 6`
> 不同数据集分布不同，必须用 `pickSoftThreshold` 实测。

---

## 7.5 构建模块并关联表型

```r
net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE,
  numericLabels = TRUE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# 示例：关联处理因子（按你的 colData 字段替换）
trait <- as.numeric(colData(dds)[rownames(datExpr), "condition"] == "AvrRpm1")
moduleTraitCor <- bicor(MEs, trait, use = "pairwise.complete.obs")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))
```

解释方式：

- 相关系数绝对值越大，模块与表型关联越强。
- `p < 0.05` 是初筛标准，建议结合多重检验和生物学一致性。

---

## 7.6 Hub Gene 提取（模块内关键基因）

```r
moduleOfInterest <- "turquoise"
inModule <- moduleColors == moduleOfInterest

kME <- signedKME(datExpr, MEs)
hub <- rownames(kME)[inModule][order(kME[inModule, paste0("kME", moduleOfInterest)], decreasing = TRUE)]
head(hub, 20)
```

> 💡 **Best Practice**
> Hub gene 不是“因果基因”，它是网络中心节点。
> 建议结合文献、功能注释和实验验证进一步筛选。

---

## 7.7 常见翻车点清单

1. 把强批次效应数据直接拿去做网络，模块会按 batch 分组。
2. 用 Pearson 相关性导致离群点主导网络结构。
3. 基因过滤太松（噪声过多）或太严（信息丢失）。
4. 把模块-表型相关当作因果结论，过度解读。

---

## 7.8 本章小结

WGCNA 适合“样本量较大 + 目标是系统层解释”的场景。  
标准流程是：`vst/rlog -> 选 power -> 建模块 -> 关联表型 -> 提 hub`。  
只要你把样本质量、批次效应和相关性度量控制好，WGCNA 是非常有教学价值和发表价值的进阶方法。

---

## 7.9 本教程实跑代码与结果（PRJDB11848）

本章对应的实跑代码在下游一体化脚本中：

```bash
Rscript scripts/04_downstream_ch4_to_ch9.R
```

结果文件：

- `validation_run_downstream/results/ch7/wgcna_soft_threshold.csv`
- `validation_run_downstream/results/ch7/wgcna_modules.csv`
- `validation_run_downstream/results/ch7/module_trait_cor.csv`
- `validation_run_downstream/results/ch7/module_trait_p.csv`

说明：

- 为保证教程运行时长可控，实跑脚本对构网输入限制为高变 `5000` 基因。