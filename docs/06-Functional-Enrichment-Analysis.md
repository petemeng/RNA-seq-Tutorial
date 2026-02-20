# Bulk RNA-seq 数据分析最佳实践 (Part 6)

## 第六章：功能富集分析 (Functional Enrichment Analysis)

> 核心原则：富集分析不是“列通路清单”，而是建立可验证的生物学解释框架。

---

## 6.1 ORA 与 GSEA 的本质区别

### ORA（Over-Representation Analysis）

- 输入：显著基因列表（离散集合）
- 统计：超几何检验/Fisher 精确检验
- 特点：直观，但依赖阈值

### GSEA（Gene Set Enrichment Analysis）

- 输入：全基因排序向量（连续信息）
- 统计：富集得分在排序列表上的偏移
- 特点：对“弱但一致”的变化更敏感

> Best Practice：报告中常同时给出 ORA 与 GSEA，增强结论稳健性。

---

## 6.2 富集前的关键步骤：ID 转换与背景集定义

### 1) ID 转换

```r
library(clusterProfiler)
library(org.Hs.eg.db)

sig_gene_symbols <- rownames(subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1))

id_map <- bitr(
  sig_gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

gene_entrez <- unique(id_map$ENTREZID)
```

### 2) 背景集（universe）

背景集应为“本次实验中有机会被检测到的基因”，而不是全基因组。

```r
universe_symbols <- rownames(res_df)
universe_map <- bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_entrez <- unique(universe_map$ENTREZID)
```

> 避坑指南：背景集错误会直接改变富集显著性，是高风险环节。

---

## 6.3 ORA 实战（GO）

```r
ego <- enrichGO(
  gene = gene_entrez,
  universe = universe_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

head(as.data.frame(ego), 10)
```

可视化：

```r
dotplot(ego, showCategory = 20)
barplot(ego, showCategory = 20)
```

---

## 6.4 GSEA 实战（GO）

构建排序向量（示例用收缩后的 `log2FC`）：

```r
gene_rank <- res_df$log2FoldChange
names(gene_rank) <- rownames(res_df)
gene_rank <- sort(gene_rank, decreasing = TRUE)

rank_map <- bitr(names(gene_rank), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
rank_map <- rank_map[!duplicated(rank_map$ENTREZID), ]

gene_rank_entrez <- gene_rank[rank_map$SYMBOL]
names(gene_rank_entrez) <- rank_map$ENTREZID
gene_rank_entrez <- sort(gene_rank_entrez, decreasing = TRUE)
```

运行 GSEA：

```r
gse <- gseGO(
  geneList = gene_rank_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

head(as.data.frame(gse), 10)
```

可视化：

```r
ridgeplot(gse)
# gseaplot2(gse, geneSetID = 1)
```

---

## 6.5 KEGG 与物种代码

做 KEGG 时要确认物种缩写：

- 人类：`hsa`
- 小鼠：`mmu`
- 拟南芥：`ath`

```r
ekk <- enrichKEGG(
  gene = gene_entrez,
  universe = universe_entrez,
  organism = "hsa",
  pAdjustMethod = "BH"
)
```

> 避坑指南：物种代码错了会导致结果异常少或完全为空。

---

## 6.6 如何解释富集结果更可靠

1. 看 `NES/富集方向`，不仅看显著性。
2. 合并同义/高度重叠通路，避免重复叙事。
3. 与实验背景对照，优先解释与研究问题直接相关的条目。
4. 输出“候选通路 + 关键基因 + 证据图”三联表。

---

## 6.7 常见翻车点

1. 用非常小的 DEG 列表做 ORA，结果不稳定。
2. 忽略 ID 映射损失，导致有效基因数过低。
3. 用全基因组当 universe，夸大显著性。
4. 只贴图不回到关键基因表达变化。

---

## 6.8 本章检查清单

- ORA 与 GSEA 至少完成一种，并保留参数。
- ID 转换率与背景集定义可追溯。
- 富集结果包含方向、显著性、代表基因。
- 结果图与表可支撑论文叙事。

基础流程到此闭环完成。
