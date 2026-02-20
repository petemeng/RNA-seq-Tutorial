# Bulk RNA-seq 数据分析最佳实践 (Part 6)

## 第六章：功能富集分析（拟南芥）

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

> Best Practice：在报告中同时给出 ORA 与 GSEA，增强结论稳健性。

---

## 6.2 富集前关键步骤：ID 转换与背景集定义

拟南芥常见基因 ID 为 TAIR（如 `AT1G01010`）。

### 1) ID 转换（TAIR -> ENTREZ）

```r
library(clusterProfiler)
library(org.At.tair.db)

sig_gene_tair <- rownames(subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1))

id_map <- bitr(
  sig_gene_tair,
  fromType = "TAIR",
  toType = "ENTREZID",
  OrgDb = org.At.tair.db
)

gene_entrez <- unique(id_map$ENTREZID)
```

### 2) 背景集（universe）

背景集应定义为“本次实验中可被检测到的基因集合”，不要直接用全基因组。

```r
universe_tair <- rownames(res_df)
universe_map <- bitr(
  universe_tair,
  fromType = "TAIR",
  toType = "ENTREZID",
  OrgDb = org.At.tair.db
)
universe_entrez <- unique(universe_map$ENTREZID)
```

> 避坑指南：背景集定义错误会直接改变显著性，属于高风险步骤。

---

## 6.3 ORA 实战（GO）

```r
ego <- enrichGO(
  gene = gene_entrez,
  universe = universe_entrez,
  OrgDb = org.At.tair.db,
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

构建排序向量（示例使用收缩后的 `log2FoldChange`）：

```r
gene_rank <- res_df$log2FoldChange
names(gene_rank) <- rownames(res_df)
gene_rank <- sort(gene_rank, decreasing = TRUE)

rank_map <- bitr(
  names(gene_rank),
  fromType = "TAIR",
  toType = "ENTREZID",
  OrgDb = org.At.tair.db
)
rank_map <- rank_map[!duplicated(rank_map$ENTREZID), ]

gene_rank_entrez <- gene_rank[rank_map$TAIR]
names(gene_rank_entrez) <- rank_map$ENTREZID
gene_rank_entrez <- sort(gene_rank_entrez, decreasing = TRUE)
```

运行 GSEA：

```r
gse <- gseGO(
  geneList = gene_rank_entrez,
  OrgDb = org.At.tair.db,
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

## 6.5 KEGG（拟南芥）

拟南芥 KEGG 物种代码是 `ath`：

```r
ekk <- enrichKEGG(
  gene = gene_entrez,
  universe = universe_entrez,
  organism = "ath",
  pAdjustMethod = "BH"
)
```

> 避坑指南：物种代码填错会导致结果异常少或为空。

---

## 6.6 如何解释富集结果更可靠

1. 解释时同时关注显著性和方向（NES / enrichment direction）。
2. 合并高度重叠条目，避免重复叙事。
3. 回到拟南芥免疫背景，优先解释与 flg22 响应相关通路。
4. 输出“候选通路 + 关键基因 + 证据图”三联表。

---

## 6.7 常见翻车点

1. 用过小 DEG 列表做 ORA，结果不稳定。
2. 忽略 TAIR -> ENTREZ 映射损失，导致有效基因过少。
3. universe 设置不合理，夸大显著性。
4. 只贴通路图，不回到关键基因表达变化。

---

## 6.8 本章检查清单

- ORA 与 GSEA 至少完成一种，并保留参数。
- TAIR 转换率与背景集定义可追溯。
- 结果包含方向、显著性、代表基因。
- 图表和表格可直接支撑论文叙事。

基础流程到此闭环完成。

---

## 6.9 本教程实跑代码与结果（PRJDB11848）

本章对应的实跑代码在下游一体化脚本中：

```bash
Rscript scripts/04_downstream_ch4_to_ch9.R
Rscript scripts/05_generate_case_figures.R
```

结果文件：

- `validation_run_downstream/results/ch6/GO_ORA_WT.csv`
- `validation_run_downstream/results/ch6/GO_GSEA_WT.csv`

GO ORA Top terms：

![GO ORA](assets/validated_case/ch6_go_ora_top15.png)

GO GSEA Top terms：

![GO GSEA](assets/validated_case/ch6_go_gsea_top15.png)
