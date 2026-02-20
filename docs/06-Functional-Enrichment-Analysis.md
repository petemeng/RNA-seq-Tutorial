# 6. 功能富集分析 (Functional Enrichment Analysis)

## 6.1 Motivation

差异表达分析给你一份基因列表——1,541 个显著 DEG。但这份列表本身不回答生物学问题。你需要知道的是：这些基因在功能上是否存在系统性的偏向？它们是否集中在某些生物学通路或过程中？这就是功能富集分析要回答的问题。

富集分析有两种根本不同的思路。ORA（Over-Representation Analysis）问的是"我的 DEG 列表中，某个通路的基因是否比随机预期多？"——它依赖于你预先设定的 DEG 阈值。GSEA（Gene Set Enrichment Analysis）问的是"在按 fold change 排序的全基因列表中，某个通路的基因是否集中出现在排序的顶端或底端？"——它不需要预先设定阈值，对"弱但一致"的信号更敏感 [^subramanian2005]。两种方法各有优势，在同一数据集上同时使用可以增强结论的稳健性。

本章将在 PRJDB11848 的 WT AvrRpm1 vs mock DEG 结果上，分别完成 GO 的 ORA 和 GSEA 分析，并讨论结果解读的关键原则。

---

## 6.2 环境与数据 (Environment and Data)

本章在 R 中完成，需要以下包：

- clusterProfiler（ORA 和 GSEA 分析）
- org.At.tair.db（拟南芥注释数据库）
- AnnotationDbi（ID 转换）
- ggplot2（可视化）

输入：[第 5 章](05-Differential-Expression-Analysis.md) 生成的 DEG 结果表（含 log2FoldChange 和 padj）

---

## 6.3 富集前关键步骤：ID 转换与背景集

拟南芥的基因 ID 通常是 TAIR 格式（如 `AT1G01010`），但 clusterProfiler 的部分函数需要 ENTREZ ID。ID 转换是富集分析的必要前置步骤。

### ID 转换（TAIR → ENTREZ）

```r
library(clusterProfiler)
library(org.At.tair.db)

# 显著 DEG 的 TAIR ID
sig_gene_tair <- rownames(subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1))

id_map <- bitr(
  sig_gene_tair,
  fromType = "TAIR",
  toType   = "ENTREZID",
  OrgDb    = org.At.tair.db
)

gene_entrez <- unique(id_map$ENTREZID)
```

### 背景集（universe）定义

背景集应定义为"本次实验中可被检测到的基因集合"（即通过低表达过滤后的所有基因），而不是全基因组。背景集定义错误会直接改变富集的显著性。

```r
universe_tair <- rownames(res_df)  # 所有检测到的基因
universe_map <- bitr(
  universe_tair,
  fromType = "TAIR",
  toType   = "ENTREZID",
  OrgDb    = org.At.tair.db
)
universe_entrez <- unique(universe_map$ENTREZID)
```

!!! warning
    背景集定义是 ORA 中最容易出错的步骤。如果你用全基因组（~27,000 基因）作为背景，但实际只检测到 15,000 个基因，富集的显著性会被人为夸大。始终使用实验中实际检测到的基因作为背景集。

---

## 6.4 ORA 实战（GO）

ORA 使用超几何检验判断 DEG 列表中某个 GO term 的基因是否显著过多：

```r
ego <- enrichGO(
  gene         = gene_entrez,
  universe     = universe_entrez,
  OrgDb        = org.At.tair.db,
  keyType      = "ENTREZID",
  ont          = "BP",              # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

head(as.data.frame(ego), 10)
```

<figure markdown>
  ![GO ORA](assets/validated_case/ch6_go_ora_top15.png){ width="680" }
  <figcaption>Fig 6.1 — WT 中 AvrRpm1 vs mock 显著 DEG 的 GO Biological Process ORA 结果（Top 15）。富集到的通路以防御响应和免疫相关过程为主，与 AvrRpm1 触发的 ETI 免疫响应一致。</figcaption>
</figure>

---

## 6.5 GSEA 实战（GO）

GSEA 不需要预先筛选 DEG，而是使用全基因的排序向量。我们用收缩后的 log2FoldChange 作为排序指标：

```r
# 构建排序向量
gene_rank <- res_df$log2FoldChange
names(gene_rank) <- rownames(res_df)
gene_rank <- sort(gene_rank, decreasing = TRUE)

# ID 转换
rank_map <- bitr(
  names(gene_rank),
  fromType = "TAIR",
  toType   = "ENTREZID",
  OrgDb    = org.At.tair.db
)
rank_map <- rank_map[!duplicated(rank_map$ENTREZID), ]

gene_rank_entrez <- gene_rank[rank_map$TAIR]
names(gene_rank_entrez) <- rank_map$ENTREZID
gene_rank_entrez <- sort(gene_rank_entrez, decreasing = TRUE)
```

运行 GSEA：

```r
gse <- gseGO(
  geneList      = gene_rank_entrez,
  OrgDb         = org.At.tair.db,
  ont           = "BP",
  keyType       = "ENTREZID",
  pAdjustMethod = "BH",
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  verbose       = FALSE
)

head(as.data.frame(gse), 10)
```

<figure markdown>
  ![GO GSEA](assets/validated_case/ch6_go_gsea_top15.png){ width="680" }
  <figcaption>Fig 6.2 — WT 中 AvrRpm1 vs mock 的 GO GSEA 结果（Top 15）。正的 NES（Normalized Enrichment Score）表示该通路的基因整体上调，负的 NES 表示整体下调。免疫和防御相关通路显著富集在上调端，与 ORA 结果相互印证。</figcaption>
</figure>

---

## 6.6 KEGG 富集（拟南芥）

拟南芥的 KEGG 物种代码是 `ath`：

```r
ekk <- enrichKEGG(
  gene         = gene_entrez,
  universe     = universe_entrez,
  organism     = "ath",
  pAdjustMethod = "BH"
)
```

!!! warning "避坑指南"
    物种代码填错会导致结果异常少或为空。拟南芥是 `ath`，不是 `at` 或 `arabidopsis`。可以在 [KEGG Organism](https://www.genome.jp/kegg/catalog/org_list.html) 页面查询正确的物种代码。

---

## 6.7 如何解读富集结果

富集分析的结果需要谨慎解读，避免过度叙事：

1. **同时关注显著性和方向**：GSEA 的 NES 告诉你通路是整体上调还是下调；ORA 只告诉你"过多"，不告诉你方向。
2. **合并高度重叠的条目**：GO 的层级结构导致很多相似的 term 同时显著（如 "defense response"、"response to biotic stimulus"、"immune response"）。关注独立的生物学主题，而不是逐条列举。
3. **回到具体基因**：通路层面的结论需要回到关键基因的表达变化来支撑。输出"候选通路 + 关键基因 + 表达证据"三联表是最有说服力的呈现方式。
4. **ORA 和 GSEA 交叉验证**：两种方法都支持的结论更可靠。

!!! tip "Recommendation"
    在报告中同时给出 ORA 和 GSEA 结果，增强结论稳健性。ORA 直观易解释，GSEA 对弱信号更敏感——两者互补。

---

## 6.8 常见翻车点

!!! warning "常见翻车点"
    1. **DEG 列表太小做 ORA**：只有几十个 DEG 时，ORA 结果极不稳定，换一个阈值结论就变了。此时 GSEA 是更好的选择。
    2. **忽略 ID 转换损失**：TAIR → ENTREZ 的映射不是 100% 的，部分基因会丢失。检查转换率，如果损失超过 20%，需要评估对结果的影响。
    3. **背景集设置不合理**：用全基因组而非实验检测到的基因作为背景，会夸大显著性。
    4. **只贴通路图不回到基因**：审稿人会问"这个通路里具体哪些基因变了？变化方向是什么？"——你需要准备好答案。

---

## 6.9 Key Takeaways

!!! abstract "本章要点"
    1. **富集分析是从基因列表到生物学解释的桥梁**——它回答的是"这些基因在功能上是否存在系统性偏向"。
    2. **ORA 和 GSEA 是互补的两种方法**——ORA 依赖阈值但直观，GSEA 不依赖阈值且对弱信号敏感，同时使用增强结论稳健性。
    3. **背景集定义是 ORA 最关键的参数**——必须使用实验中实际检测到的基因，而非全基因组 [^subramanian2005]。
    4. **ID 转换损失需要监控**——TAIR → ENTREZ 的映射不完美，损失过大会影响结果可靠性。
    5. **富集结果需要回到具体基因**——通路层面的结论必须有关键基因的表达变化作为支撑。
    6. **合并重叠的 GO term**——GO 的层级结构会产生大量冗余条目，关注独立的生物学主题而非逐条列举。

---

## 6.10 References

[^subramanian2005]: Subramanian A, Tamayo P, Mootha VK, et al. (2005) Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *PNAS* 102(43):15545-15550. DOI: [10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)

[^yu2012]: Yu G, Wang LG, Han Y, He QY. (2012) clusterProfiler: an R package for comparing biological themes among gene clusters. *OMICS* 16(5):284-287. DOI: [10.1089/omi.2011.0118](https://doi.org/10.1089/omi.2011.0118)

---

## 6.11 Contributors

### Authors

- *(待补充)*

### Reviewers

- *(待补充)*
