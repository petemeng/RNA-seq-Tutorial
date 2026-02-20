# 📘 Bulk RNA-seq 进阶专题 (Part 9)

## 第九章：多因素实验设计 —— 交互作用与可解释对比

真实实验很少只有一个变量。常见组合包括：

- 基因型（WT/Mut）
- 处理（mock/AvrRpm1）
- 时间（0h/1h/6h）
- 批次（Batch1/Batch2）

多因素模型的目标不是“把变量都塞进去”，而是回答清楚：
主效应是什么？交互效应是什么？我到底在比较什么？

---

## 9.1 核心概念：主效应 vs 交互效应

以两因素为例：

```text
design = ~ Genotype + Treatment + Genotype:Treatment
```

- **Genotype 主效应**：在参考处理条件下，Mut 与 WT 的差异。
- **Treatment 主效应**：在参考基因型下，AvrRpm1 与 mock 的差异。
- **Interaction 交互效应**：
  `(Mut_AvrRpm1 - Mut_mock) - (WT_AvrRpm1 - WT_mock)`

> 💡 交互显著 = “差异的差异”显著。
> 它表示某处理效应依赖于另一因素的水平。

---

## 9.2 什么时候必须建交互项？

如果你的科学问题是：

1. “突变体对 AvrRpm1 的响应是否不同于 WT？”
2. “药物效应是否依赖基因型/性别/组织？”

那就必须包含交互项，而不是只做分组拆分比较。

---

## 9.3 实战演练：DESeq2 多因素建模

```r
library(DESeq2)

colData$genotype <- relevel(factor(colData$genotype), ref = "WT")
colData$condition <- relevel(factor(colData$condition), ref = "mock")
colData$batch <- factor(colData$batch)

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = colData,
  design = ~ batch + genotype + condition + genotype:condition
)

dds <- DESeq(dds)
resultsNames(dds)
```

> ⚠️ 避坑指南：一定先看 `resultsNames(dds)`。
> 交互项名字会随因子命名和参考水平变化，不要硬写固定字符串。

---

## 9.4 如何正确提取结果？

```r
# 1) WT 中的处理效应（参考基因型）
res_treat_in_wt <- results(dds, name = "condition_AvrRpm1_vs_mock")

# 2) 交互效应：clf 对处理的额外响应
res_interaction <- results(dds, name = "genotypeclf.conditionAvrRpm1")

# 3) clf 中的处理总效应 = 主效应 + 交互项
res_treat_in_clf <- results(
  dds,
  list(c("condition_AvrRpm1_vs_mock", "genotypeclf.conditionAvrRpm1"))
)
```

解释逻辑：

- `res_treat_in_wt` 讲 WT 的处理响应。
- `res_interaction` 讲 Mut 相对 WT 的“额外变化”。
- `res_treat_in_mut` 才是 Mut 的处理总响应。

---

## 9.5 模型质量检查与可视化建议

1. 先看设计平衡性：每个 cell（如 `WT-mock`、`WT-AvrRpm1`、`clf-mock`、`clf-AvrRpm1`）至少建议 3 个重复。
2. 画交互图（interaction plot）检查方向是否符合预期。
3. 对显著交互基因画分组箱线图，避免“统计显著但生物学不可解释”。

示例（单基因可视化）：

```r
library(ggplot2)
library(tibble)

vsd <- vst(dds, blind = FALSE)
gene_id <- "AT1G01010"

df_plot <- tibble(
  expr = assay(vsd)[gene_id, ],
  genotype = colData(dds)$genotype,
  treatment = colData(dds)$treatment
)

ggplot(df_plot, aes(x = treatment, y = expr, color = genotype, group = genotype)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point", size = 3)
```

---

## 9.6 常见翻车点清单

1. **混杂设计**：例如 `batch` 与 `treatment` 完全重合，模型无法识别真实效应。
2. **参考水平没设定**：导致结果名称和解释完全跑偏。
3. **把主效应当总效应解读**：在含交互模型中这是常见错误。
4. **拆组分别跑 DESeq2**：损失统计效率，也难以统一解释。

---

## 9.7 本章小结

多因素分析的核心不是“复杂”，而是“可解释”。
只要你明确参考水平、对比逻辑和交互含义，DESeq2 可以稳定回答大多数实验设计问题。

---

## 9.8 本教程实跑代码与结果（PRJDB11848）

### 代码：Chapter 9 实跑片段（对应 `scripts/04_downstream_ch4_to_ch9.R`）

```r
# 多因素模型结果提取
res_treat_in_wt <- results(dds, name = "condition_AvrRpm1_vs_mock")
res_interaction <- results(dds, name = "genotypeclf.conditionAvrRpm1")
res_treat_in_clf <- results(dds, list(c("condition_AvrRpm1_vs_mock", "genotypeclf.conditionAvrRpm1")))

write.csv(as.data.frame(res_treat_in_wt), "validation_run_downstream/results/ch9/res_treat_in_wt.csv")
write.csv(as.data.frame(res_interaction), "validation_run_downstream/results/ch9/res_interaction.csv")
write.csv(as.data.frame(res_treat_in_clf), "validation_run_downstream/results/ch9/res_treat_in_clf.csv")
```

### 代码：交互项火山图（对应 `scripts/05_generate_case_figures.R`）

```r
library(readr)
library(ggplot2)

inter <- read_csv("validation_run_downstream/results/ch9/res_interaction.csv", show_col_types = FALSE)
inter$neglog10padj <- -log10(pmax(inter$padj, 1e-300))
inter$sig <- ifelse(!is.na(inter$padj) & inter$padj < 0.05 & abs(inter$log2FoldChange) > 1, "Sig", "NS")

p_inter <- ggplot(inter, aes(log2FoldChange, neglog10padj, color = sig)) +
  geom_point(alpha = 0.6, size = 1.0) +
  scale_color_manual(values = c("Sig" = "#6a3d9a", "NS" = "grey70")) +
  theme_bw(base_size = 12) +
  labs(title = "Interaction Volcano (clf:AvrRpm1)", x = "log2FoldChange", y = "-log10(padj)")

grDevices::png("docs/assets/validated_case/ch9_interaction_volcano.png", width = 6.8, height = 5.2, units = "in", res = 150)
print(p_inter)
grDevices::dev.off()
```

### 代码：验收命令

```bash
ls artifacts/prjdb11848/results/ch9/res_treat_in_wt.csv \
   artifacts/prjdb11848/results/ch9/res_interaction.csv \
   artifacts/prjdb11848/results/ch9/res_treat_in_clf.csv

wc -l artifacts/prjdb11848/results/ch9/res_treat_in_wt.csv \
      artifacts/prjdb11848/results/ch9/res_interaction.csv \
      artifacts/prjdb11848/results/ch9/res_treat_in_clf.csv

awk -F, 'NR>1 && $7!="" && $7<0.05 && (($3+0>1)||($3+0<-1)) {n++} END{print n+0}' artifacts/prjdb11848/results/ch9/res_interaction.csv
```

### 输出结果

```text
/home/data/t060551/Codex/RNA-seq-Tutorial/artifacts/prjdb11848/results/ch9/res_interaction.csv
/home/data/t060551/Codex/RNA-seq-Tutorial/artifacts/prjdb11848/results/ch9/res_treat_in_clf.csv
/home/data/t060551/Codex/RNA-seq-Tutorial/artifacts/prjdb11848/results/ch9/res_treat_in_wt.csv
  17453 /home/data/t060551/Codex/RNA-seq-Tutorial/artifacts/prjdb11848/results/ch9/res_treat_in_wt.csv
  17453 /home/data/t060551/Codex/RNA-seq-Tutorial/artifacts/prjdb11848/results/ch9/res_interaction.csv
  17453 /home/data/t060551/Codex/RNA-seq-Tutorial/artifacts/prjdb11848/results/ch9/res_treat_in_clf.csv
  52359 total
56
```

交互项火山图：

![Interaction Volcano](assets/validated_case/ch9_interaction_volcano.png)
