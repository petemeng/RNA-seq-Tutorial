# Bulk RNA-seq 数据分析最佳实践 (Part 5)

## 第五章：差异表达分析 (Differential Expression Analysis)

本章按固定格式组织：先讲原理，再给代码，再给代码输出结果。

---

## 5.1 原理：DEG 的核心是“带离散度建模的统计检验”

- RNA-seq 是计数型数据，通常存在过离散
- DESeq2 用负二项模型估计 size factor、dispersion，再做检验
- 结果解释要同时看显著性（`padj`）和效应量（`log2FoldChange`）

---

## 5.2 代码：先构建 `dds`（含 `tximport`）

```r
library(readr)
library(tximport)
library(DESeq2)

base_dir <- "validation_run_downstream"

samples <- read_csv(file.path(base_dir, "metadata", "samplesheet.csv"), show_col_types = FALSE)
samples$genotype <- factor(samples$genotype, levels = c("WT", "clf"))
samples$condition <- factor(samples$condition, levels = c("mock", "AvrRpm1"))
samples$time <- factor(samples$time, levels = c("0h", "0.5h", "3h"))
samples$batch <- factor(samples$batch)

files <- file.path(base_dir, "data", "quant", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

tx2gene <- read_tsv(
  file.path(base_dir, "metadata", "tx2gene.tsv"),
  col_names = c("TXNAME", "GENEID"),
  show_col_types = FALSE
)

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = as.data.frame(samples),
  design = ~ genotype + time + condition + genotype:condition
)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)
```

---

## 5.3 代码：运行 DESeq2 并导出 DEG

```r
results_names <- resultsNames(dds)
writeLines(results_names, file.path(base_dir, "results/ch5/results_names.txt"))

# WT: AvrRpm1 vs mock（LFC shrink）
coef_wt <- grep("^condition_", results_names, value = TRUE)[1]
res_wt <- lfcShrink(dds, coef = coef_wt, type = "apeglm")
res_wt_df <- as.data.frame(res_wt)
res_wt_df <- res_wt_df[order(res_wt_df$padj), ]
write.csv(res_wt_df, file.path(base_dir, "results/ch5/DEG_WT_AvrRpm1_vs_mock_shrunk.csv"))

deg_wt <- subset(res_wt_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg_wt, file.path(base_dir, "results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv"))

# 交互项（clf 相对 WT 的额外处理响应）
res_inter <- results(dds, name = "genotypeclf.conditionAvrRpm1")
res_inter_df <- as.data.frame(res_inter)
res_inter_df <- res_inter_df[order(res_inter_df$padj), ]
write.csv(res_inter_df, file.path(base_dir, "results/ch5/DEG_interaction_clf_vs_WT.csv"))

deg_inter <- subset(res_inter_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg_inter, file.path(base_dir, "results/ch5/DEG_interaction_clf_vs_WT_sig.csv"))
```

---

## 5.4 代码：火山图和 MA 图（直接绘图代码）

```r
library(ggplot2)

fig_dir <- "docs/assets/validated_case"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(plot_obj, file, width, height, dpi = 150) {
  grDevices::png(filename = file, width = width, height = height, units = "in", res = dpi)
  print(plot_obj)
  grDevices::dev.off()
}

deg <- read_csv(file.path(base_dir, "results/ch5/DEG_WT_AvrRpm1_vs_mock_shrunk.csv"), show_col_types = FALSE)
deg$neglog10padj <- -log10(pmax(deg$padj, 1e-300))
deg$sig <- ifelse(!is.na(deg$padj) & deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, "DEG", "NS")

# Volcano
p_vol <- ggplot(deg, aes(log2FoldChange, neglog10padj, color = sig)) +
  geom_point(alpha = 0.6, size = 1.1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(values = c("DEG" = "#d73027", "NS" = "grey70")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Volcano Plot: WT AvrRpm1 vs mock",
    x = "log2FoldChange",
    y = "-log10(adjusted p-value)"
  )

save_png(p_vol, file.path(fig_dir, "ch5_volcano_wt.png"), width = 6.8, height = 5.2)

# MA
p_ma <- ggplot(deg, aes(baseMean, log2FoldChange, color = sig)) +
  geom_point(alpha = 0.6, size = 1.0) +
  scale_x_log10() +
  scale_color_manual(values = c("DEG" = "#d73027", "NS" = "grey70")) +
  theme_bw(base_size = 12) +
  labs(title = "MA Plot: WT AvrRpm1 vs mock", x = "baseMean (log10)", y = "log2FoldChange")

save_png(p_ma, file.path(fig_dir, "ch5_ma_wt.png"), width = 6.8, height = 5.2)
```

---

## 5.5 代码：检查 DEG 数量和结果头部

```bash
# 显著 DEG 数量（WT）
tail -n +2 artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv | wc -l

# 交互项显著 DEG 数量
tail -n +2 artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT_sig.csv | wc -l

# 查看 WT 显著 DEG 前 5 行
head -n 6 artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv

# 查看交互项显著 DEG 前 5 行
head -n 6 artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT_sig.csv
```

### 输出结果

```text
1541
56
"","baseMean","log2FoldChange","lfcSE","pvalue","padj"
"AT5G15950",1058.46429885084,-1.45376652095534,0.132731770663756,5.16733272722286e-29,4.50074680541111e-25
"AT1G09935",24.1197681031312,3.74698979819632,0.364417164046629,2.44079011892465e-26,1.41728546238891e-22
"AT3G17420",88.8579253411927,1.63508829634183,0.162226020945527,4.51659669259174e-25,1.57358228769896e-21
"AT1G09310",2457.14852891811,-1.39998771718618,0.142523808393732,6.84344059187218e-24,1.7030390730059e-20
"AT1G64060",266.611658384728,1.18063546805895,0.122747998195324,7.40067456565989e-23,1.61149688667244e-19
"","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
"AT1G67070",49.4803711208944,-1.5755619595937,0.313552297126098,-5.02487774458903,5.03753339073799e-07,0.00168686843122252
"AT3G61060",127.599719325857,-1.22285834313158,0.253060089891281,-4.83228447305515,1.34975200454455e-06,0.00251098864578772
"AT3G55560",39.2719988293596,-1.19171445089373,0.253284246997521,-4.70504764911572,2.53806668243077e-06,0.00265592815399615
"AT4G38470",206.993199756586,-1.0452717288759,0.221874001988044,-4.71110503939182,2.46377201395389e-06,0.00265592815399615
"AT5G13400",140.757962415624,1.19134842011623,0.250851013492529,4.7492270552525,2.04195583734993e-06,0.00265592815399615
```

---

## 5.6 代码：检查图片文件

```bash
ls docs/assets/validated_case/ch5_*.png
```

### 输出结果

```text
docs/assets/validated_case/ch5_ma_wt.png
docs/assets/validated_case/ch5_volcano_wt.png
```

火山图：

![Volcano WT](assets/validated_case/ch5_volcano_wt.png)

MA 图：

![MA WT](assets/validated_case/ch5_ma_wt.png)

---

## 5.7 本章产物（网页可见）

- `artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_shrunk.csv`
- `artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv`
- `artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT.csv`
- `artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT_sig.csv`
- `docs/assets/validated_case/ch5_volcano_wt.png`
- `docs/assets/validated_case/ch5_ma_wt.png`
