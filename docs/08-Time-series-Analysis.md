# 8. 时间序列分析 (Time-series Analysis)

## 8.1 Motivation

时间序列数据的核心价值在于"轨迹信息"。如果你的实验包含三个或更多时间点，只做 T1 vs T0、T2 vs T0 这种两两比较（pairwise comparison），你会丢失动态过程中的整体模式——一个"先升后降"的基因在 T1 vs T0 和 T2 vs T0 中可能都不显著，但它的时间轨迹整体上是显著变化的。

时间序列分析的核心工具是似然比检验（Likelihood Ratio Test, LRT）。与 Wald 检验不同，LRT 不是比较两个条件的差异，而是比较两个模型的拟合优度：一个包含时间效应的完整模型（full model）和一个去掉时间效应的简化模型（reduced model）。如果完整模型显著优于简化模型，就说明该基因的表达确实随时间（或时间与处理的交互）发生了显著变化。

本章将在 PRJDB11848 的 WT 子集上使用 DESeq2 的 LRT 检验时间相关的动态基因，鉴定出 **2,927 个显著动态基因**（padj < 0.05），并通过聚类热图展示不同的时间响应模式。

---

## 8.2 环境与数据 (Environment and Data)

本章在 R 中完成，需要以下包：

- DESeq2（LRT 检验）
- pheatmap（聚类热图）
- ggplot2（可视化）

输入：计数矩阵 + 包含 `time` 和 `condition` 列的 colData
输出：动态基因列表、时间模式聚类热图

---

## 8.3 什么时候必须用时间序列模型

满足以下任意条件，就不建议只做 pairwise 比较：

1. **时间点 ≥ 3**：两两比较的数量随时间点增加而爆炸，且无法捕捉整体趋势。
2. **你关心非单调变化**：如"先升后降"的脉冲型响应、"延迟启动"的阶梯型响应。
3. **你需要比较不同处理在时间维度上的响应差异**：这是 condition × time 交互效应，pairwise 比较无法直接回答。

---

## 8.4 设计公式与 LRT

LRT 的关键是定义 full model 和 reduced model。两者的差异决定了你在检验什么：

| 生物学问题 | Full model (design) | Reduced model | 检验内容 |
| :--- | :--- | :--- | :--- |
| 单条件下的时间变化 | `~ time` | `~ 1` | 基因是否随时间变化 |
| 两条件共同时间效应 | `~ condition + time` | `~ condition` | 控制处理后，时间是否有效应 |
| 条件-时间交互 | `~ condition + time + condition:time` | `~ condition + time` | 处理的时间响应模式是否不同 |

对于 PRJDB11848 的 WT 子集，我们检验 condition × time 交互效应：

```r
library(DESeq2)

# 设置因子水平
colData$time      <- factor(colData$time, levels = c("0h", "0.5h", "3h"))
colData$condition <- relevel(factor(colData$condition), ref = "mock")

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = colData,
  design    = ~ condition + time + condition:time
)

# LRT：检验交互项是否显著
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ condition + time)
res_lrt <- results(dds_lrt)

sig_dynamic <- subset(as.data.frame(res_lrt), padj < 0.05)
nrow(sig_dynamic)  # PRJDB11848 结果：2927
```

!!! note
    LRT 的 padj < 0.05 表示该基因存在显著的时间交互动态——即它对 AvrRpm1 的响应模式随时间变化。但 LRT 只告诉你"是否动态"，不告诉你"在哪个时间点变化最大"。要定位具体时间点的效应，需要用 Wald 检验提取特定 contrast。

!!! warning "避坑指南"
    `time` 默认建议当作 **factor**（分类变量），除非你有明确的理由假设线性趋势。将时间当作 numeric 会强制拟合线性关系，错过非线性轨迹（如脉冲型响应）。

---

## 8.5 提取具体时间点效应

LRT 筛选出动态基因后，用 Wald 检验提取特定时间点的处理效应：

```r
# 同一设计下运行 Wald 检验
dds_wald <- DESeq(dds)
resultsNames(dds_wald)

# 在参考时间点（0h）下的处理主效应
res_treat_0h <- results(dds_wald, name = "condition_AvrRpm1_vs_mock")

# 处理在 3h 的总效应 = 主效应 + interaction(3h)
res_treat_3h <- results(
  dds_wald,
  list(c("condition_AvrRpm1_vs_mock", "conditionAvrRpm1.time3h"))
)
```

!!! warning
    在含交互项的模型中，"主效应"只对应参考水平（这里是 0h）。`condition_AvrRpm1_vs_mock` 是 AvrRpm1 在 0h 时相对于 mock 的效应，不是所有时间点的平均效应。要获取其他时间点的处理效应，必须将主效应和对应的交互项相加。

---

## 8.6 动态模式可视化

对 LRT 显著的动态基因做聚类热图，可以直观展示不同的时间响应模式：

```r
library(pheatmap)

vsd <- vst(dds_wald, blind = FALSE)
sel_genes <- rownames(sig_dynamic)

mat <- assay(vsd)[sel_genes, ]
mat_z <- t(scale(t(mat)))  # 按基因做 z-score 标准化

anno <- as.data.frame(colData(dds_wald)[, c("time", "condition")])
pheatmap(
  mat_z,
  show_rownames  = FALSE,
  annotation_col = anno,
  cluster_cols   = FALSE,   # 按时间顺序排列，不聚类
  kmeans_k       = 6        # 可选：k-means 分组识别动态模式
)
```

<figure markdown>
  ![LRT padj histogram](assets/validated_case/ch8_lrt_padj_hist.png){ width="580" }
  <figcaption>Fig 8.1 — WT 子集 LRT 检验的 padj 分布直方图。在 0 附近有明显的富集峰，表明存在大量真实的动态基因。padj < 0.05 的基因共 2,927 个。</figcaption>
</figure>

!!! tip "Recommendation"
    先按时间顺序排列样本（不对列聚类），再对基因做 k-means 分组，可以得到"瞬时诱导型"（0.5h 升高、3h 回落）、"持续上调型"（0.5h 和 3h 都升高）、"延迟响应型"（0.5h 不变、3h 才变化）等直观的动态模块。

---

## 8.7 常见翻车点

!!! warning "常见翻车点"
    1. **每个时间点只有 1-2 个重复**：方差估计不稳定，动态结论不可靠。时间序列实验每个时间点至少需要 3 个生物学重复。
    2. **没有 0h 基线**：没有基线时间点，所有的"变化"都缺乏参照，解释非常困难。
    3. **把时间当 numeric 强行拟合线性**：如果真实轨迹是"先升后降"，线性模型会完全错过这个信号。除非有明确的线性假设，否则将时间作为 factor 处理。
    4. **把重复测量当独立样本**：如果同一个体在多个时间点取样（纵向设计），样本之间不独立，需要 mixed model 或配对设计来处理。

---

## 8.8 Key Takeaways

!!! abstract "本章要点"
    1. **时间序列分析的核心是 LRT，不是 pairwise 比较**——LRT 直接检验"基因表达是否随时间变化"，避免了多次两两比较的统计效率损失。
    2. **LRT 筛选"是否动态"，Wald 定位"何时变化"**——两步结合才能完整回答时间序列问题。
    3. **交互模型中的主效应只对应参考时间点**——要获取其他时间点的处理效应，必须将主效应和交互项相加。
    4. **时间应作为 factor 而非 numeric 处理**——除非你有明确的线性假设，否则 factor 编码能捕捉任意形状的时间轨迹。
    5. **动态模式聚类是解读时间序列结果的关键可视化**——k-means 分组可以识别"瞬时诱导""持续上调""延迟响应"等生物学上有意义的模式。
    6. **每个时间点至少 3 个生物学重复**——重复数不足时，时间序列分析的结论不可靠。

---

## 8.9 References

[^love2014]: Love MI, Huber W, Anders S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550. DOI: [10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8)

---

## 8.10 Contributors

### Authors

- *(待补充)*

### Reviewers

- *(待补充)*
