# 11. 故障排查与复现保障 (Troubleshooting and Reproducibility)

## 11.1 Motivation

RNA-seq 分析流程涉及多个环节——从数据下载、比对定量到差异分析和功能富集——任何一步出错都可能导致结果偏差甚至完全无效。系统性的故障排查能力是区分"能跑通流程"和"能产出可靠结论"的关键分水岭。在实际项目中，我们反复观察到：一些看似微小的错误（如物种参考基因组不匹配）可以浪费数周的分析时间，而设计矩阵不满秩等问题则会悄无声息地产出错误结果，直到审稿人或合作者质疑时才被发现。

除了排查具体错误，可复现性 (Reproducibility) 是计算生物学的基本要求。一次分析如果无法被他人（或未来的自己）重现，其科学价值将大打折扣。本章汇总了我们在本教程实测过程中遇到的高频问题及其解决方案，并在最后给出一份复现保障的最小清单。我们的目标是：让你在遇到同类问题时能快速定位原因，而不是在 Stack Overflow 上漫无目的地搜索。

## 11.2 高优先级排查：物种与参考是否匹配

!!! warning "陷阱：物种与参考基因组不匹配"

    在本项目实测中，最先遇到的严重问题是：

    - 下载到的 SRR（`SRR8694017-4020`）实际物种是 `Gossypium`（棉花）
    - 却使用了拟南芥 `TAIR10` 参考
    - 直接导致超低比对率与"无 DEG"

    **这类错误不会报错，只会产出"看起来正常但结果全错"的输出。**

建议在下载前先核验 ENA 元数据：

```bash
for r in SRR8694017 SRR8694018 SRR8694019 SRR8694020; do
  wget -qO- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${r}&result=read_run&fields=run_accession,scientific_name,study_accession,sample_accession"
done
```

## 11.3 设计矩阵不满秩（DESeq2）

!!! warning "报错：the model matrix is not full rank"

    报错特征：

    ```text
    the model matrix is not full rank
    ```

    常见原因：

    - `batch` 与 `condition` 完全共线（例如一个批次只含一个条件）。

    处理策略：

    1. 优先修正实验设计（最佳）。
    2. 若当前数据无法分离批次与处理，降阶到可估计模型（如 `~ condition`）。
    3. 在结果报告中显式写明该限制。

## 11.4 `lfcShrink(type="apeglm")` 版本兼容

!!! tip "lfcShrink 版本兼容"

    在新版本 DESeq2 中，`apeglm` 推荐配 `coef`，不是 `contrast`。

```r
coef_name <- grep("^condition_", resultsNames(dds), value = TRUE)[1]
res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
```

## 11.5 GSEA 并行端口受限

!!! tip "GSEA 并行端口受限的解决方案"

    在受限环境中，`BiocParallel` 可能无法开端口，导致 `gseGO()` 失败。
    可强制串行来绕过此问题。

```r
gse <- gseGO(
  geneList = gene_rank_entrez,
  OrgDb = org.At.tair.db,
  ont = "BP",
  keyType = "ENTREZID",
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
```

## 11.6 WGCNA 运行时间过长

!!! tip "WGCNA 运行时间优化"

    大规模基因直接构网可能耗时非常久。
    教学与实战折中策略：

    1. 先过滤低表达基因。
    2. 只保留高变基因子集（例如 5k）。
    3. 在报告中声明筛选规则，保证可审计性。

## 11.7 复现最小清单

!!! abstract "复现清单"

    建议每次分析都固定输出：

    1. `samplesheet.csv`
    2. 软件版本清单
    3. 关键日志（下载、定量、R）
    4. `sessionInfo()`
    5. 中间结果矩阵与最终统计表

    这样才能做到"可复现"而不是"可运行一次"。

## 11.8 Key Takeaways

!!! abstract "本章要点 (Key Takeaways)"

    1. **物种-参考匹配是第零步**：下载数据后、比对之前，务必用 ENA 元数据核验物种信息，物种错误不会报错但会让所有下游分析失效。
    2. **设计矩阵问题要在建模前排查**：`batch` 与 `condition` 共线会导致不满秩，DESeq2 会报错；如果无法修正实验设计，应降阶模型并在报告中声明限制。
    3. **关注工具版本变更**：DESeq2 的 `lfcShrink` 接口在不同版本间有行为差异，始终查阅当前版本文档并使用 `coef` 参数配合 `apeglm`。
    4. **环境受限时回退到串行**：`BiocParallel` 在某些计算环境中无法正常并行，使用 `SerialParam()` 是最稳妥的兜底方案。
    5. **可复现 ≠ 可运行一次**：每次分析都应输出样本表、软件版本、关键日志、`sessionInfo()` 和中间结果，这是科学计算的基本纪律。
