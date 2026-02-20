# 10. 可复现实战教程：PRJDB11848 (Validated Case Study)

## 10.1 Motivation

教程章节（Ch0-9）逐步讲解了 RNA-seq 分析的每个环节，但读者最终需要回答一个问题：把这些步骤串起来，在真实数据上能跑通吗？结果和教程里展示的一致吗？端到端复现（end-to-end reproduction）是检验教程可信度的唯一标准——如果读者按照说明操作却得到不同的数字，那教程的价值就大打折扣。

本章提供一个"能直接照做"的完整流程：从下载原始 FASTQ 到生成差异基因列表和富集分析结果，全部使用公开数据集 PRJDB11848（拟南芥 AvrRpm1 免疫响应）。每一步都给出实际运行的命令、预期产物和真实结果图。你可以将自己的运行结果与本章的验收数字逐一比对，确认复现成功。

与教学章节不同，本章不再解释"为什么这么做"——原理请回到对应章节查阅。这里只关注"怎么跑"和"跑出来应该是什么样"。

!!! note "本章定位"
    本章是"能直接照做"的版本：给出已实跑成功的代码、每一步会产出的文件、以及真实结果图。原理讲解请参见 Ch0-9 对应章节。

## 10.2 代码与结果的网页可见位置

本章同时区分两类路径：

- 运行目录（本地执行时会写入）：`validation_run_downstream/`
- 已上传快照（GitHub 网页可直接查看）：`artifacts/prjdb11848/`

快速入口：

- 产物目录：`artifacts/prjdb11848/`
- 文件清单：`artifacts/prjdb11848/FILELIST.txt`
- 校验文件：`artifacts/prjdb11848/CHECKSUMS.sha256`

## 10.3 实战数据与设计

- 数据集：`PRJDB11848`（`Arabidopsis thaliana`）
- 样本数：36（PE，72 个 FASTQ）
- 设计因子：
  - `genotype`: `WT` / `clf`
  - `condition`: `mock` / `AvrRpm1`
  - `time`: `0h` / `0.5h` / `3h`
  - `replicate`: `1..3`

## 10.4 运行前准备

在仓库根目录：

```bash
cd RNA-seq-Tutorial
```

确保工具可用：

- `wget`
- `salmon`
- `Rscript`
- R 包：`DESeq2`, `tximport`, `clusterProfiler`, `org.At.tair.db`, `WGCNA`, `pheatmap`, `BiocParallel`

!!! warning "磁盘空间与路径依赖"
    完整运行需要约 20 GB 磁盘空间（72 个 FASTQ 文件 + 定量结果 + 下游分析产物）。所有脚本假设从仓库根目录 `RNA-seq-Tutorial/` 执行，请勿在子目录中运行。

## 10.5 Step 1: 生成标准样本表

```bash
bash scripts/01_prepare_prjdb11848_samplesheet.sh
```

产物：

- `validation_run_downstream/metadata/prjdb11848.tsv`
- `validation_run_downstream/metadata/samplesheet.csv`

检查点：

- `samplesheet.csv` 行数应为 `37`（含表头）
- 字段包含 `sample_id,run,genotype,condition,time,replicate,layout,fastq_1,fastq_2,batch`

## 10.6 Step 2: 下载 FASTQ

```bash
bash scripts/02_download_fastq.sh
```

产物：

- `validation_run_downstream/data/raw_data/*.fastq.gz`（72 个）
- `validation_run_downstream/logs/download_fastq.log`

检查点：

```bash
ls validation_run_downstream/data/raw_data/*.fastq.gz | wc -l
```

预期输出：`72`

## 10.7 Step 3: Salmon 定量

```bash
bash scripts/03_quantify_salmon.sh
```

产物：

- `validation_run_downstream/data/quant/<sample_id>/quant.sf`
- `validation_run_downstream/logs/salmon_quant_<sample_id>.log`

检查点：

```bash
ls validation_run_downstream/data/quant | wc -l
```

预期输出：`36`

本次实跑 mapping rate 统计（36 样本）：

- `min = 63.953%`
- `mean = 95.231%`
- `max = 97.471%`

## 10.8 Step 4: Chapter 4-9 一体化下游分析

```bash
Rscript scripts/04_downstream_ch4_to_ch9.R
```

产物目录：

```text
validation_run_downstream/results/
├── ch4/
├── ch5/
├── ch6/
├── ch7/
├── ch8/
├── ch9/
└── R_sessionInfo.txt
```

## 10.9 Step 5: 生成教程展示图

```bash
Rscript scripts/05_generate_case_figures.R
```

产物：

- `docs/assets/validated_case/ch4_pca.png`
- `docs/assets/validated_case/ch5_volcano_wt.png`
- `docs/assets/validated_case/ch5_ma_wt.png`
- `docs/assets/validated_case/ch6_go_ora_top15.png`
- `docs/assets/validated_case/ch6_go_gsea_top15.png`
- `docs/assets/validated_case/ch8_lrt_padj_hist.png`
- `docs/assets/validated_case/ch9_interaction_volcano.png`

## 10.10 真实结果展示（本仓库实跑）

### Chapter 4: PCA

<figure markdown>
![PCA](assets/validated_case/ch4_pca.png)
<figcaption>Fig. 10.1 PCA 双标图。样本按 condition（颜色）和 genotype（形状）分组，PC1 主要分离 AvrRpm1 与 mock 处理。</figcaption>
</figure>

关键文件：

- `validation_run_downstream/results/ch4/PCA_data.csv`
- `validation_run_downstream/results/ch4/sample_distance_heatmap.pdf`

### Chapter 5: DEG（WT: AvrRpm1 vs mock）

统计结果：

- 显著 DEG（`padj < 0.05` 且 `|log2FC| > 1`）：`1541`
- 交互项显著 DEG（`genotype:condition`）：`56`

<figure markdown>
![Volcano WT](assets/validated_case/ch5_volcano_wt.png)
<figcaption>Fig. 10.2 火山图（WT: AvrRpm1 vs mock）。红色点为显著 DEG（padj < 0.05 且 |log2FC| > 1），共 1541 个基因。</figcaption>
</figure>

<figure markdown>
![MA WT](assets/validated_case/ch5_ma_wt.png)
<figcaption>Fig. 10.3 MA 图（WT: AvrRpm1 vs mock）。横轴为 baseMean（log10 标度），纵轴为 log2FoldChange，红色点为显著 DEG。</figcaption>
</figure>

结果文件：

- `validation_run_downstream/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv`
- `validation_run_downstream/results/ch5/DEG_interaction_clf_vs_WT_sig.csv`

### Chapter 6: 富集分析

<figure markdown>
![GO ORA](assets/validated_case/ch6_go_ora_top15.png)
<figcaption>Fig. 10.4 GO ORA Top 15 富集条目。横轴为 -log10(adjusted p-value)，展示 WT DEG 中最显著富集的生物学过程。</figcaption>
</figure>

<figure markdown>
![GO GSEA](assets/validated_case/ch6_go_gsea_top15.png)
<figcaption>Fig. 10.5 GO GSEA Top 15 富集条目。横轴为 NES（Normalized Enrichment Score），红色为上调通路，蓝色为下调通路。</figcaption>
</figure>

结果文件：

- `validation_run_downstream/results/ch6/GO_ORA_WT.csv`
- `validation_run_downstream/results/ch6/GO_GSEA_WT.csv`

### Chapter 7: WGCNA

输出文件：

- `validation_run_downstream/results/ch7/wgcna_soft_threshold.csv`
- `validation_run_downstream/results/ch7/wgcna_modules.csv`
- `validation_run_downstream/results/ch7/module_trait_cor.csv`
- `validation_run_downstream/results/ch7/module_trait_p.csv`

!!! note "运行时长控制"
    为保证运行时长可控，实跑时在脚本中限制为高变 5000 基因构网。

### Chapter 8: 时间序列（WT 子集 LRT）

统计结果：

- 显著动态基因：`2927`

<figure markdown>
![LRT padj](assets/validated_case/ch8_lrt_padj_hist.png)
<figcaption>Fig. 10.6 LRT adjusted p-value 分布直方图（WT 子集）。左侧峰值表明大量基因的表达随时间发生显著动态变化。</figcaption>
</figure>

结果文件：

- `validation_run_downstream/results/ch8/time_series_LRT_WT.csv`
- `validation_run_downstream/results/ch8/time_series_LRT_WT_sig.csv`

### Chapter 9: 多因素交互

<figure markdown>
![Interaction Volcano](assets/validated_case/ch9_interaction_volcano.png)
<figcaption>Fig. 10.7 交互项火山图（clf:AvrRpm1）。紫色点为交互效应显著的基因，即 clf 突变体对 AvrRpm1 处理的响应与 WT 存在显著差异的基因。</figcaption>
</figure>

结果文件：

- `validation_run_downstream/results/ch9/res_treat_in_wt.csv`
- `validation_run_downstream/results/ch9/res_interaction.csv`
- `validation_run_downstream/results/ch9/res_treat_in_clf.csv`

## 10.11 本章对应脚本（完整代码）

- `scripts/01_prepare_prjdb11848_samplesheet.sh`
- `scripts/02_download_fastq.sh`
- `scripts/03_quantify_salmon.sh`
- `scripts/04_downstream_ch4_to_ch9.R`
- `scripts/05_generate_case_figures.R`

脚本全文请见第 12 章：`docs/12-Validated-Case-Study-Full-Scripts.md`。

## 10.12 复现验收清单

你可以按下面最小清单确认"完整复现成功"：

1. `quant` 目录样本数为 `36`
2. `DEG_WT_AvrRpm1_vs_mock_sig.csv` 行数（去表头）为 `1541`
3. `time_series_LRT_WT_sig.csv` 行数（去表头）为 `2927`
4. `docs/assets/validated_case/` 下有 7 张结果图
5. `validation_run_downstream/results/R_sessionInfo.txt` 存在

!!! abstract "本章要点 (Key Takeaways)"
    1. 端到端复现是检验教程可信度的金标准——如果读者跑出不同的数字，教程就需要修正。
    2. 所有脚本从仓库根目录执行，由 `samplesheet.csv` 驱动，避免硬编码路径和样本名。
    3. 关键验收数字（DEG 1541、LRT 2927、样本 36、图 7 张）是快速判断复现是否成功的锚点。
    4. `artifacts/prjdb11848/` 提供已上传的结果快照，即使不运行也可以在 GitHub 网页上直接查看和下载。
    5. 复现不等于"代码能跑"——还需要比对输出数字、检查图表模式、验证 `sessionInfo()` 中的软件版本。

## 10.13 References

[^sakata2023]: Sakata N, Stonoha-Arber C, Wang D, et al. *Plant Cell Physiol.* 2023. PRJDB11848 数据集：拟南芥 WT 与 clf 突变体对 AvrRpm1 的免疫响应时间序列转录组。
