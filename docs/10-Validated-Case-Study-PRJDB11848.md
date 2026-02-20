# Chapter 10: 实测完整案例（PRJDB11848，拟南芥）

> 目标：给出一套“从下载到 Chapter 9 下游结果”的可直接复现实战流程。

---

## 10.1 为什么新增这一章

教程示例必须经真实数据验证。  
本章使用 ENA 项目 `PRJDB11848`（`Arabidopsis thaliana`，36 样本）实跑并固化脚本。

设计结构：

- `genotype`: `WT` / `clf`
- `condition`: `mock` / `AvrRpm1`
- `time`: `0h` / `0.5h` / `3h`
- `replicate`: `1..3`

---

## 10.2 一键式执行顺序

在仓库根目录运行：

```bash
# 1) 生成样本表
bash scripts/01_prepare_prjdb11848_samplesheet.sh

# 2) 下载 72 个 FASTQ（36 样本 x PE）
bash scripts/02_download_fastq.sh

# 3) Salmon 定量
bash scripts/03_quantify_salmon.sh

# 4) Chapter 4-9 下游分析
Rscript scripts/04_downstream_ch4_to_ch9.R
```

> 默认工作目录为 `validation_run_downstream/`，可通过参数或环境变量调整。

---

## 10.3 实测运行结果（2026-02-20）

### 定量质量

- Salmon 成功样本数：`36/36`
- Mapping rate 统计（来自 `validation_run_downstream/logs/salmon_quant_*.log`）：
  - `min = 63.953%`
  - `mean = 95.231%`
  - `max = 97.471%`

> 说明：少数样本映射率偏低但整体可分析；建议结合 MultiQC 再做样本级排查。

### Chapter 5: DEG

- WT 中 `AvrRpm1 vs mock` 显著 DEG（`padj < 0.05 & |log2FC| > 1`）：`1541`
- 交互项 `genotype:condition` 显著 DEG：`56`

Top DEG（按 padj）示例：

| Gene | log2FC | padj |
| :--- | :---: | :---: |
| AT5G15950 | -1.454 | 4.50e-25 |
| AT1G09935 | 3.747 | 1.42e-22 |
| AT3G17420 | 1.635 | 1.57e-21 |
| AT1G09310 | -1.400 | 1.70e-20 |
| AT1G64060 | 1.181 | 1.61e-19 |

结果文件：

- `validation_run_downstream/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv`
- `validation_run_downstream/results/ch5/DEG_interaction_clf_vs_WT_sig.csv`

### Chapter 6: 富集

GO ORA top 条目示例：

- `secondary metabolic process`
- `photosynthesis`
- `photosynthesis, light reaction`
- `S-glycoside biosynthetic process`

GO GSEA top 条目示例：

- `photosynthesis, light reaction`（NES 负向）
- `photosynthesis`（NES 负向）
- `cellular response to hypoxia`（NES 正向）

结果文件：

- `validation_run_downstream/results/ch6/GO_ORA_WT.csv`
- `validation_run_downstream/results/ch6/GO_GSEA_WT.csv`

### Chapter 7: WGCNA

已输出：

- `wgcna_soft_threshold.csv`
- `wgcna_modules.csv`
- `module_trait_cor.csv`
- `module_trait_p.csv`

> 为控制运行时长，实跑使用高变 `5000` 基因子集构网（脚本内已注明）。

### Chapter 8: 时间序列

- WT 子集 LRT 显著动态基因：`2927`
- 结果：`validation_run_downstream/results/ch8/time_series_LRT_WT_sig.csv`

### Chapter 9: 多因素交互

已输出三类核心结果：

- `res_treat_in_wt.csv`
- `res_interaction.csv`
- `res_treat_in_clf.csv`

---

## 10.4 结果文件总览

关键产物目录：

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

运行记录：

- `validation_run_downstream/results/RUN_RECORD_2026-02-20_downstream_full.md`

---

## 10.5 教学建议

如果你要写论文方法部分，建议直接引用本章脚本与输出路径，并在补充材料中附上：

1. `samplesheet.csv`
2. 关键日志（下载、quant、R）
3. `R_sessionInfo.txt`
