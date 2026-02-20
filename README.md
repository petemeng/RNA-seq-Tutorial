# RNA-seq Data Analysis Best Practices

一套可复现、可审计、可发表的 Bulk RNA-seq 教程。

## 本仓库现在包含
- 完整教程文档（Chapter 0-9）
- 实测跑通案例（拟南芥 `PRJDB11848`，36 样本）
- 可直接执行的脚本（下载、定量、Chapter 4-9 下游）
- 真实结果展示（DEG、富集、WGCNA、时间序列、多因素交互）

## 关键更新（2026-02-20）
- 将实跑代码与结果整合进原有章节（`docs/03` 到 `docs/09`）
- 新增可运行脚本：`scripts/`

## 快速开始
```bash
cd RNA-seq-Tutorial

# 1) 准备样本表（PRJDB11848）
bash scripts/01_prepare_prjdb11848_samplesheet.sh

# 2) 下载 FASTQ
bash scripts/02_download_fastq.sh

# 3) Salmon 定量
bash scripts/03_quantify_salmon.sh

# 4) 下游 4-9 章
Rscript scripts/04_downstream_ch4_to_ch9.R

# 5) 生成教程展示图（火山图等）
Rscript scripts/05_generate_case_figures.R
```

## 文档预览
```bash
pip install mkdocs-material
mkdocs serve
```

## 关键文档入口
- 基础流程与实跑内容（定量到富集）：`docs/03-Alignment-and-Quantification.md` 到 `docs/06-Functional-Enrichment-Analysis.md`
- 进阶流程与实跑内容（WGCNA/时间序列/交互）：`docs/07-WGCNA.md` 到 `docs/09-Multi-factor-Design.md`
