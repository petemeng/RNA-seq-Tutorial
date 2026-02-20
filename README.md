# RNA-seq Tutorial (sc-best-practices Style)

一个以"可运行 + 可验收 + 可追溯"为核心的 Bulk RNA-seq 教程仓库。

在线文档：<https://petemeng.github.io/RNA-seq-Tutorial/>

## 5 步跑通（PRJDB11848）

```bash
cd RNA-seq-Tutorial

# 1) 生成样本表
bash scripts/01_prepare_prjdb11848_samplesheet.sh

# 2) 下载 FASTQ
bash scripts/02_download_fastq.sh

# 3) Salmon 定量
bash scripts/03_quantify_salmon.sh

# 4) Chapter 4-9 下游分析
Rscript scripts/04_downstream_ch4_to_ch9.R

# 5) 生成教程展示图
Rscript scripts/05_generate_case_figures.R
```

## 验收数字

| 指标 | 预期值 |
|------|--------|
| WT 显著 DEG（padj < 0.05, \|log2FC\| > 1） | 1541 |
| 时间序列 LRT 显著基因 | 2927 |
| Salmon 定量样本数 | 36 |
| 教程展示图 | 7 张 |

## 产物快照

已上传至 `artifacts/prjdb11848/`，包含结果文件 20 个、日志 39 个，可通过 `CHECKSUMS.sha256` 校验完整性。

## 本地预览文档

```bash
pip install mkdocs-material
mkdocs serve
```

详细教程内容请访问[在线文档](https://petemeng.github.io/RNA-seq-Tutorial/)或查看 `docs/index.md`。
