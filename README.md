# RNA-seq Tutorial (scRNA-best-practice Style)

一个以“可运行 + 可验收 + 可追溯”为核心的 Bulk RNA-seq 教程仓库。

## 你现在可以在网页上直接看到什么

- 分步运行代码：`scripts/`
- 本次实跑结果快照：`artifacts/prjdb11848/`
- 教程图像结果：`docs/assets/validated_case/`
- 端到端实跑文档：`docs/10-Validated-Case-Study-PRJDB11848.md`
- 完整脚本源码：`docs/12-Validated-Case-Study-Full-Scripts.md`

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

## 结果如何验收

- DEG 数量（WT AvrRpm1 vs mock）：
```bash
tail -n +2 validation_run_downstream/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv | wc -l
```
预期：`1541`

- 时间序列 LRT 显著基因数：
```bash
tail -n +2 validation_run_downstream/results/ch8/time_series_LRT_WT_sig.csv | wc -l
```
预期：`2927`

## 已上传的“本次运行”产物目录

- 产物根目录：`artifacts/prjdb11848/`
- 文件清单：`artifacts/prjdb11848/FILELIST.txt`
- SHA256 校验：`artifacts/prjdb11848/CHECKSUMS.sha256`
- 环境版本：`artifacts/prjdb11848/tool_versions.txt`

## 文档预览

```bash
pip install mkdocs-material
mkdocs serve
```
