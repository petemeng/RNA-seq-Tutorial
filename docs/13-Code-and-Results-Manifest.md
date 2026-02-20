# Chapter 13: 代码与结果索引（Runbook -> Scripts -> Artifacts）

本页用于把“步骤、代码、结果文件”一一对齐，确保网页版可直接审计。

---

## 13.1 已上传产物快照

- 产物根目录（GitHub）：[artifacts/prjdb11848](https://github.com/petemeng/RNA-seq-Tutorial/tree/main/artifacts/prjdb11848)
- 文件清单：[FILELIST.txt](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/FILELIST.txt)
- 完整校验：[CHECKSUMS.sha256](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/CHECKSUMS.sha256)
- 版本信息：[tool_versions.txt](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/tool_versions.txt)

当前快照规模：

- 结果文件：20 个
- 日志文件：39 个
- 总体积：约 15 MB

---

## 13.2 分步映射（scRNA-best-practice 风格）

| Step | 命令 | 代码 | 运行目录产物 | 已上传产物 |
| :--- | :--- | :--- | :--- | :--- |
| 1 | `bash scripts/01_prepare_prjdb11848_samplesheet.sh` | `scripts/01_prepare_prjdb11848_samplesheet.sh` | `validation_run_downstream/metadata/` | `artifacts/prjdb11848/metadata/` |
| 2 | `bash scripts/02_download_fastq.sh` | `scripts/02_download_fastq.sh` | `validation_run_downstream/logs/download_fastq.log` | `artifacts/prjdb11848/logs/download_fastq.log` |
| 3 | `bash scripts/03_quantify_salmon.sh` | `scripts/03_quantify_salmon.sh` | `validation_run_downstream/logs/salmon_quant_*.log` | `artifacts/prjdb11848/logs/salmon_quant_*.log` |
| 4 | `Rscript scripts/04_downstream_ch4_to_ch9.R` | `scripts/04_downstream_ch4_to_ch9.R` | `validation_run_downstream/results/ch4..ch9/` | `artifacts/prjdb11848/results/ch4..ch9/` |
| 5 | `Rscript scripts/05_generate_case_figures.R` | `scripts/05_generate_case_figures.R` | `docs/assets/validated_case/*.png` | `docs/assets/validated_case/*.png` |

---

## 13.3 关键结果文件（直接看网页可下）

差异表达：

- [DEG_WT_AvrRpm1_vs_mock_sig.csv](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv)
- [DEG_interaction_clf_vs_WT_sig.csv](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/results/ch5/DEG_interaction_clf_vs_WT_sig.csv)

富集分析：

- [GO_ORA_WT.csv](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/results/ch6/GO_ORA_WT.csv)
- [GO_GSEA_WT.csv](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/results/ch6/GO_GSEA_WT.csv)

时间序列：

- [time_series_LRT_WT_sig.csv](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/results/ch8/time_series_LRT_WT_sig.csv)

多因素交互：

- [res_treat_in_wt.csv](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/results/ch9/res_treat_in_wt.csv)
- [res_interaction.csv](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/results/ch9/res_interaction.csv)
- [res_treat_in_clf.csv](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/results/ch9/res_treat_in_clf.csv)

---

## 13.4 快速验收命令

```bash
# WT 显著 DEG 数量
tail -n +2 artifacts/prjdb11848/results/ch5/DEG_WT_AvrRpm1_vs_mock_sig.csv | wc -l
# 预期: 1541

# WT 时间序列 LRT 显著基因数量
tail -n +2 artifacts/prjdb11848/results/ch8/time_series_LRT_WT_sig.csv | wc -l
# 预期: 2927

# 结果清单与校验
wc -l artifacts/prjdb11848/FILELIST.txt
sha256sum -c artifacts/prjdb11848/CHECKSUMS.sha256
```
