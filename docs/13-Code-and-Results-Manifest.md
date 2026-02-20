# 13. 代码与结果清单 (Code and Results Manifest)

## 13.1 Motivation

可重复性 (Reproducibility) 是计算生物学的基石。仅仅提供代码和最终结论是不够的——我们需要一份完整的审计轨迹 (Audit Trail)，将每一个分析步骤与其对应的脚本、输入数据和输出产物一一对齐。这样做的好处是双重的：首先，任何人（包括未来的你）都可以沿着清单逐步验证结果是否可复现；其次，同行评审者可以精确定位某个结论来自哪段代码、哪个中间文件，从而对计算结果进行实质性的审查而非仅凭信任。本章就是这样一份清单——它把前面所有章节的"步骤、代码、结果文件"串联起来，确保网页版可直接审计。

## 13.2 已上传产物快照

- 产物根目录（GitHub）：[artifacts/prjdb11848](https://github.com/petemeng/RNA-seq-Tutorial/tree/main/artifacts/prjdb11848)
- 文件清单：[FILELIST.txt](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/FILELIST.txt)
- 完整校验：[CHECKSUMS.sha256](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/CHECKSUMS.sha256)
- 版本信息：[tool_versions.txt](https://github.com/petemeng/RNA-seq-Tutorial/blob/main/artifacts/prjdb11848/tool_versions.txt)

当前快照规模：

- 结果文件：20 个
- 日志文件：39 个
- 总体积：约 15 MB

## 13.3 分步映射

| Step | 命令 | 代码 | 运行目录产物 | 已上传产物 |
| :--- | :--- | :--- | :--- | :--- |
| 1 | `bash scripts/01_prepare_prjdb11848_samplesheet.sh` | `scripts/01_prepare_prjdb11848_samplesheet.sh` | `validation_run_downstream/metadata/` | `artifacts/prjdb11848/metadata/` |
| 2 | `bash scripts/02_download_fastq.sh` | `scripts/02_download_fastq.sh` | `validation_run_downstream/logs/download_fastq.log` | `artifacts/prjdb11848/logs/download_fastq.log` |
| 3 | `bash scripts/03_quantify_salmon.sh` | `scripts/03_quantify_salmon.sh` | `validation_run_downstream/logs/salmon_quant_*.log` | `artifacts/prjdb11848/logs/salmon_quant_*.log` |
| 4 | `Rscript scripts/04_downstream_ch4_to_ch9.R` | `scripts/04_downstream_ch4_to_ch9.R` | `validation_run_downstream/results/ch4..ch9/` | `artifacts/prjdb11848/results/ch4..ch9/` |
| 5 | `Rscript scripts/05_generate_case_figures.R` | `scripts/05_generate_case_figures.R` | `docs/assets/validated_case/*.png` | `docs/assets/validated_case/*.png` |

## 13.4 关键结果文件

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

## 13.5 快速验收命令

!!! tip "验证方法"

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

## 13.6 Key Takeaways

!!! abstract "本章要点 (Key Takeaways)"

    1. **清单即审计轨迹**：每个分析步骤都应可追溯到具体脚本和输出文件。
    2. **SHA256 校验**确保结果文件未被意外修改。
    3. **关键数字**（DEG 1541、LRT 2927）是快速判断复现是否成功的锚点。
