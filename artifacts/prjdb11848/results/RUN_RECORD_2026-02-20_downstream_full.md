# 下游全流程运行记录（Chapter 4-9）

## 背景修正
- 原先 `SRR8694017-4020` 经 ENA 核验为棉花（Gossypium），与拟南芥参考不匹配，导致低映射和无 DEG。
- 本次改用有效拟南芥数据集：`PRJDB11848`（36 样本；genotype=WT/clf，condition=mock/AvrRpm1，time=0h/0.5h/3h，rep=1..3）。

## 已执行内容
1. 下载 36 样本（72 FASTQ）到 `validation_run_downstream/data/raw_data`。
2. Salmon 定量（36/36）：`validation_run_downstream/data/quant/*/quant.sf`。
3. MultiQC 汇总：`validation_run_downstream/data/multiqc/multiqc_report.html`。
4. 运行第 4-9 章一体化脚本：`validation_run_downstream/scripts/run_ch4_to_ch9.R`。

## 关键结果
- Chapter 4（标准化+EDA）
  - `validation_run_downstream/results/ch4/PCA_data.csv`
  - `validation_run_downstream/results/ch4/sample_distance_heatmap.pdf`
- Chapter 5（DEA）
  - WT: AvrRpm1 vs mock 显著 DEG：`1541`（文件含表头共 1542 行）
  - 交互项 clf:AvrRpm1 显著 DEG：`56`（文件含表头共 57 行）
- Chapter 6（富集）
  - GO ORA：`validation_run_downstream/results/ch6/GO_ORA_WT.csv`
  - GO GSEA：`validation_run_downstream/results/ch6/GO_GSEA_WT.csv`
- Chapter 7（WGCNA）
  - `wgcna_soft_threshold.csv`
  - `wgcna_modules.csv`
  - `module_trait_cor.csv`
  - `module_trait_p.csv`
- Chapter 8（时间序列 LRT, WT 子集）
  - 显著动态基因：`2927`（文件含表头共 2928 行）
- Chapter 9（多因素交互）
  - `res_treat_in_wt.csv`
  - `res_interaction.csv`
  - `res_treat_in_clf.csv`

## 运行中做的兼容修复
- GSEA 在该环境下因 BiocParallel 端口限制失败，改为 `SerialParam()` 后通过。
- WGCNA 为保证时效，使用高变 5000 基因子集（其余流程不变）。

## 质量信号
- Salmon 映射率约 96%-97%（示例见 `validation_run_downstream/logs/salmon_quant_*.log`），与拟南芥参考匹配良好。
