# Bulk RNA-seq 数据分析最佳实践 (Part 2)

## 第二章：数据获取与质量控制 (Data Acquisition & QC)

本章按固定格式组织：先讲原理，再给代码，再给代码输出结果。

---

## 2.1 原理：QC 的目标是判断数据“是否可用于统计分析”

- 下载是否完整（样本数、FASTQ 文件数）
- 运行是否稳定（是否需要重试、是否有失败）
- 定量前质量是否稳定（MultiQC 汇总、映射率分布）

---

## 2.2 代码：检查元数据与下载结果

```bash
# 1) 样本表行数（含表头）
wc -l artifacts/prjdb11848/metadata/samplesheet.csv

# 2) 下载成功的 FASTQ 记录数（download log 中 URL -> 文件 [1]）
rg -n 'fastq.gz" \[1\]$' artifacts/prjdb11848/logs/download_fastq.log | wc -l

# 3) 下载阶段出现的 502 重试次数
rg -n 'ERROR 502' artifacts/prjdb11848/logs/download_fastq.log | wc -l
```

### 输出结果

```text
37 /home/data/t060551/Codex/RNA-seq-Tutorial/artifacts/prjdb11848/metadata/samplesheet.csv
72
3
```

解释：

- `37` 行 = `36` 个样本 + `1` 行表头
- 成功下载 `72` 个 FASTQ（36 对 PE）
- 发生过 `3` 次 `502`，但重试后全部成功

---

## 2.3 代码：检查 MultiQC 汇总和 Salmon 映射率

```bash
# 1) 查看 MultiQC 关键输出
rg -n 'Found 36 meta reports|Report\\s+:' artifacts/prjdb11848/logs/multiqc.log -N

# 2) 汇总 36 个样本的 mapping rate（min/mean/max）
rg -N 'Mapping rate =' artifacts/prjdb11848/logs/salmon_quant_*.log \
  | sed -E 's/.*Mapping rate = ([0-9.]+)%.*/\1/' \
  | awk 'BEGIN{min=100;max=0;sum=0;n=0} {x=$1; if(x<min)min=x; if(x>max)max=x; sum+=x; n++} END{printf "n=%d min=%.3f mean=%.3f max=%.3f\n",n,min,sum/n,max}'
```

### 输出结果

```text
|            salmon | Found 36 meta reports
|           multiqc | Report      : validation_run_downstream/data/multiqc/multiqc_report.html
n=36 min=63.953 mean=95.231 max=97.471
```

解释：

- MultiQC 识别到 36 个 Salmon 报告
- 映射率整体较高（均值 95.231%）
- 最低样本约 63.953%，需要在后续分析中重点关注

---

## 2.4 本章产物（网页可见）

- `artifacts/prjdb11848/metadata/samplesheet.csv`
- `artifacts/prjdb11848/logs/download_fastq.log`
- `artifacts/prjdb11848/logs/multiqc.log`
- `artifacts/prjdb11848/logs/salmon_quant_*.log`

---

## 2.5 小结

本章给出的 QC 结果已经是实跑输出，不是示例文本。下一章进入比对与定量的实现细节。
