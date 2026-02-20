# Bulk RNA-seq 数据分析最佳实践 (Part 2)

## 第二章：数据获取与质量控制 (Data Acquisition & QC)

> 核心原则：QC 不是“走流程”，而是决定数据是否可用于统计推断。

---

## 2.1 从样本表驱动下载，而不是手敲命令

假设你已经有 `metadata/samplesheet.csv`，并且其中 `fastq_1/fastq_2` 是可下载 URL。

目录初始化：

```bash
mkdir -p data/{raw_data,clean_data,qc/fastqc_raw,qc/fastqc_clean,fastp_reports,multiqc}
```

示例批量下载脚本：

```bash
#!/usr/bin/env bash
set -euo pipefail

awk -F',' 'NR>1 {print $6"\n"$7}' metadata/samplesheet.csv | sed '/^$/d' | while read -r url; do
  echo "Downloading: ${url}"
  wget -c "$url" -P data/raw_data/
done
```

> 避坑指南：下载脚本请保留日志（`tee logs/download.log`），后续追溯非常有用。

---

## 2.2 FastQC：先判断数据是否“可救”

对原始 FASTQ 运行 FastQC：

```bash
fastqc -t 8 -o data/qc/fastqc_raw data/raw_data/*.fastq.gz
```

重点观察四类指标：

1. `Per base sequence quality`：末端轻微下降常见，但大面积 <Q20 需要关注。
2. `Adapter content`：明显上升说明需要去接头。
3. `Per sequence GC content`：异常双峰可能提示污染或混样。
4. `Overrepresented sequences`：判断接头、rRNA 或外源污染。

---

## 2.3 fastp 清洗：默认参数优先，避免过度裁剪

### 单样本示例

```bash
fastp \
  -i data/raw_data/SRR8694017_1.fastq.gz \
  -I data/raw_data/SRR8694017_2.fastq.gz \
  -o data/clean_data/SRR8694017_1.clean.fastq.gz \
  -O data/clean_data/SRR8694017_2.clean.fastq.gz \
  --thread 8 \
  --detect_adapter_for_pe \
  -h data/fastp_reports/SRR8694017.fastp.html \
  -j data/fastp_reports/SRR8694017.fastp.json
```

### 多样本批处理（samplesheet 驱动）

```bash
#!/usr/bin/env bash
set -euo pipefail

while IFS=',' read -r sample_id run condition replicate layout fastq1 fastq2 batch time; do
  [ "$sample_id" = "sample_id" ] && continue

  f1="data/raw_data/$(basename "$fastq1")"
  f2="data/raw_data/$(basename "$fastq2")"

  fastp \
    -i "$f1" -I "$f2" \
    -o "data/clean_data/${sample_id}_1.clean.fastq.gz" \
    -O "data/clean_data/${sample_id}_2.clean.fastq.gz" \
    --thread 8 \
    --detect_adapter_for_pe \
    -h "data/fastp_reports/${sample_id}.html" \
    -j "data/fastp_reports/${sample_id}.json"
done < metadata/samplesheet.csv
```

> Best Practice：先用默认参数跑一轮，再根据报告调整，不要一上来就激进过滤。

---

## 2.4 清洗后复检 + MultiQC 汇总

```bash
fastqc -t 8 -o data/qc/fastqc_clean data/clean_data/*.fastq.gz

multiqc data \
  --outdir data/multiqc \
  --filename multiqc_report.html
```

MultiQC 是你给 PI 和同事汇报质量的主入口。

---

## 2.5 如何判定“通过 QC”

没有绝对阈值，但可以用以下经验标准做初筛：

| 指标 | 建议范围 | 说明 |
| :--- | :--- | :--- |
| Q30 比例 | > 80%（常见） | 具体与平台相关 |
| Adapter 含量 | 清洗后显著下降 | 不追求绝对 0 |
| 序列长度分布 | 与建库预期一致 | 异常截短需排查 |
| 样本间总 reads | 量级接近 | 极端低深度样本单独评估 |

> 避坑指南：QC 的目标是“保证可分析性”，不是把图修到全绿。

---

## 2.6 常见问题与应对

1. **某样本 reads 数极低**：优先排查下载完整性和建库失败，不要盲目继续。
2. **Adapter 清不干净**：确认接头来源，必要时手动指定 `--adapter_sequence`。
3. **GC 曲线异常**：考虑污染、混样或物种不匹配。
4. **重复率偏高**：在高表达组织中可能正常，需结合文库复杂度解释。

---

## 2.7 本章检查清单

- 原始与清洗后 FastQC 都完成。
- MultiQC 报告已生成并保存。
- 低质量或异常样本有处理记录（剔除或保留理由）。
- 清洗后的 FASTQ 文件命名规则一致。

完成后进入比对/定量章节。

---

## 2.8 本章实跑代码与结果（PRJDB11848）

本章在实跑流程中对应的脚本：

```bash
bash scripts/02_download_fastq.sh
```

已上传快照：

- `artifacts/prjdb11848/logs/download_fastq.log`
- `artifacts/prjdb11848/logs/multiqc.log`
- `artifacts/prjdb11848/CHECKSUMS.sha256`

快速验收：

```bash
ls artifacts/prjdb11848/logs/salmon_quant_*.log | wc -l
# 预期: 36
```

说明：本次验证 run 以 Salmon 为主流程，未额外执行 fastp 输出快照；如需裁剪，可按本章命令在本地补跑。

