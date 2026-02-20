# 2. 质量控制与数据清洗 (Quality Control and Preprocessing)

## 2.1 Motivation

质量控制（Quality Control, QC）不是"走流程"，而是决定数据是否可用于统计推断的关键判断步骤。一个接头没去干净的样本会在比对时产生大量 soft-clipping，降低有效比对率；一个 GC 含量异常的样本可能存在外源污染，如果不在此阶段发现，它会在 PCA 图上表现为离群点，而你可能会错误地将其归因于生物学差异。

QC 的核心目标不是"把所有指标修到全绿"，而是回答两个问题：第一，这份数据的质量是否足以支撑下游的统计推断？第二，是否存在需要特殊处理的样本（如极低深度、严重污染）？过度清洗和不清洗一样危险——激进的质量过滤会丢失真实的生物学信号，尤其是低表达基因的 reads。

本章将使用 FastQC 评估原始数据质量，用 fastp 进行接头去除和质量过滤，用 MultiQC 生成汇总报告，并建立"通过 QC"的判断标准。

---

## 2.2 环境与数据 (Environment and Data)

本章需要以下工具（均已包含在 [第 1 章](01-Environment-and-References.md) 的 `rnaseq_env` 环境中）：

- FastQC（质量评估）
- fastp（接头去除 + 质量过滤）
- MultiQC（多样本报告汇总）

输入文件：`data/raw_data/*.fastq.gz`（原始测序数据）
输出文件：`data/clean_data/*.fastq.gz`（清洗后数据）、QC 报告

准备目录结构：

```bash
mkdir -p data/{raw_data,clean_data,qc/fastqc_raw,qc/fastqc_clean,fastp_reports,multiqc}
```

---

## 2.3 FastQC：先判断数据是否"可救"

对原始 FASTQ 运行 FastQC，获取每个样本的质量概览：

```bash
fastqc -t 8 -o data/qc/fastqc_raw data/raw_data/*.fastq.gz
```

重点观察四类指标：

1. **Per base sequence quality**：末端轻微下降是正常的（Illumina 测序特性），但大面积低于 Q20 需要关注。
2. **Adapter content**：曲线明显上升说明存在接头残留，需要去接头处理。
3. **Per sequence GC content**：异常双峰可能提示外源污染（如细菌 DNA 混入）或混样错误。
4. **Overrepresented sequences**：高频序列可能是接头、rRNA 或外源污染的信号。

!!! note
    FastQC 的"PASS/WARN/FAIL"标签是基于通用标准设定的，不一定适用于所有实验类型。例如，RNA-seq 数据在"Per sequence GC content"和"Sequence Duplication Levels"上经常显示 WARN 或 FAIL，这通常是正常的。关键是理解每个指标的含义，而不是追求全绿。

---

## 2.4 fastp 清洗：默认参数优先，避免过度裁剪

fastp 是一个高效的 FASTQ 预处理工具，集成了接头检测、质量过滤和长度过滤功能 [^chen2018]。

### 单样本示例

```bash
fastp \
  -i data/raw_data/DRR302386_1.fastq.gz \
  -I data/raw_data/DRR302386_2.fastq.gz \
  -o data/clean_data/DRR302386_1.clean.fastq.gz \
  -O data/clean_data/DRR302386_2.clean.fastq.gz \
  --thread 8 \
  --detect_adapter_for_pe \    # 自动检测 PE 接头
  -h data/fastp_reports/DRR302386.fastp.html \
  -j data/fastp_reports/DRR302386.fastp.json
```

### 多样本批处理（samplesheet 驱动）

```bash
#!/usr/bin/env bash
set -euo pipefail

tail -n +2 metadata/samplesheet.csv | while IFS=',' read -r sample_id run geno cond time rep layout fq1 fq2 batch; do
  f1="data/raw_data/${run}_1.fastq.gz"
  f2="data/raw_data/${run}_2.fastq.gz"

  fastp \
    -i "$f1" -I "$f2" \
    -o "data/clean_data/${run}_1.clean.fastq.gz" \
    -O "data/clean_data/${run}_2.clean.fastq.gz" \
    --thread 8 \
    --detect_adapter_for_pe \
    -h "data/fastp_reports/${run}.html" \
    -j "data/fastp_reports/${run}.json"
done
```

!!! tip "Recommendation"
    先用 fastp 默认参数跑一轮，检查报告后再决定是否需要调整。默认参数已经包含了合理的质量过滤（Q15）和长度过滤（15bp），对大多数数据集足够。不要一上来就设置激进的过滤阈值——过度裁剪会丢失低表达基因的 reads，降低下游检测灵敏度。

---

## 2.5 清洗后复检 + MultiQC 汇总

对清洗后的数据再次运行 FastQC，然后用 MultiQC 将所有报告汇总为一份交互式 HTML：

```bash
# 清洗后 FastQC
fastqc -t 8 -o data/qc/fastqc_clean data/clean_data/*.fastq.gz

# MultiQC 汇总（自动发现 FastQC 和 fastp 报告）
multiqc data \
  --outdir data/multiqc \
  --filename multiqc_report.html
```

MultiQC 报告是你向 PI 和合作者汇报数据质量的主入口。它将所有样本的关键指标整合在同一张图上，让你一眼看出哪些样本偏离了整体分布。

---

## 2.6 如何判定"通过 QC"

没有绝对的阈值标准，但以下经验值可以作为初筛参考：

| 指标 | 建议范围 | 说明 |
| :--- | :--- | :--- |
| Q30 比例 | > 80% | 具体与测序平台相关 |
| Adapter 含量 | 清洗后显著下降 | 不追求绝对为零 |
| 序列长度分布 | 与建库预期一致 | 异常截短需排查 |
| 样本间总 reads | 量级接近 | 极端低深度样本单独评估 |

!!! note
    QC 的目标是"保证数据可分析性"，不是把所有指标修到完美。一个 Q30 = 78% 的样本如果其他指标正常，通常仍然可用；但一个 reads 数只有其他样本 1/10 的样本，即使质量指标全绿，也需要认真评估是否保留。

---

## 2.7 常见翻车点

!!! warning "常见翻车点"
    1. **某样本 reads 数极低**：优先排查下载完整性（MD5 校验）和建库失败，不要盲目继续。极低深度样本会在下游分析中产生不稳定的估计。
    2. **Adapter 清不干净**：fastp 的自动检测偶尔会遗漏非标准接头。如果清洗后 adapter content 仍然偏高，手动指定 `--adapter_sequence` 参数。
    3. **GC 曲线异常**：双峰分布通常提示污染。检查是否混入了其他物种的 DNA，或者文库制备过程中是否有交叉污染。
    4. **重复率偏高**：在高表达组织（如肝脏）或低起始量文库中，较高的重复率可能是正常的。需要结合文库复杂度（如 Picard 的 EstimateLibraryComplexity）综合判断。

---

## 2.8 Key Takeaways

!!! abstract "本章要点"
    1. **QC 是判断步骤，不是清洗步骤**——它的核心目标是评估数据是否可用于统计推断，而不是把所有指标修到全绿。
    2. **先用默认参数，再按需调整**——fastp 的默认参数对大多数数据集已经足够，过度裁剪反而会丢失信号。
    3. **MultiQC 是质量沟通的主入口**——它将所有样本的关键指标整合在一起，让异常样本一目了然。
    4. **GC 异常和极低深度是最需要警惕的信号**——前者可能提示污染，后者会导致下游估计不稳定。
    5. **FastQC 的 PASS/FAIL 不是金标准**——RNA-seq 数据在某些指标上天然会显示 WARN，关键是理解指标含义而非追求全绿。
    6. **所有 QC 决策都要有记录**——哪些样本被剔除、为什么剔除、用了什么过滤参数，这些信息是可复现性的一部分。

---

## 2.9 References

[^chen2018]: Chen S, Zhou Y, Chen Y, Gu J. (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884-i890. DOI: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

[^ewels2016]: Ewels P, Magnusson M, Lundin S, Käller M. (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* 32(19):3047-3048. DOI: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)

