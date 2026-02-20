# 📖 Bulk RNA-seq 数据分析最佳实践 (Part 2)

## 第二章：数据获取与质量控制 (Data Acquisition & QC)

### 2.1 优雅地获取原始数据 (Downloading Raw Data)

我们在论文中常看到数据存放在 NCBI 的 **SRA (Sequence Read Archive)** 数据库中。许多老旧教程会教你使用 `sratoolkit` 中的 `fastq-dump` 来下载 `.sra` 文件并转化为 `.fastq`。

> 🛑 **避坑指南 (Pitfall)：放弃 fastq-dump！**
> `fastq-dump` 速度奇慢，且经常断点失败。
> **业界最佳实践：使用 ENA (European Nucleotide Archive) 数据库。** ENA 与 NCBI SRA 数据是实时同步的，但 ENA 直接提供打包好的 `.fastq.gz` 文件下载链接！

为了演示，我们从 Airway 数据集中选取 2 个样本（双端测序 Paired-end）：
* `SRR1039508` (Control 组)
* `SRR1039509` (Treated 组)

```bash
# 构建项目目录结构
mkdir -p ~/RNAseq_project/{raw_data,qc,cleandata}
cd ~/RNAseq_project/raw_data

# 使用 wget 下载 ENA 上的 FASTQ 文件
# Control Sample
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/008/SRR1039508/SRR1039508_2.fastq.gz

# Treated Sample
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039509/SRR1039509_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/009/SRR1039509/SRR1039509_2.fastq.gz
```

### 2.2 质量控制 (Quality Control)

拿到数据后，第一步是使用 `FastQC` 查看数据质量。

> 📘 **核心概念：Phred 质量得分 (Phred Quality Score, Q-score)**
> 测序仪每测出一个碱基 (A/T/C/G)，都会给它一个“自信度”打分。
> 公式：$Q = -10 \times \log_{10}(P)$，其中 P 是测错的概率。
> *   **Q20** = 错误率 1% (准确率 99%)
> *   **Q30** = 错误率 0.1% (准确率 99.9%)
> **金标准：** 优秀的测序数据，绝大多数碱基的质量得分应在 **Q30 以上**。

```bash
# 运行 FastQC
# -t 4 表示使用 4 个 CPU 线程
fastqc -t 4 -o ../qc *.fastq.gz
```

运行后，FastQC 会生成 HTML 报告。你需要重点关注以下三个指标：
1. **Per base sequence quality**: 纵坐标就是 Q-score。箱线图如果掉入红色区域 (Q < 20)，说明测序质量极差。
2. **Adapter Content**: 检查是否残留了测序接头 (Adapter)。
3. **Sequence Duplication Levels**: Bulk RNA-seq 中一定程度的重复是正常的。

### 2.3 数据清洗与过滤 (Data Trimming)

传统的做法是使用 `Trimmomatic` 去接头。但今天，我们推荐使用**当前最先进的工具：`fastp`**。
它由 C++ 编写，速度极快，一步到位完成过滤、去接头、甚至生成清洗前后的双份网页报告。

```bash
cd ~/RNAseq_project/
mkdir -p cleandata fastp_reports

# 对 Control 样本进行清洗
fastp \
  -i raw_data/SRR1039508_1.fastq.gz -I raw_data/SRR1039508_2.fastq.gz \
  -o cleandata/SRR1039508_1_clean.fq.gz -O cleandata/SRR1039508_2_clean.fq.gz \
  --thread 4 \
  --detect_adapter_for_pe \
  -h fastp_reports/SRR1039508_fastp.html \
  -j fastp_reports/SRR1039508_fastp.json
```

> 💡 **Best Practice：过度裁剪 (Over-trimming) 是有害的！**
> 新手常犯的错误是为了让 FastQC 报告好看，设置极其严格的过滤条件。
> **现代比对软件 (如 STAR) 具有“软裁切 (Soft-clipping)”功能**，它们足够聪明，能在比对时自动忽略少量低质量碱基。因此，`fastp` 使用默认参数去除接头和极低质量碱基即可，**保留更多数据意味着保留更多生物学真实信号**。

### 2.4 化繁为简：聚合报告 (Aggregating Reports with MultiQC)

想象一下，如果你有 100 个样本，你要点开 100 个 FastQC 的 HTML 吗？
**`MultiQC`** 解决了这个痛点。它会扫描目录下所有支持的日志文件，并将它们合并成一个精美、可交互的单一网页报告。

```bash
# 回到项目根目录并运行
multiqc . -o multiqc_report/
```

执行完毕后，打开 `multiqc_report/multiqc_report.html`。你可以在一张图里看到所有样本的测序深度、Q30 比例等。
