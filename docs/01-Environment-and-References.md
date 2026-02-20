# 1. 环境与参考基因组准备 (Environment and References)

## 1.1 Motivation

RNA-seq 分析结果对软件版本和参考文件高度敏感。同一份数据用不同版本的 STAR 比对，唯一比对率可能相差数个百分点；不同版本的 GTF 注释文件包含的基因数量和转录本结构都不同；甚至 R/Bioconductor 的小版本更新都可能改变 DESeq2 的默认行为。如果不在分析开始前锁定环境，你的结果就无法被他人（包括三个月后的你自己）复现。

更隐蔽的问题是参考文件不一致。基因组 FASTA 来自 Ensembl release 58，GTF 来自 release 55，转录本序列来自 TAIR——这种混源会导致染色体命名不匹配、基因 ID 映射失败、计数矩阵出现大量零值，而且这些错误往往不会产生明确的报错信息，只会在下游分析中表现为"结果看起来不太对"。

本章建立两个基础：第一，用 Conda/Mamba 创建隔离的分析环境并导出环境文件；第二，从同一数据源的同一 release 下载配套的参考基因组、注释文件和转录本序列，并构建比对和定量所需的索引。

---

## 1.2 环境与数据 (Environment and Data)

本章需要以下工具：

- Conda 或 Mamba（环境管理）
- R + BiocManager（下游分析包安装）
- wget（下载参考文件）
- samtools（FASTA 索引）
- STAR 和 Salmon（构建比对/定量索引）

---

## 1.3 上游环境：Mamba 一次安装

我们推荐用 Mamba（Conda 的高性能替代）创建独立环境，将所有上游工具锁定在同一环境中：

```bash
mamba create -n rnaseq_env -c conda-forge -c bioconda \
  fastqc=0.12.1 \
  multiqc=1.19 \
  fastp=0.23.4 \
  star=2.7.11b \
  salmon=1.10.2 \
  samtools=1.19.2 \
  subread=2.0.6 \
  sra-tools=3.1.1

conda activate rnaseq_env
```

安装完成后立即导出环境文件：

```bash
mkdir -p envs
conda env export -n rnaseq_env > envs/rnaseq_env.yml
```

!!! tip "Recommendation"
    任何图表和结论都应可通过 `rnaseq_env.yml` 在新机器上复现。将这个文件纳入版本控制，与分析代码一起提交。

---

## 1.4 下游环境：R/Bioconductor 包

下游统计分析在 R 中完成。以拟南芥为例，需要安装以下核心包：

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "tximport",
  "apeglm",
  "IHW",
  "clusterProfiler",
  "AnnotationDbi",
  "org.At.tair.db",      # 拟南芥注释数据库
  "EnhancedVolcano",
  "pheatmap"
))
```

!!! warning "避坑指南"
    同一项目不要混用多个 R library 路径。常见场景是本地用个人 library 安装了新版包，但服务器上用的是系统 library 的旧版——这会导致"本机能跑、服务器报错"。建议用 `renv` 或在脚本开头记录 `sessionInfo()`。

---

## 1.5 下载参考文件：同源同版本

参考文件必须从同一数据源的同一 release 下载，确保基因组、注释和转录本序列完全配套。本教程使用 Ensembl Plants release 58 的拟南芥 TAIR10 参考：

```bash
mkdir -p data/reference/arabidopsis_ens58
cd data/reference/arabidopsis_ens58

# 基因组 FASTA
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# 基因注释 GTF
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-58/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gtf.gz

# 转录本序列（Salmon 定量用）
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
```

!!! warning
    只能在同一来源、同一 release 内配对 FASTA + GTF + 转录本序列。混用不同来源（如 TAIR 官网的基因组 + Ensembl 的 GTF）会导致染色体命名不一致（`Chr1` vs `1`）、基因 ID 映射失败等问题。

---

## 1.6 参考文件一致性检查

下载后先做基本验证，确认文件完整且相互兼容：

```bash
# 解压
pigz -d *.gz

# 构建 FASTA 索引
samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

# 检查 GTF 格式和基因数量
grep -v '^#' Arabidopsis_thaliana.TAIR10.58.gtf | head -n 3
grep -v '^#' Arabidopsis_thaliana.TAIR10.58.gtf | awk '$3=="gene"{n++} END{print "gene lines:", n}'
```

### 常见兼容性问题

!!! warning "常见翻车点"
    1. **染色体命名不一致**：FASTA 用 `1, 2, 3...`，GTF 用 `Chr1, Chr2, Chr3...`——STAR 比对时会找不到参考序列。
    2. **基因组版本与注释版本不匹配**：基因组是 TAIR10 但 GTF 是 Araport11，坐标体系不同会导致大量 reads 无法定位到基因。
    3. **转录本 FASTA 与 GTF 的 transcript ID 对不上**：Salmon 定量后用 tximport 汇总时会丢失大量转录本。

---

## 1.7 构建 STAR 与 Salmon 索引

### STAR 基因组索引

```bash
mkdir -p data/reference/star_index

STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir data/reference/star_index \
  --genomeFastaFiles data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  --sjdbGTFfile data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.58.gtf \
  --sjdbOverhang 149   # read_length - 1，PE150 用 149
```

### Salmon 转录本索引

```bash
mkdir -p data/reference/salmon_index

salmon index \
  -t data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.cdna.all.fa \
  -i data/reference/salmon_index \
  -k 31
```

!!! note
    `sjdbOverhang` 建议取 `read_length - 1`。对于 PE150 数据设为 149；如果你的数据是 PE100，应设为 99。这个参数影响剪接位点检测的灵敏度。

---

## 1.8 记录参考版本清单

在项目目录中创建一份版本清单，记录所有参考文件和工具版本：

```text
Species: Arabidopsis thaliana
Genome: Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
Annotation: Arabidopsis_thaliana.TAIR10.58.gtf
Transcriptome: Arabidopsis_thaliana.TAIR10.cdna.all.fa
Source: Ensembl Plants release 58
STAR: 2.7.11b
Salmon: 1.10.2
```

!!! tip "Recommendation"
    将这份清单保存为 `results/reference_manifest.txt` 并纳入版本控制。它是复现分析、回复审稿意见和团队协作的关键证据。

---

## 1.9 常见翻车点

!!! warning "常见翻车点"
    1. **混用 TAIR 与其他数据库注释**：TAIR 官网的基因 ID 格式、Ensembl 的基因 ID 格式和 NCBI 的基因 ID 格式各不相同，混用会导致下游 ID 映射混乱。
    2. **复用旧索引到新项目**：上一个项目用 Ensembl release 55 建的索引，新项目的 GTF 是 release 58——索引和注释不匹配会产生隐性错误。每个项目应重新构建索引。
    3. **未记录版本**：三个月后审稿人问"你用的什么版本的注释"，你答不上来。
    4. **参考文件未校验就开始比对**：FASTA 下载不完整或 GTF 格式损坏，比对阶段才报错，返工成本极高。

---

## 1.10 Key Takeaways

!!! abstract "本章要点"
    1. **环境隔离是可复现性的基础**——用 Conda/Mamba 创建独立环境，导出 `rnaseq_env.yml`，确保任何人都能重建完全相同的分析环境。
    2. **参考文件必须同源同版本**——基因组 FASTA、GTF 注释和转录本序列必须来自同一数据源的同一 release，混源是隐性错误的主要来源。
    3. **下载后立即验证**——检查染色体命名一致性、基因数量合理性、转录本 ID 匹配性，不要等到比对阶段才发现问题。
    4. **索引按项目构建**——不要跨项目复用旧索引，版本不匹配的风险远大于重建索引的时间成本。
    5. **版本清单是必备产物**——`reference_manifest.txt` 和 `sessionInfo()` 是复现性三件套的核心组成部分。

---

## 1.11 References

[^patro2017]: Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods* 14:417-419. DOI: [10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)

[^dobin2013]: Dobin A, Davis CA, Schlesinger F, et al. (2013) STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29(1):15-21. DOI: [10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

