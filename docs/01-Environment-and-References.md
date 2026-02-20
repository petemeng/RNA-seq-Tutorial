# Bulk RNA-seq 数据分析最佳实践 (Part 1)

## 第一章：环境与参考基因组准备 (Environment & References)

> 核心原则：可重复性来自“固定版本 + 固定参考 + 固定流程”，而不是“这次能跑通”。

---

## 1.1 为什么环境管理是第一步

RNA-seq 的结果对软件版本很敏感：

- 不同版本 STAR/Salmon 可能带来不同的比对率或定量偏差。
- 不同 R/Bioconductor 版本可能改变统计结果或默认参数。

因此，任何教程级分析都应该：

1. 用 `mamba/conda` 建立独立环境。
2. 明确写出工具版本。
3. 保留环境导出文件。

---

## 1.2 上游环境：Mamba 一次安装

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

导出环境（必须留档）：

```bash
mkdir -p envs
conda env export -n rnaseq_env > envs/rnaseq_env.yml
```

> Best Practice：任何图表和结论都应该可以在新机器上仅靠 `rnaseq_env.yml` 复现。

---

## 1.3 下游环境：R/Bioconductor 推荐包

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "tximport",
  "apeglm",
  "IHW",
  "clusterProfiler",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "org.At.tair.db",
  "EnhancedVolcano",
  "pheatmap"
))
```

> 避坑指南：同一个项目内不要混用多个 R library 路径，避免“本机能跑，服务器报错”。

---

## 1.4 参考基因组下载：只用配套文件

### 人类示例（GENCODE GRCh38）

```bash
mkdir -p data/reference/human_gencode_v45
cd data/reference/human_gencode_v45

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.transcripts.fa.gz
```

### 拟南芥示例（Ensembl Plants）

```bash
mkdir -p data/reference/arabidopsis_ens58
cd data/reference/arabidopsis_ens58

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-58/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gtf.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
```

> 只能在同一来源、同一 release 内配对 FASTA + GTF + 转录本序列。不要混源。

---

## 1.5 参考文件一致性检查

下载完成后，先做完整性和格式检查。

```bash
# 解压
pigz -d *.gz

# FASTA 建索引
samtools faidx GRCh38.primary_assembly.genome.fa

# GTF 行数和注释类型快速检查
grep -v '^#' gencode.v45.annotation.gtf | head -n 3
grep -v '^#' gencode.v45.annotation.gtf | awk '$3=="gene"{n++} END{print "gene lines:", n}'
```

### 常见兼容性问题

1. FASTA 染色体名带 `chr`，GTF 不带（或反过来）。
2. 基因组版本与注释版本不一致（如 GRCh37 + gencode v45）。
3. 转录本 FASTA 和 GTF 的 transcript ID 不一致。

---

## 1.6 构建 STAR 与 Salmon 索引

### STAR 基因组索引

```bash
mkdir -p data/reference/star_index
STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir data/reference/star_index \
  --genomeFastaFiles data/reference/human_gencode_v45/GRCh38.primary_assembly.genome.fa \
  --sjdbGTFfile data/reference/human_gencode_v45/gencode.v45.annotation.gtf \
  --sjdbOverhang 149
```

`sjdbOverhang` 建议取 `read_length - 1`（例如 PE150 用 149）。

### Salmon 转录本索引

```bash
mkdir -p data/reference/salmon_index
salmon index \
  -t data/reference/human_gencode_v45/gencode.v45.transcripts.fa \
  -i data/reference/salmon_index \
  -k 31
```

---

## 1.7 为什么“参考版本号”必须写进结果

建议在结果目录保留 `results/reference_manifest.txt`：

```text
Species: Homo sapiens
Genome: GRCh38.primary_assembly.genome.fa
Annotation: gencode.v45.annotation.gtf
Transcriptome: gencode.v45.transcripts.fa
STAR: 2.7.11b
Salmon: 1.10.2
```

这份文件是后续复现、审稿回复和团队协作的关键证据。

---

## 1.8 常见翻车点

1. 混用 UCSC/Ensembl 命名体系，导致大量 read 无法注释到基因。
2. 把旧项目残留索引直接复用到新数据。
3. 没记录工具版本，后续无法定位结果差异来源。
4. 未校验参考文件就开始比对，后期代价极高。

---

## 1.9 本章检查清单

- `rnaseq_env.yml` 已保存。
- FASTA/GTF/转录本文件来源一致、版本一致。
- STAR 与 Salmon 索引构建成功。
- `reference_manifest.txt` 已写入项目。

完成后再进入 QC 与清洗章节，后续流程会稳定很多。
