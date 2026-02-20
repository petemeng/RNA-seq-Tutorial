# Bulk RNA-seq 数据分析最佳实践 (Part 1)

## 第一章：环境与参考基因组准备（拟南芥主线）

> 核心原则：可重复性来自“固定版本 + 固定参考 + 固定流程”。

---

## 1.1 为什么环境管理必须放在第一步

RNA-seq 结果对软件版本高度敏感：

- 不同版本 STAR/Salmon 可能影响比对率与定量结果。
- 不同 R/Bioconductor 版本可能改变默认行为或统计输出。

因此教程级分析必须做到：

1. 用 `mamba/conda` 建立独立环境。
2. 明确记录工具版本。
3. 导出环境文件，保证异地复现。

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

导出环境：

```bash
mkdir -p envs
conda env export -n rnaseq_env > envs/rnaseq_env.yml
```

> Best Practice：任何图表和结论都应可通过 `rnaseq_env.yml` 在新机器复现。

---

## 1.3 下游环境：R/Bioconductor 包（拟南芥）

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",
  "tximport",
  "apeglm",
  "IHW",
  "clusterProfiler",
  "AnnotationDbi",
  "org.At.tair.db",
  "EnhancedVolcano",
  "pheatmap"
))
```

> 避坑指南：同一项目不要混用多个 R library 路径，避免“本机能跑、服务器报错”。

---

## 1.4 下载拟南芥参考文件（TAIR10 / Ensembl Plants）

```bash
mkdir -p data/reference/arabidopsis_ens58
cd data/reference/arabidopsis_ens58

wget ftp://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-58/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.58.gtf.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-58/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
```

> 只能在同一来源、同一 release 内配对 FASTA + GTF + 转录本序列，避免混源。

---

## 1.5 参考文件一致性检查

```bash
# 解压
pigz -d *.gz

# FASTA 建索引
samtools faidx Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

# GTF 快速检查
grep -v '^#' Arabidopsis_thaliana.TAIR10.58.gtf | head -n 3
grep -v '^#' Arabidopsis_thaliana.TAIR10.58.gtf | awk '$3=="gene"{n++} END{print "gene lines:", n}'
```

### 常见兼容性问题

1. 染色体命名不一致（如 `Chr1` vs `1`）。
2. 基因组版本与注释版本不一致。
3. 转录本 FASTA 和 GTF 的 transcript ID 对不上。

---

## 1.6 构建 STAR 与 Salmon 索引（拟南芥）

### STAR 基因组索引

```bash
mkdir -p data/reference/star_index
STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir data/reference/star_index \
  --genomeFastaFiles data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  --sjdbGTFfile data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.58.gtf \
  --sjdbOverhang 149
```

### Salmon 转录本索引

```bash
mkdir -p data/reference/salmon_index
salmon index \
  -t data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.cdna.all.fa \
  -i data/reference/salmon_index \
  -k 31
```

`sjdbOverhang` 建议取 `read_length - 1`（例如 PE150 用 149）。

---

## 1.7 建议记录参考版本清单

在 `results/reference_manifest.txt` 写入：

```text
Species: Arabidopsis thaliana
Genome: Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
Annotation: Arabidopsis_thaliana.TAIR10.58.gtf
Transcriptome: Arabidopsis_thaliana.TAIR10.cdna.all.fa
STAR: 2.7.11b
Salmon: 1.10.2
```

这份清单是复现、审稿回复和团队协作的关键证据。

---

## 1.8 常见翻车点

1. 混用 TAIR 与其他数据库注释，导致 ID 映射混乱。
2. 复用旧索引到新项目，造成隐性错误。
3. 未记录版本，后续无法解释结果差异。
4. 参考文件未校验就开始比对，返工成本极高。

---

## 1.9 本章检查清单

- `rnaseq_env.yml` 已保存。
- TAIR10 FASTA/GTF/cDNA 文件来源一致。
- STAR 与 Salmon 索引已成功构建。
- `reference_manifest.txt` 已写入结果目录。

完成后进入 QC 与清洗章节。
