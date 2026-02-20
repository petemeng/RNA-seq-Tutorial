# 📖 Bulk RNA-seq 数据分析最佳实践 (Part 1)

## 引言 (Introduction)

虽然单细胞 RNA 测序 (scRNA-seq) 近年来大放异彩，但传统的 **Bulk RNA-seq (群体细胞转录组测序)** 依然是生物学和临床研究中的中流砥柱。由于其极高的测序深度和较低的噪声，Bulk RNA-seq 依然是发现差异表达基因 (DEG)、可变剪切 (Alternative Splicing) 以及构建基因共表达网络的最可靠手段。

本教程旨在提供一个**标准化、可重复、符合当前业界最佳实践**的 Bulk RNA-seq 分析全流程。

### 🗺️ 标准分析工作流 (The Standard Workflow)

一个标准的 Bulk RNA-seq 流程通常分为两个阶段：**上游分析 (Linux/Bash)** 和 **下游分析 (R/Bioconductor)**。

1. **实验设计与数据获取 (Experimental Design & Data Acquisition)**
2. **质量控制与数据清洗 (Quality Control & Trimming)** -> *FastQC, MultiQC, fastp*
3. **序列比对与定量 (Alignment & Quantification)** -> *STAR (比对) / Salmon (免比对定量)*
4. **探索性数据分析 (Exploratory Data Analysis, EDA)** -> *PCA, 层次聚类 (R)*
5. **差异表达分析 (Differential Expression Analysis)** -> *DESeq2 / edgeR (R)*
6. **功能富集分析 (Functional Enrichment Analysis)** -> *GO, KEGG, GSEA (R)*

---

## 第一章：启航 —— 环境与数据准备 (Environment & References)

在直接下载 Fastq 数据之前，我们必须先解决生物信息学中最令人头疼的问题：**环境配置与版本控制**。

> 📘 **核心概念：可重复性 (Reproducibility)**
> 为了确保你的分析在三年后依然能得出相同的结果，**永远不要在系统的全局环境中直接安装软件**。请使用包管理器 (如 Conda) 构建独立的虚拟环境。

### 1.1 搭建基础计算环境 (Setting up the Environment)

我们将使用 `mamba` 来管理上游分析的工具。

```bash
# 创建环境并安装基础工具
mamba create -n rnaseq_env -c bioconda -c conda-forge \
    fastqc=0.12.1 multiqc=1.19 fastp=0.23.4 \
    star=2.7.11a salmon=1.10.2 samtools=1.19.2

# 激活环境
conda activate rnaseq_env
```

### 1.2 准备参考基因组与注释文件 (Preparing References)

> ⚠️ **避坑指南 (Pitfall)**
> 不要混用不同数据库（如 NCBI, UCSC, Ensembl）的参考文件！**业界最佳实践是：统一从 Ensembl 或 GENCODE 下载配套的 FASTA 和 GTF。**

这里以 **Human (GENCODE 数据库, GRCh38 版本)** 为例：

```bash
mkdir -p ~/RNAseq_project/reference && cd ~/RNAseq_project/reference
# 下载参考基因组和注释文件
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz
```
