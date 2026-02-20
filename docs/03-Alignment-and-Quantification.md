# Bulk RNA-seq 数据分析最佳实践 (Part 3)

## 第三章：序列比对与定量 (Alignment & Quantification)

> 核心原则：工具选择要由“生物学问题”驱动，而不是由“谁更快”驱动。

---

## 3.1 STAR 与 Salmon：怎么选

### 用 STAR 的典型场景

- 你要分析可变剪切、融合基因、新转录本。
- 你需要 BAM 做可视化或下游结构分析。

### 用 Salmon 的典型场景

- 你主要做差异表达与富集分析。
- 你追求速度、资源效率和偏差校正。

> Best Practice：常见方案是“Salmon 做表达分析 + STAR 保留关键样本做结构验证”。

---

## 3.2 STAR 比对标准流程（推荐 2-pass）

### 1) 构建索引

```bash
STAR --runThreadN 8 --runMode genomeGenerate \
  --genomeDir data/reference/star_index \
  --genomeFastaFiles data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  --sjdbGTFfile data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.58.gtf \
  --sjdbOverhang 149
```

### 2) 比对

```bash
STAR --runThreadN 8 \
  --genomeDir data/reference/star_index \
  --readFilesIn data/clean_data/SRR8694017_1.clean.fastq.gz data/clean_data/SRR8694017_2.clean.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix data/alignments/SRR8694017. \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \
  --quantMode GeneCounts
```

### 3) 基本质控

```bash
samtools flagstat data/alignments/SRR8694017.Aligned.sortedByCoord.out.bam > data/alignments/SRR8694017.flagstat.txt
```

重点看：

- 唯一比对率（过低通常提示参考不匹配、污染或链参数错误）
- 多重比对比例（重复序列/基因家族常见）
- 未比对原因（too many mismatches / too short 等）

---

## 3.3 Salmon 定量标准流程

### 1) 构建索引

```bash
salmon index \
  -t data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.cdna.all.fa \
  -i data/reference/salmon_index \
  -k 31
```

### 2) 定量（PE 示例）

```bash
salmon quant \
  -i data/reference/salmon_index \
  -l A \
  -1 data/clean_data/SRR8694017_1.clean.fastq.gz \
  -2 data/clean_data/SRR8694017_2.clean.fastq.gz \
  --validateMappings \
  --gcBias \
  --seqBias \
  -p 8 \
  -o data/quant/SRR8694017
```

> `-l A` 会自动推断文库链特异性，适合教程和未知文库类型场景。

---

## 3.4 链特异性：最常见的致命坑

链类型填错会让表达量系统性偏差，严重时接近“全基因低表达”。

可用两种方法确认：

1. 询问测序平台（最可靠）。
2. 用工具自动推断（Salmon `-l A` 或 RSeQC）。

常见链类型标记：

- `ISR`：常见于 dUTP 类反向链文库
- `ISF`：正向链文库
- `IU`：非链特异

---

## 3.5 计数矩阵准备：从转录本到基因层

如果你用 Salmon，建议在 R 里用 `tximport` 汇总到基因层再进 DESeq2。

```r
library(tximport)
library(readr)

samples <- read_csv("metadata/samplesheet.csv")
files <- file.path("data/quant", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id

# tx2gene 需提前准备：两列 transcript_id, gene_id
tx2gene <- read_tsv("metadata/tx2gene.tsv", col_names = c("TXNAME", "GENEID"))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
```

> Best Practice：用 `tximport` 处理长度偏差，比直接读 TPM 稳定。

---

## 3.6 如何判断比对/定量是否合格

经验上可先做以下检查：

1. 样本间映射率是否一致（异常样本需单独排查）。
2. 多数样本是否达到预期测序深度。
3. 是否存在极端低复杂度样本（大量重复/少量基因占据 reads）。
4. 计数矩阵中是否有大比例全零样本或全零基因。

---

## 3.7 常见翻车点

1. FASTA 和 GTF 版本不配，导致注释失败或计数偏差。
2. 链特异性参数错误，造成全局表达偏移。
3. 混用 STAR geneCounts 与 Salmon quant 结果做同一统计模型。
4. 忽略比对日志，直接进入差异分析。

---

## 3.8 本章检查清单

- 每个样本都有可追溯的比对/定量日志。
- 已确认文库链类型或使用自动推断并记录。
- 计数矩阵准备方式明确（STAR counts 或 Salmon+tximport）。
- 异常样本有处理结论（剔除、保留、重测建议）。

完成后进入标准化与 EDA 章节。

---

## 3.9 本教程实跑代码与结果（PRJDB11848）

以下命令已在教程环境实跑通过：

```bash
# 生成样本表
bash scripts/01_prepare_prjdb11848_samplesheet.sh

# 下载 FASTQ
bash scripts/02_download_fastq.sh

# Salmon 定量
bash scripts/03_quantify_salmon.sh
```

对应产物：

- `validation_run_downstream/metadata/samplesheet.csv`（36 样本，含表头共 37 行）
- `validation_run_downstream/data/raw_data/*.fastq.gz`（72 个文件）
- `validation_run_downstream/data/quant/<sample_id>/quant.sf`（36 个样本目录）
- `validation_run_downstream/logs/salmon_quant_<sample_id>.log`

本次实跑 mapping rate（36 样本）：

- `min = 63.953%`
- `mean = 95.231%`
- `max = 97.471%`
