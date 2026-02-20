# 3. 序列比对与定量 (Alignment and Quantification)

## 3.1 Motivation

比对与定量是将原始测序 reads 转化为基因表达计数矩阵的关键步骤——这个矩阵是后续所有统计分析的输入。这一步的选择直接影响你能检测到多少基因、表达量估计的准确性以及下游差异分析的统计功效。

这里有一个根本性的方法选择：传统的"比对到基因组 → 计数"流程（STAR + featureCounts）和基于准比对的轻量级定量流程（Salmon）。两者各有适用场景，选择应由你的生物学问题驱动，而不是由"谁更快"驱动。如果你需要 BAM 文件做可变剪切分析或 IGV 可视化，STAR 是必须的；如果你主要做差异表达和富集分析，Salmon 在速度、资源消耗和偏差校正方面都有优势 [^patro2017]。

本章将分别介绍两种流程的标准操作，并在 PRJDB11848 数据集上使用 Salmon 完成定量。36 个样本的平均 mapping rate 为 95.2%，确认数据质量良好。

---

## 3.2 环境与数据 (Environment and Data)

本章需要以下工具（均已包含在 [第 1 章](01-Environment-and-References.md) 的 `rnaseq_env` 环境中）：

- STAR（基因组比对）
- Salmon（转录本准比对定量）
- samtools（BAM 处理）
- featureCounts / Subread（基因计数）
- R + tximport（Salmon 结果汇总）

输入文件：`data/clean_data/*.fastq.gz`（清洗后的 FASTQ）
输出文件：BAM 文件（STAR 流程）或 `quant.sf`（Salmon 流程）

---

## 3.3 工具选择：STAR vs Salmon

| 维度 | STAR + featureCounts | Salmon |
| :--- | :--- | :--- |
| 原理 | 比对到基因组 → 按注释计数 | 准比对到转录本 → 直接定量 |
| 速度 | 较慢（需写 BAM） | 快 5-10 倍 |
| 内存 | 高（人类 ~30GB） | 低（~8GB） |
| 偏差校正 | 无内置 | GC bias、序列偏差、长度偏差 |
| 输出 | BAM + 计数矩阵 | 转录本 TPM + count |
| 适用场景 | 需要 BAM 的分析（剪切、可视化） | 差异表达、富集分析 |

!!! tip "Recommendation"
    对于以差异表达为主要目标的分析，我们推荐使用 **Salmon** [^patro2017]。它不仅速度快，还内置了 GC bias 和序列偏差校正，配合 tximport 可以直接输入 DESeq2。常见的折中方案是"Salmon 做表达定量 + STAR 保留关键样本做结构验证"。

---

## 3.4 STAR 比对标准流程

### 构建索引

如果你在 [第 1 章](01-Environment-and-References.md) 已经构建了 STAR 索引，可以跳过此步。

```bash
STAR --runThreadN 8 --runMode genomeGenerate \
  --genomeDir data/reference/star_index \
  --genomeFastaFiles data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  --sjdbGTFfile data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.58.gtf \
  --sjdbOverhang 149
```

### 比对（2-pass 模式）

```bash
STAR --runThreadN 8 \
  --genomeDir data/reference/star_index \
  --readFilesIn data/clean_data/DRR302386_1.clean.fastq.gz data/clean_data/DRR302386_2.clean.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix data/alignments/DRR302386. \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \       # 2-pass 模式提高新剪接位点检测灵敏度
  --quantMode GeneCounts       # 同时输出基因计数
```

### 比对质量检查

```bash
samtools flagstat data/alignments/DRR302386.Aligned.sortedByCoord.out.bam \
  > data/alignments/DRR302386.flagstat.txt
```

重点关注：

- **唯一比对率**：过低（< 70%）通常提示参考基因组不匹配、外源污染或链参数错误。
- **多重比对比例**：基因家族丰富的物种中较高是正常的。
- **未比对 reads**：检查 STAR 日志中的未比对原因（too many mismatches / too short）。

---

## 3.5 Salmon 定量标准流程

### 构建索引

```bash
salmon index \
  -t data/reference/arabidopsis_ens58/Arabidopsis_thaliana.TAIR10.cdna.all.fa \
  -i data/reference/salmon_index \
  -k 31
```

### 定量（PE 示例）

```bash
salmon quant \
  -i data/reference/salmon_index \
  -l A \                        # 自动推断文库链特异性
  -1 data/clean_data/DRR302386_1.clean.fastq.gz \
  -2 data/clean_data/DRR302386_2.clean.fastq.gz \
  --validateMappings \          # 启用选择性比对验证
  --gcBias \                    # GC 偏差校正
  --seqBias \                   # 序列偏差校正
  -p 8 \
  -o data/quant/DRR302386
```

!!! note
    `-l A` 让 Salmon 自动推断文库的链特异性类型，适合教程场景和未知文库类型。推断结果会记录在 `aux_info/lib_format_counts.json` 中，建议检查确认。

---

## 3.6 链特异性：最常见的致命坑

链特异性参数填错会让表达量产生系统性偏差，严重时接近"全基因低表达"。

确认链特异性的两种方法：

1. **询问测序平台**（最可靠）：建库方法决定了链特异性类型。
2. **工具自动推断**：Salmon `-l A` 或 RSeQC 的 `infer_experiment.py`。

常见链类型标记：

| Salmon 标记 | 含义 | 常见建库方法 |
| :--- | :--- | :--- |
| `ISR` | 反向链文库 | dUTP、Illumina TruSeq Stranded |
| `ISF` | 正向链文库 | 较少见 |
| `IU` | 非链特异 | 早期 RNA-seq |

!!! warning
    链类型填错不会让 pipeline 报错——它只会让你的计数矩阵系统性偏低。如果你发现所有样本的 mapping rate 都异常低，首先检查链特异性参数。

---

## 3.7 计数矩阵准备：从转录本到基因层

Salmon 输出的是转录本级别的定量结果。进入 DESeq2 之前，需要用 `tximport` 将转录本计数汇总到基因层 [^soneson2015]：

```r
library(tximport)
library(readr)

samples <- read_csv("metadata/samplesheet.csv")
files <- file.path("data/quant", samples$run, "quant.sf")
names(files) <- samples$sample_id

# tx2gene：两列映射表（transcript_id → gene_id）
tx2gene <- read_tsv("metadata/tx2gene.tsv", col_names = c("TXNAME", "GENEID"))

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"  # 校正转录本长度偏差
)
```

!!! tip "Recommendation"
    使用 `tximport` 配合 `countsFromAbundance = "lengthScaledTPM"` 处理转录本长度偏差，比直接读取 TPM 或原始 count 更适合 DESeq2 的统计模型 [^soneson2015]。

---

## 3.8 PRJDB11848 实跑结果

本教程案例数据使用 Salmon 完成定量：

```bash
bash scripts/03_quantify_salmon.sh
```

36 个样本的 mapping rate 统计：

| 指标 | 值 |
| :--- | :--- |
| 最小值 | 63.95% |
| 平均值 | 95.23% |
| 最大值 | 97.47% |

!!! note
    一个样本的 mapping rate 为 63.95%，明显低于其他样本。在实际项目中应进一步排查原因（如污染、文库质量）。本教程保留该样本以展示完整流程，但在正式分析中建议评估其对下游结果的影响。

---

## 3.9 常见翻车点

!!! warning "常见翻车点"
    1. **FASTA 和 GTF 版本不配**：基因组是 TAIR10 但注释是 Araport11，坐标体系不同导致大量 reads 无法定位到基因。参见 [第 1 章](01-Environment-and-References.md) 的一致性检查。
    2. **链特异性参数错误**：造成全局表达量系统性偏低，但不会报错。
    3. **混用不同定量工具的结果**：STAR geneCounts 和 Salmon quant 的计数含义不同（前者是原始计数，后者经过偏差校正），不能混入同一个统计模型。
    4. **忽略比对日志**：直接跳到差异分析，错过了 mapping rate 异常低的样本。

---

## 3.10 Key Takeaways

!!! abstract "本章要点"
    1. **工具选择由生物学问题驱动**——需要 BAM 做结构分析用 STAR，主要做差异表达用 Salmon。
    2. **Salmon + tximport 是当前 bulk RNA-seq 差异表达分析的推荐定量方案**——速度快、内置偏差校正、与 DESeq2 无缝衔接 [^patro2017]。
    3. **链特异性是最常见的致命错误**——填错不会报错，但会让计数矩阵系统性偏低。务必通过建库信息或自动推断确认。
    4. **比对/定量日志是质控的重要组成部分**——mapping rate 异常的样本需要在进入差异分析前排查原因。
    5. **tximport 的长度校正优于直接使用 TPM**——它在保留计数数据统计特性的同时校正了转录本长度偏差。

---

## 3.11 References

[^patro2017]: Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. *Nature Methods* 14:417-419. DOI: [10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)

[^dobin2013]: Dobin A, Davis CA, Schlesinger F, et al. (2013) STAR: ultrafast universal RNA-seq aligner. *Bioinformatics* 29(1):15-21. DOI: [10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

[^soneson2015]: Soneson C, Love MI, Robinson MD. (2015) Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research* 4:1521. DOI: [10.12688/f1000research.7563.2](https://doi.org/10.12688/f1000research.7563.2)

---

## 3.12 Contributors

### Authors

- *(待补充)*

### Reviewers

- *(待补充)*
