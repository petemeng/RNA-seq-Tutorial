# 0. 数据获取与元数据管理 (Data Acquisition)

## 0.1 Motivation

一个 RNA-seq 项目的分析上限，往往不是由算法决定的，而是由元数据的质量决定的。样本标签错了，后面所有的统计推断都在回答一个错误的问题；分组信息缺失了，批次效应就无法校正；重复数不够，差异检验的统计功效就不足以发现真实的生物学信号。这些问题一旦进入下游流程，几乎不可能被自动检测到——你会得到一个"看起来正常"但结论完全错误的分析结果。

元数据错误的危害在于它的隐蔽性。一个标签互换的样本不会让 pipeline 报错，它只会让你的 PCA 图上出现一个"离群点"，而你可能会选择删掉它而不是追查原因。一个遗漏的批次变量不会阻止 DESeq2 运行，它只会让你的差异基因列表里混入大量假阳性。这类错误的修复成本极高——往往意味着从头返工。

本章建立两个基础：第一，如何从公共数据库高效获取原始测序数据并验证完整性；第二，如何构建一份结构化的元数据表（samplesheet），使其成为驱动整个分析流程的"单一真实来源"（Single Source of Truth）。后续所有章节的代码都将从这份 samplesheet 读取样本信息，而不是在脚本中硬编码样本名。

---

## 0.2 先定义问题，再下载数据

在打开终端敲下第一条 `wget` 之前，先回答三个问题：

1. **你要比较什么？** 明确主效应（例如 AvrRpm1 vs mock）和可能的交互效应（例如 genotype × condition）。
2. **你的最小生物学重复数是多少？** 每组至少 3 个生物学重复是基本要求；Schurch 等人的 benchmark 表明，要可靠检测中等效应量的差异基因，每组 6 个重复更为理想 [^schurch2016]。
3. **你的主要干扰因素是什么？** 批次、时间点、组织来源、文库制备方法——这些都需要记录在元数据中，以便后续建模时纳入。

!!! tip "Recommendation"
    在下载任何数据之前，先写一份简短的"分析设计说明"（1 页即可），明确实验问题、分组、重复数和已知的混杂因素。这份文档将帮助你在后续分析中保持方向感，也是与合作者沟通的基础。

---

## 0.3 数据来源：SRA、ENA、GEO

公共 RNA-seq 数据通常通过三个入口获取：GEO（Gene Expression Omnibus）、SRA（Sequence Read Archive）和 ENA（European Nucleotide Archive）[^leinonen2011]。它们之间的关系是：

- 论文通常提供 **GEO 编号**（如 `GSE123456`），GEO 页面会链接到对应的 SRA 项目（如 `SRP/PRJNA`）。
- **SRA** 和 **ENA** 的数据是同步的，但下载体验不同。
- 对于日本 DDBJ 提交的数据（如本教程使用的 PRJDB11848），三个数据库均可检索到。

!!! tip "Recommendation"
    我们推荐优先使用 **ENA** 进行批量下载，原因有三：（1）ENA 直接提供 `.fastq.gz` 的 FTP/HTTP 链接，无需像 SRA 那样先下载 `.sra` 再转换；（2）元数据表结构清晰，可直接导出为 TSV 用于构建 samplesheet；（3）对脚本化批量下载更友好。

### 从 ENA 获取下载链接

访问 [ENA Browser](https://www.ebi.ac.uk/ena/browser/)，输入项目编号（如 `PRJDB11848`），在结果页面选择 "TSV" 格式导出，勾选 `run_accession`、`sample_title`、`fastq_ftp` 等字段。导出的文件就是构建 samplesheet 的起点。

---

## 0.4 元数据：单一真实来源

整个分析流程应由一份 `samplesheet.csv` 驱动，而不是在每个脚本里硬编码样本名。这份表就是项目的"单一真实来源"——所有下游脚本从它读取样本信息、分组信息和文件路径。

### 最小字段规范

| 列名 | 是否必填 | 说明 |
| :--- | :--- | :--- |
| `sample_id` | 是 | 样本唯一标识，不可重复 |
| `run` | 是 | SRR/ERR/DRR 编号 |
| `genotype` | 视设计 | 基因型（如 WT / clf） |
| `condition` | 是 | 主要处理因素（如 mock / AvrRpm1） |
| `time` | 视设计 | 时间点（如 0h / 0.5h / 3h） |
| `replicate` | 是 | 生物学重复编号 |
| `layout` | 是 | `PAIRED` 或 `SINGLE` |
| `fastq_1` | 是 | R1 文件路径或 URL |
| `fastq_2` | PE 必填 | R2 文件路径或 URL |
| `batch` | 建议 | 文库批次 / 上机批次 |

### 本教程案例：PRJDB11848

本教程使用的案例数据来自 PRJDB11848，这是一个拟南芥（*Arabidopsis thaliana*）免疫转录组实验，采用 **2 × 2 × 3 因子设计**：

- **基因型**：WT（野生型）、clf（*CURLY LEAF* 突变体，表观遗传调控因子）
- **处理**：mock（对照）、AvrRpm1（细菌效应蛋白，触发 ETI 免疫响应）
- **时间点**：0h、0.5h、3h
- **每组 3 个生物学重复**，共 2 × 2 × 3 × 3 = **36 个样本**

实际 samplesheet 示例（前 6 行）：

```csv
sample_id,run,genotype,condition,time,replicate,layout,fastq_1,fastq_2,batch
DRR302386,DRR302386,WT,mock,0h,2,PAIRED,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302386/DRR302386_1.fastq.gz,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302386/DRR302386_2.fastq.gz,B1
DRR302389,DRR302389,WT,mock,0.5h,2,PAIRED,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302389/DRR302389_1.fastq.gz,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302389/DRR302389_2.fastq.gz,B1
DRR302397,DRR302397,WT,AvrRpm1,0.5h,1,PAIRED,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302397/DRR302397_1.fastq.gz,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302397/DRR302397_2.fastq.gz,B1
DRR302415,DRR302415,clf,AvrRpm1,0.5h,1,PAIRED,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302415/DRR302415_1.fastq.gz,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302415/DRR302415_2.fastq.gz,B1
DRR302411,DRR302411,clf,mock,3h,3,PAIRED,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302411/DRR302411_1.fastq.gz,http://ftp.sra.ebi.ac.uk/vol1/fastq/DRR302/DRR302411/DRR302411_2.fastq.gz,B1
```

!!! note
    完整的 36 行 samplesheet 见 `artifacts/prjdb11848/metadata/samplesheet.csv`。这份文件是后续第 2–9 章所有分析脚本的输入起点。

---

## 0.5 实战：批量下载与校验

### 准备项目目录

```bash
mkdir -p RNAseq_project/{raw_data,metadata,logs}
cd RNAseq_project
```

### 从 ENA 批量下载

如果你已经从 ENA 导出了包含 `fastq_ftp` 列的 TSV 文件，可以直接提取 URL 进行下载：

```bash
# 从 samplesheet 提取 R1 和 R2 的下载链接
tail -n +2 metadata/samplesheet.csv | while IFS=',' read -r sid run geno cond time rep layout fq1 fq2 batch; do
  echo "Downloading ${run} ..."
  wget -c "${fq1}" -P raw_data/    # -c 支持断点续传
  wget -c "${fq2}" -P raw_data/
done
```

!!! warning "避坑指南"
    ENA 的 FTP 路径并非所有 run 都遵循统一的目录层级规则。对于 10 位数的 run accession（如 DRR3023856），路径中会多一层子目录。**最可靠的做法是直接使用 ENA 导出的真实 URL**，而不是自己拼接路径。

### 完整性校验

下载完成后，必须验证文件完整性。如果 ENA 提供了 md5 校验文件：

```bash
cd raw_data
md5sum -c md5.txt 2>&1 | tee ../logs/md5_check.log

# 检查是否有失败项
if grep -q "FAILED" ../logs/md5_check.log; then
  echo "ERROR: 存在校验失败的文件，请重新下载"
  grep "FAILED" ../logs/md5_check.log
  exit 1
fi
```

!!! warning
    任何 `FAILED` 都必须重新下载。不要带着损坏的 FASTQ 进入下游流程——它不会让 pipeline 报错，但会产生不可预测的结果。

### 本教程的实跑脚本

本教程案例数据的 samplesheet 生成脚本：

```bash
bash scripts/01_prepare_prjdb11848_samplesheet.sh
```

验收：

```bash
wc -l artifacts/prjdb11848/metadata/samplesheet.csv
# 预期: 37（含表头，36 个样本）
```

---

## 0.6 数据组织规范

一个清晰的目录结构是可复现分析的基础 [^noble2009]。我们推荐以下布局：

```text
project_root/
├── data/
│   ├── raw_data/              # 原始 FASTQ，设为只读
│   ├── clean_data/            # fastp 清洗后的数据
│   ├── alignments/            # BAM/CRAM 比对文件
│   ├── quant/                 # Salmon/featureCounts 定量输出
│   └── reference/             # 参考基因组 FASTA/GTF/索引
├── metadata/
│   ├── samplesheet.csv        # 单一真实来源
│   └── design_notes.md        # 分析设计说明
├── scripts/                   # 分析脚本，按步骤编号
│   ├── 01_qc.sh
│   ├── 02_align.sh
│   └── 03_deseq2.R
├── logs/                      # 运行日志
└── results/                   # 分析结果（表格、图表）
```

!!! tip "Recommendation"
    下载完成后立即将 `raw_data/` 设为只读（`chmod -R a-w raw_data/`），防止后续操作意外覆盖原始文件。原始数据是不可再生的——一旦覆盖，只能重新下载。

---

## 0.7 常见翻车点

!!! warning "常见翻车点"
    1. **样本名与文件名不一致**：samplesheet 里写的是 `WT_mock_0h_1`，但 FASTQ 文件名是 `DRR302385_1.fastq.gz`。下游脚本按 sample_id 查找文件时会找不到。解决方案：samplesheet 中同时记录 sample_id 和 run accession，用 run 定位文件。
    2. **技术重复当生物学重复**：同一个 RNA 样本测了两个 lane 不是两个生物学重复。技术重复应在比对后合并，不能作为独立样本进入差异分析，否则会严重高估统计显著性 [^schurch2016]。
    3. **遗漏批次信息**：如果样本分两批上机但 samplesheet 里没有 `batch` 列，后续模型就无法校正批次效应，可能导致大量假阳性。
    4. **只保留处理组数据**：有些人只下载处理组的 FASTQ，不下载对照组。差异分析需要对照组，这是基本前提。
    5. **跳过 MD5 校验**：下载过程中网络中断导致文件截断，但文件名和大小看起来"差不多正常"。不做校验就进入 QC，后面会出现各种随机报错，排查成本极高。

---

## 0.8 Key Takeaways

!!! abstract "本章要点"
    1. **元数据质量决定分析上限**——标签错误、分组缺失等问题不会让 pipeline 报错，但会让结论完全错误，且修复成本极高。
    2. **先定义问题再下载数据**——明确实验设计（比较什么、重复数多少、有哪些混杂因素）是一切分析的起点。
    3. **用一份 samplesheet 驱动全流程**——所有下游脚本从 `samplesheet.csv` 读取样本信息，不在代码中硬编码样本名。
    4. **推荐 ENA 作为首选下载源**——直接提供 FASTQ 链接，元数据导出方便，对批量下载友好。
    5. **下载后必须做 MD5 校验**——文件截断不会让 pipeline 报错，但会产生不可预测的下游结果。
    6. **技术重复 ≠ 生物学重复**——同一 RNA 样本的多个 lane 是技术重复，应在比对后合并，不能作为独立样本。
    7. **目录结构一旦确定就不要改**——路径变动会导致所有脚本中的文件引用失效，是可复现性的大敌。

---

## 0.9 References

[^schurch2016]: Schurch NJ, Schofield P, Gierliński M, et al. (2016) How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22(6):839-851. DOI: [10.1261/rna.053959.115](https://doi.org/10.1261/rna.053959.115)

[^leinonen2011]: Leinonen R, Akhtar R, Birney E, et al. (2011) The European Nucleotide Archive. *Nucleic Acids Research* 39(Database issue):D28-D31. DOI: [10.1093/nar/gkq967](https://doi.org/10.1093/nar/gkq967)

[^noble2009]: Noble WS. (2009) A quick guide to organizing computational biology projects. *PLoS Computational Biology* 5(7):e1000424. DOI: [10.1371/journal.pcbi.1000424](https://doi.org/10.1371/journal.pcbi.1000424)

---

## 0.10 Contributors

### Authors

- *(待补充)*

### Reviewers

- *(待补充)*
