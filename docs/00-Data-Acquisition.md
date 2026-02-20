# Chapter 0: 数据获取与元数据管理 (Data Acquisition)

> 核心原则：代码可以重写，数据很难重来；而元数据错误会让整个项目不可挽回。

---

## 0.1 先定义问题，再下载数据

在真正执行下载之前，先写清楚三个问题：

1. 你要比较什么（例如 `Mock vs flg22`）？
2. 你的最小生物学重复数是多少（建议每组 `n >= 3`）？
3. 你的主要干扰因素是什么（批次、时间、组织、基因型）？

> Best Practice：先写一份 1 页的“分析假设说明”，再开始做任何命令行操作。

---

## 0.2 数据来源：SRA、ENA、GEO 怎么选

### 为什么推荐 ENA 优先

- ENA 与 SRA 数据同步，但通常提供更直接的 `.fastq.gz` 下载链接。
- 对批量下载更友好，适合脚本化处理。
- 元数据表结构清晰，方便与 `samplesheet` 对接。

### 典型来源关系

- 论文常给 `GSE`（GEO）编号。
- GEO 页面通常链接到 `SRP/PRJNA/SRR`。
- 实际下载常由 ENA 或 SRA 完成。

---

## 0.3 元数据是分析管线的“单一真实来源”

建议用一个 `metadata/samplesheet.csv` 驱动全流程，而不是把样本名写死在脚本里。

### 最小字段模板

| 列名 | 是否必填 | 说明 |
| :--- | :--- | :--- |
| `sample_id` | 是 | 样本唯一标识，不要重复 |
| `run` | 是 | SRR/ERR/DRR 编号 |
| `condition` | 是 | 主要比较因素（如 Mock/Flg22） |
| `replicate` | 是 | 生物学重复编号 |
| `layout` | 是 | `PE` 或 `SE` |
| `fastq_1` | 是 | R1 文件路径或 URL |
| `fastq_2` | PE 必填 | R2 文件路径或 URL |
| `batch` | 建议 | 文库批次/上机批次 |
| `time` | 视设计 | 时间点（如 0h/1h/6h） |

示例：

```csv
sample_id,run,condition,replicate,layout,fastq_1,fastq_2,batch,time
Mock_R1,SRR8694017,Mock,1,PE,raw_data/SRR8694017_1.fastq.gz,raw_data/SRR8694017_2.fastq.gz,B1,0h
Mock_R2,SRR8694018,Mock,2,PE,raw_data/SRR8694018_1.fastq.gz,raw_data/SRR8694018_2.fastq.gz,B1,0h
Flg22_R1,SRR8694019,Flg22,1,PE,raw_data/SRR8694019_1.fastq.gz,raw_data/SRR8694019_2.fastq.gz,B2,1h
Flg22_R2,SRR8694020,Flg22,2,PE,raw_data/SRR8694020_1.fastq.gz,raw_data/SRR8694020_2.fastq.gz,B2,1h
```

---

## 0.4 实战：用 ENA 批量下载并校验

先准备下载目录：

```bash
mkdir -p RNAseq_project/{raw_data,metadata,logs}
cd RNAseq_project
```

如果你已经有 run 列表（每行一个 SRR）：

```bash
# runs.txt 示例
# SRR8694017
# SRR8694018
# SRR8694019
# SRR8694020
```

用 ENA 的 FTP 路径批量下载（示意）：

```bash
while read -r RUN; do
  echo "Downloading ${RUN} ..."
  wget -c "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${RUN:0:6}/${RUN: -1}/${RUN}/${RUN}_1.fastq.gz" -P raw_data/
  wget -c "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${RUN:0:6}/${RUN: -1}/${RUN}/${RUN}_2.fastq.gz" -P raw_data/
done < runs.txt
```

> 避坑指南：不同 SRR 的目录层级并不总是统一，建议优先从 ENA 网页导出真实 URL，再批量下载。

### 完整性校验

```bash
# 使用官方提供的 md5 列表文件（示例名为 md5.txt）
cd raw_data
md5sum -c md5.txt | tee ../logs/md5_check.log
```

任何 `FAILED` 都要重新下载，不要带病进入下游流程。

---

## 0.5 数据组织规范

```text
project_root/
├── data/
│   ├── raw_data/              # 原始 FASTQ，只读
│   ├── clean_data/            # fastp 后的数据
│   ├── alignments/            # BAM/CRAM
│   ├── quant/                 # Salmon/featureCounts 输出
│   └── reference/             # FASTA/GTF/索引
├── metadata/
│   ├── samplesheet.csv
│   └── design_notes.md
├── scripts/
│   ├── 01_qc.sh
│   ├── 02_align.sh
│   └── 03_deseq2.R
├── logs/
└── results/
```

> Best Practice：`raw_data` 设为只读，避免无意覆盖原始文件。

---

## 0.6 常见翻车点

1. 样本名和 FASTQ 文件名不一致，导致下游错配。
2. 把技术重复当生物学重复，导致显著性虚高。
3. 忽略 `batch` 字段，后续模型无法校正批次效应。
4. 只保留处理组，不保留对照组原始数据与元信息。
5. 下载后不做 MD5 校验，后面出现随机报错。

---

## 0.7 本章检查清单

在进入下一章前，确认：

- `samplesheet.csv` 字段完整、无重复 `sample_id`。
- 每个样本都能在磁盘上找到对应 FASTQ。
- 所有 FASTQ 通过 MD5 校验。
- 项目目录结构已经固定，后续不再随意改路径。

如果以上四项有任意一项未完成，不建议继续做 QC 与比对。
