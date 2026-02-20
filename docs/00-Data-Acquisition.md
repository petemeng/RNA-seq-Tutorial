# Chapter 0: 数据获取与元数据管理 (Data Acquisition)

> **核心原则**：代码是廉价的，数据是昂贵的，而**元数据 (Metadata) 是不可或缺的**。

---

## 1. 数据来源：从 SRA 到 ENA

虽然 NCBI 的 **SRA (Sequence Read Archive)** 是最常用的数据库，但其实际下载体验往往不如欧洲生物信息学研究所的 **ENA (European Nucleotide Archive)**。

### 1.1 为什么优先选择 ENA？
1.  **直接访问 FASTQ**：ENA 提供直接的 `.fastq.gz` 下载链接，无需使用复杂的 `fastq-dump` 转换。
2.  **速度优势**：支持 `Aspera (ascp)` 极速下载。
3.  **元数据清晰**：ENA 的表格索引更适合程序化处理。

---

## 2. 实战演练：下载测试数据

我们将使用一组经典的拟南芥免疫响应数据（例如 PRJNA526071，包含 Mock 和 flg22 处理）。为了教程流畅，我们选取其中 4 个样本。

### 2.1 推荐工具：fasterq-dump (SRA-Toolkit)
如果你必须从 SRA 下载，请放弃过时的 `fastq-dump`，改用多线程的 `fasterq-dump`。

```bash
# 1. 安装 SRA-Toolkit (建议通过 Conda)
# conda install -c bioconda sra-tools

# 2. 下载并转换为 FASTQ (以 SRR8694017 为例)
# --split-files: 自动处理双端测序
# --include-technical: 排除技术序列
fasterq-dump --split-files --include-technical SRR8694017 -e 8 -O raw_data/
```

### 2.2 验证完整性：MD5 Checksum
**绝不要跳过这一步。** 网络波动可能导致文件损坏，产生无法预知的报错。

```bash
# 下载官方提供的 MD5 文件并验证
md5sum -c md5.txt
```

---

## 3. 最佳实践：构建样本表 (Sample Sheet)

**这是 scRNA-best-practice 风格的灵魂：** 不要把样本信息硬编码在脚本里，而是维护一个 `samplesheet.csv`。

| sample_id | group | treatment | time_point | fastq_1 | fastq_2 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| Mock_R1 | Control | Mock | 0h | mock_r1_1.fq.gz | mock_r1_2.fq.gz |
| flg22_R1 | Treated | flg22 | 1h | flg22_r1_1.fq.gz | flg22_r1_2.fq.gz |

### 为什么这样做？
1.  **可读性**：任何人看到这个表都能理解你的实验设计。
2.  **自动化**：下游的 Snakemake 或 Nextflow 流程可以直接读取此表进行批处理。
3.  **一致性**：确保 DGE 分析中的 `Condition` 分组与原始数据严格对应。

---

## 4. 目录结构规范

一个整洁的项目目录是成功的一半：

```text
project_root/
├── data/
│   ├── raw_data/          # 原始 FASTQ (只读)
│   ├── processed_data/    # 比对与定量结果
│   └── reference/         # Genome, GTF, Index
├── metadata/
│   └── samplesheet.csv    # 样本信息表
├── scripts/
│   └── analysis.R         # 分析代码
└── results/               # 图表、表格输出
```
