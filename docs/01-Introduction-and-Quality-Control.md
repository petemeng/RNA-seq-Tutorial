# RNA-seq Best Practices: 从原始序列到生物学发现

> **致读者**：本教程不只是告诉你运行哪个工具，而是解释**为什么要这么做**。我们遵循主流学术界（如 NF-Core, Bioconductor）公认性 Gold Standard。

---

## 1. 实验设计：统计效能的起点 (Experimental Design)

在任何测序开始前，最核心的不是技术，而是**统计效能 (Statistical Power)**。

### 1.1 生物学重复 (Biological Replicates)
*   **原则**：绝对不要做没有重复的测序。
*   **最佳实践**：至少 **n=3**。对于临床样本或异质性强的植物野外样本，建议 **n≥6**。
*   **逻辑**：差异表达分析（如 DESeq2）本质上是比较“组间差异”与“组内方差”。没有重复，我们就无法估计组内方差，任何观察到的差异都可能是随机噪声。

### 1.2 测序深度与读长 (Sequencing Depth & Read Length)
*   **差异表达 (DGE)**：10M - 20M reads/sample 即可，单端 (Single-end) 50-75bp 足够。
*   **可变剪接/新转录本发现**：建议 >50M reads，且必须使用双端 (Paired-end) 100-150bp。

---

## 2. 原始数据质量控制：不仅仅是 FastQC

### 2.1 核心评价指标 (Key Metrics)

| 指标 | 阈值 (Best Practice) | 说明 |
| :--- | :--- | :--- |
| **Per Base Quality** | > Q30 (99.9% 准确度) | 碱基识别准确率。后期下降是测序反应化学衰减的正常现象。 |
| **Adapter Content** | < 0.1% | 接头残留会干扰比对，必须在比对前去除。 |
| **GC Content** | 接近物种分布 (如拟南芥 ~36%) | 异常峰可能意味着污染或过度的 PCR 扩增。 |
| **Sequence Duplication** | 视实验而定 | 高表达基因会导致自然重复，但极高的重复率可能暗示文库复杂度不足。 |

### 2.2 工具链：FastQC + MultiQC
不要一个一个看报告。使用 `MultiQC` 将所有样本的质控结果汇总到一张交互式图表中。

```bash
# 典型的质控脚本
fastqc -t 8 -o qc_results/ data/*.fastq.gz
multiqc qc_results/ -o multiqc_report/
```

### 2.3 深度洞察：为什么关注 GC 偏好性？
GC 含量偏好性 (GC-bias) 是导致样本间非生物学差异的主要原因之一。在下游标准化中，DESeq2 等工具会尝试校正文库大小，但对于严重的 GC 偏好，可能需要使用 `cqn` (Conditional Quantile Normalization) 进行专门处理。

---

## 3. 数据修剪 (Trimming) 的权衡

### 3.1 最佳工具建议：fastp
传统的 `Trimmomatic` 参数复杂且速度慢。目前学术界更倾向于使用 **`fastp`**，它能同时完成：
1. 自动接头识别 (Auto-adapter detection)
2. 低质量碱基切除 (Sliding window trimming)
3. 简单的质控统计

### 3.2 什么时候不建议过度修剪？
如果你使用的是 **Pseudo-alignment** 工具（如 Salmon 或 Kallisto），适度的修剪（仅去接头和极低质量碱基）通常比剧烈的硬切除更好。过度修剪缩短了 Read 长度，会降低转录本定量时解决多重比对 (Multi-mapping) 的能力。
