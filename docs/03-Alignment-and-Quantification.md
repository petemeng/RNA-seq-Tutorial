# 📖 Bulk RNA-seq 数据分析最佳实践 (Part 3)

## 第三章：序列比对 —— “ATCG”的归位 (Alignment)

在质量控制之后，我们需要将数亿条短序列 (Reads) 比对到参考基因组上。目前，业界有两个最主流的选择：**STAR**（传统比对派）和 **Salmon**（免比对定量派）。

### 3.1 内存之王：STAR (Spliced Transcripts Alignment to a Reference)

如果你要研究**可变剪切 (Alternative Splicing)**、**基因融合 (Gene Fusion)** 或**未注释的新转录本**，STAR 是无可争议的“金标准”。

> 📘 **核心概念：剪接位点感知的比对 (Spliced-aware Alignment)**
> 与 DNA 比对不同，RNA Reads 往往跨越外显子 (Exon)。STAR 能够自动识别剪接位点 (Splice Junctions)，即使一条 Read 的前一半在 1 号外显子，后一半在几万个碱基之外的 2 号外显子，它也能精准定位。

#### ⚠️ 避坑指南：STAR 的“吞金”属性
*   **现象**：STAR 在比对人类或小鼠基因组时，通常需要 **30GB - 40GB** 内存。内存不足会导致进程被 Linux 内核直接杀死 (Out of Memory, OOM)。
*   **解决方案**：如果你的服务器内存有限，可以使用 `--genomeSAindexNbases` 参数减小索引大小，或者在比对时使用共享内存模式。

#### 💡 Best Practice：开启 2-pass 模式
在单次比对中，STAR 可能漏掉一些非典型的剪接位点。
**最佳实践：** 开启 `--twopassMode Basic`。它会先进行一轮比对寻找剪接位点，然后将这些位点加入索引，再进行第二轮比对。这能显著提升剪接点识别的灵敏度和准确度。

```bash
# 构建索引 (以人类 GRCh38 为例)
STAR --runThreadN 8 --runMode genomeGenerate \
     --genomeDir ./star_index \
     --genomeFastaFiles ./reference/GRCh38.fa \
     --sjdbGTFfile ./reference/gencode.v45.annotation.gtf \
     --sjdbOverhang 149

# 比对实战 (开启 2-pass)
STAR --runThreadN 8 --genomeDir ./star_index \
     --readFilesIn cleandata/SRR1039508_1_clean.fq.gz cleandata/SRR1039508_2_clean.fq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix ./star_out/SRR1039508_ \
     --outSAMtype BAM SortedByCoordinate \
     --twopassMode Basic \
     --quantMode GeneCounts 
```

---

### 3.2 速度之王：Salmon (Pseudo-alignment)

如果你的目标只是**差异表达分析 (DGE)**，那么 **Salmon** 是目前的“新宠”。

> 📘 **核心概念：准比对 (Pseudo-alignment) 与期望最大化 (EM) 算法**
> Salmon 不逐个碱基寻找 Reads 的精确位置，而是寻找 Reads 属于哪个转录本的可能性。它使用 EM 算法来解决**多重比对 (Multi-mapping)** 问题。

#### 💡 Why Salmon is Better for DGE?
1.  **极速**：比 STAR 快 10 倍以上，且内存消耗极低。
2.  **偏差校正**：Salmon 能自动校正 **GC 偏好性** 和 **位置偏好性**，这在统计学上让表达量估计更接近真实值。

```bash
# 构建转录组索引
salmon index -t ./reference/transcripts.fa -i salmon_index -k 31

# 定量实战 (开启 GC 与序列偏好校正)
salmon quant -i salmon_index -l A \
     -1 cleandata/SRR1039508_1_clean.fq.gz -2 cleandata/SRR1039508_2_clean.fq.gz \
     --gcBias --seqBias \
     -o salmon_quant/SRR1039508_quant
```

---

### 3.3 🚫 最深的一道坎：链特异性 (Strand-specificity)

这是导致无数生信新手分析彻底失败的“头号杀手”。

> ⚠️ **避坑指南 (Pitfall)**
> 现代建库方法通常是有**链特异性**的。这意味着 Reads 是有方向的。
> 如果选错了参数（如：明明是反向链 Reverse-stranded 却选了正向），**全样本表达量可能接近为 0**。

#### 如何确定链特异性？
1.  **问测序公司**：这是最保险的方法。
2.  **自动探测**：Salmon 使用 `-l A` (Automatic) 参数可以自动探测；或者使用 `check_strandness.py` 脚本预览。
