# 术语表 (Glossary)

本术语表收录教程中反复出现的核心概念，按英文字母顺序排列。每个术语格式为：**英文缩写 / 英文全称 (中文)** — 简要定义。

## A

**Admonition (标注框)** — MkDocs Material 提供的高亮提示框组件，用于标注推荐、警告、技巧等信息。本教程使用 `note`、`warning`、`tip`、`abstract` 四种类型。

**apeglm** — Approximate Posterior Estimation for Generalized Linear Models 的缩写，DESeq2 推荐的 log2 fold change 收缩方法，相比 `normal` 和 `ashr` 在保留大效应量信号方面表现更优。

## B

**Batch effect (批次效应)** — 由非生物学因素（如不同实验日期、不同操作员、不同测序 lane）引入的系统性技术变异。需要在实验设计和统计建模中加以控制。

**BH correction (Benjamini-Hochberg 校正)** — 一种控制 FDR 的多重检验校正方法，是 RNA-seq 差异分析中最常用的 p 值校正策略。

## C

**Composition bias (组成偏差)** — 当少数高表达基因在不同条件间发生剧烈变化时，基于总 count 的归一化会系统性地扭曲其余基因的表达量估计。DESeq2 的 median-of-ratios 方法专门应对此问题。

**Count matrix (计数矩阵)** — RNA-seq 定量的核心数据结构，行为基因、列为样本，每个元素表示该基因在该样本中被检测到的 read/fragment 数量。

## D

**DEG (Differentially Expressed Gene, 差异表达基因)** — 在两个或多个实验条件之间表达水平存在统计学显著差异的基因，通常以 padj < 0.05 且 |log2FC| > 1 为筛选标准。

**DESeq2** — 基于负二项分布的 RNA-seq 差异表达分析 R 包，是目前 bulk RNA-seq 领域最广泛使用的工具之一。

**Dispersion (离散度)** — 负二项分布中描述方差与均值关系的参数。DESeq2 通过 shrinkage estimation 在全基因组水平借力估计每个基因的离散度，提高低重复数下的统计稳健性。

## E

**ENA (European Nucleotide Archive, 欧洲核酸档案库)** — EMBL-EBI 维护的公共测序数据库，与 NCBI SRA 数据同步，提供更便捷的 FASTQ 直接下载链接。

## F

**FASTQ** — 存储测序 read 及其碱基质量分数的标准文本格式，每条 read 占 4 行。

**FDR (False Discovery Rate, 假发现率)** — 在所有被判定为"显著"的结果中，假阳性所占的期望比例。RNA-seq 分析中通常控制 FDR < 0.05。

**FPKM (Fragments Per Kilobase of transcript per Million mapped fragments)** — 一种归一化指标，同时校正基因长度和测序深度。适用于同一样本内基因间比较，但不适合样本间比较。

## G

**GEO (Gene Expression Omnibus)** — NCBI 维护的基因表达数据公共仓库，研究者通常在论文中提供 GEO 编号（如 GSE123456）作为数据入口。

**GO (Gene Ontology, 基因本体论)** — 描述基因功能的标准化词汇体系，分为三个子本体：Biological Process (BP)、Molecular Function (MF)、Cellular Component (CC)。

**GSEA (Gene Set Enrichment Analysis, 基因集富集分析)** — 一种不依赖硬阈值的富集分析方法，使用全部基因的排序列表（通常按 log2FC 排序）检测基因集是否在排序的两端富集。

## H

**Hub gene (枢纽基因)** — 在 WGCNA 共表达网络中，与同模块内其他基因具有高连接度的基因，通常被认为是该模块生物学功能的关键调控因子。

## I

**Interaction effect (交互效应)** — 在多因素实验设计中，一个因素的效应随另一个因素水平的变化而变化的现象。例如，基因型对处理响应的调节作用。在 DESeq2 中通过 `genotype:condition` 交互项建模。

## K

**KEGG (Kyoto Encyclopedia of Genes and Genomes)** — 整合基因组、化学和系统功能信息的数据库，常用于通路富集分析。

## L

**LFC / log2FC (log2 Fold Change)** — 两个条件间基因表达量比值的以 2 为底的对数。LFC = 1 表示表达量翻倍，LFC = -1 表示减半。

**LRT (Likelihood Ratio Test, 似然比检验)** — 通过比较完整模型与简化模型的拟合优度来检验特定因素效应的统计方法。在时间序列分析中，LRT 用于检测基因表达是否随时间发生动态变化。

## M

**Main effect (主效应)** — 在多因素实验设计中，单个因素独立于其他因素对响应变量的影响。

**Mapping rate (比对率)** — 成功比对到参考基因组或转录组的 read 占总 read 数的比例，是评估测序质量和参考选择是否正确的关键指标。

**Module eigengene (模块特征基因)** — WGCNA 中，代表一个共表达模块整体表达模式的第一主成分，用于将模块与表型性状进行关联分析。

## N

**NB distribution (Negative Binomial distribution, 负二项分布)** — RNA-seq count 数据的标准统计模型，能够处理生物学重复间的过离散（overdispersion）现象，即方差大于均值的情况。

## O

**ORA (Over-Representation Analysis, 过表示分析)** — 基于超几何检验的富集分析方法，检验 DEG 列表中某个功能类别的基因是否显著多于随机期望。需要预先设定 DEG 筛选阈值。

## P

**padj (adjusted p-value, 校正后 p 值)** — 经过多重检验校正后的 p 值，用于控制 FDR。在 DESeq2 中默认使用 BH 方法校正。

**PCA (Principal Component Analysis, 主成分分析)** — 一种降维方法，将高维基因表达数据投影到少数几个主成分上，用于可视化样本间的整体关系、检测批次效应和离群样本。

## R

**rlog (regularized log transformation, 正则化对数变换)** — DESeq2 提供的方差稳定化变换方法之一，对低 count 基因的变换更稳健，但计算速度慢于 VST，适用于小样本量（< 30）。

## S

**Salmon** — 基于 quasi-mapping 的快速转录本定量工具，无需传统比对步骤即可直接从 FASTQ 估计转录本丰度，支持 GC bias 和 sequence bias 校正。

**Size factor (大小因子)** — DESeq2 用于校正样本间测序深度差异的归一化因子，通过 median-of-ratios 方法计算，对组成偏差具有鲁棒性。

**Soft threshold (软阈值)** — WGCNA 中用于将相关性矩阵转换为邻接矩阵的幂次参数，选择标准是使网络拓扑近似无标度（scale-free topology）。

**SRA (Sequence Read Archive)** — NCBI 维护的公共测序数据档案库，存储原始测序数据（`.sra` 格式），需使用 `fasterq-dump` 转换为 FASTQ。

## T

**TPM (Transcripts Per Million)** — 一种归一化指标，先校正基因长度再校正测序深度。所有基因的 TPM 之和恒为 100 万，便于样本间比较相对丰度。

**tximport** — 将 Salmon/kallisto 等转录本级定量结果汇总为基因级 count 矩阵的 R 包，是连接定量工具与 DESeq2 的标准桥梁。

**tx2gene** — 转录本 ID 到基因 ID 的映射表，tximport 用它将转录本级丰度汇总为基因级计数。

## V

**VST (Variance Stabilizing Transformation, 方差稳定化变换)** — DESeq2 提供的变换方法，使不同表达水平的基因具有近似相同的方差，适用于 PCA、聚类、热图等需要同方差假设的下游分析。计算速度快于 rlog，推荐用于大样本量。

## W

**WGCNA (Weighted Gene Co-expression Network Analysis, 加权基因共表达网络分析)** — 一种系统生物学方法，通过构建基因间的加权共表达网络，识别高度协同表达的基因模块，并将模块与外部表型性状关联，用于发现功能相关的基因群和枢纽基因。
