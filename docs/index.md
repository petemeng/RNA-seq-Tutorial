# RNA-seq 数据分析最佳实践

## 这本教程是什么

这是一本 **benchmark 驱动、代码优先、可完整复现** 的 bulk RNA-seq 数据分析指南。我们不会列出十个工具让你自己选——而是根据 benchmark 结果和社区共识，明确推荐每一步的最佳实践，并解释为什么这样推荐。

本教程的风格对标 [Single-cell best practices](https://www.sc-best-practices.org/)（Theis Lab）：每个分析步骤都有可执行代码、真实输出和解释文字三位一体。你可以从头到尾跟着跑一遍，也可以按需查阅某个具体步骤。所有代码都在同一个公开数据集上验证通过，结果可完整复现。

## 谁适合读这本教程

本教程面向需要独立完成 RNA-seq 分析的研究生、博后和生信分析人员。我们假设你具备以下基础：

- **Linux 命令行**：能在终端中导航目录、运行脚本、使用 `wget`/`conda` 等基本工具。
- **R 语言基础**：能读懂 R 代码、安装包、运行脚本。不需要精通，但需要能跟着代码跑。
- **RNA-seq 基本概念**：知道什么是测序 reads、比对、基因表达定量。不需要深入理解统计模型——本教程会在需要时解释。

如果你是完全的初学者，建议先完成一个 Linux 命令行入门教程和 R 语言基础课程，再回来阅读本教程。

## 教程数据集：PRJDB11848

整本教程使用同一个公开数据集贯穿所有章节：**PRJDB11848**，一个拟南芥（*Arabidopsis thaliana*）免疫转录组实验。

实验采用 **2 × 2 × 3 因子设计**：两个基因型（WT 野生型、clf 突变体）× 两个处理（mock 对照、AvrRpm1 效应蛋白）× 三个时间点（0h、0.5h、3h），每组 3 个生物学重复，共 **36 个样本**，全部为双端测序（paired-end）。

我们选择这个数据集的原因是：（1）因子设计足够丰富，能覆盖从简单两组比较到多因素交互、时间序列、WGCNA 等所有进阶分析；（2）样本量适中（36 个），下载和计算成本可控；（3）数据质量良好，能产生清晰的生物学信号（WT 中 AvrRpm1 vs mock 可检出 1,541 个显著 DEG）。

## 如何使用本教程

本教程支持三种使用方式：

- **顺序学习**：从第 0 章开始，按顺序完成基础流程（Ch0–6），再根据兴趣选读进阶专题（Ch7–9）。适合系统学习 RNA-seq 分析的读者。
- **按需查阅**：直接跳到你需要的章节（如差异分析、富集分析）。每章开头的 Motivation 会帮助你理解这一步在整个流程中的位置。
- **直接复现**：如果你只想快速跑通整个流程，可以直接使用实跑与复现部分的脚本和数据，在自己的环境中完整复现所有结果。

## 教程结构

### 基础流程

| 章节 | 内容 |
| :--- | :--- |
| [0. 数据获取与元数据管理](00-Data-Acquisition.md) | 从公共数据库下载数据，构建驱动全流程的 samplesheet |
| [1. 环境与参考基因组准备](01-Environment-and-References.md) | Conda 环境配置、参考基因组和注释文件的获取与索引构建 |
| [2. 质量控制与数据清洗](02-Data-Acquisition-and-QC.md) | FastQC 质量评估、fastp 接头去除与质量过滤、MultiQC 汇总报告 |
| [3. 序列比对与定量](03-Alignment-and-Quantification.md) | HISAT2 比对、featureCounts 基因定量、比对质量评估 |
| [4. 标准化与探索性分析](04-Normalization-and-Exploratory-Analysis.md) | DESeq2 标准化、PCA 降维、样本聚类与离群值检测 |
| [5. 差异表达分析](05-Differential-Expression-Analysis.md) | DESeq2 建模、LFC 收缩、阈值筛选、MA 图与火山图 |
| [6. 功能富集分析](06-Functional-Enrichment-Analysis.md) | GO/KEGG 的 ORA 与 GSEA 分析、结果可视化与解读 |

### 进阶专题

| 章节 | 内容 |
| :--- | :--- |
| [7. WGCNA 共表达网络分析](07-WGCNA.md) | 加权基因共表达网络构建、模块鉴定与性状关联 |
| [8. 时间序列分析](08-Time-series-Analysis.md) | DESeq2 LRT 检验时间效应、时序表达模式聚类 |
| [9. 多因素与交互作用设计](09-Multi-factor-Design.md) | 交互项建模、基因型 × 处理的差异响应分析 |

### 实跑与复现

| 章节 | 内容 |
| :--- | :--- |
| [10. 端到端实战复现（PRJDB11848）](10-Validated-Case-Study-PRJDB11848.md) | 从原始数据到最终结果的完整复现步骤 |
| [11. 故障排查与复现保障](11-Troubleshooting-and-Reproducibility.md) | 常见问题解决方案、环境锁定、结果校验 |
| [12. 实测案例完整脚本源码](12-Validated-Case-Study-Full-Scripts.md) | 所有分析脚本的完整源代码 |
| [13. 代码与结果索引](13-Code-and-Results-Manifest.md) | 脚本、中间文件、最终结果的完整清单 |

## 致谢与贡献

本教程是一个持续更新的开源项目。如果你发现错误、有改进建议或想贡献新内容，欢迎通过以下方式参与：

- **报告问题**：在 [GitHub Issues](https://github.com/petemeng/RNA-seq-Tutorial/issues) 中提交 bug 或建议。
- **贡献内容**：Fork 本仓库，修改后提交 Pull Request。
- **讨论交流**：在 Issue 中发起讨论，分享你的使用经验或分析心得。

生信工具迭代很快，我们会定期检查并更新推荐的工具和方法。如果你发现某个推荐已经过时，请告诉我们。
