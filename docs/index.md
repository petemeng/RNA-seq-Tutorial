# RNA-seq Best Practices

欢迎来到 **Bulk RNA-seq 数据分析最佳实践教程**。

本教程聚焦一个目标：构建“可重复、可解释、可发表”的 RNA-seq 分析流程。  
当前版本以**拟南芥（Arabidopsis thaliana）- Mock/flg22**实验体系作为统一主线示例。  
每章都遵循同一教学逻辑：

1. 先讲为什么（统计与生物学原理）。
2. 再讲怎么做（可直接执行的代码）。
3. 最后讲怎么避免翻车（真实项目常见问题）。

---

## 适用人群

- 需要从原始 FASTQ 做到 DEG 与富集分析的研究生/科研人员
- 想把分析流程标准化为可复用管线的课题组
- 需要准备论文 Methods 或复现实验的项目成员

---

## 目录概览

### 基础流程 (The Standard Pipeline)

- [0. 数据获取与元数据管理](00-Data-Acquisition.md)
- [1. 环境与参考基因组准备](01-Environment-and-References.md)
- [2. 质量控制 (QC) 与数据清洗](02-Data-Acquisition-and-QC.md)
- [3. 序列比对与定量](03-Alignment-and-Quantification.md)
- [4. 标准化与探索性分析 (EDA)](04-Normalization-and-Exploratory-Analysis.md)
- [5. 差异表达分析 (DEA)](05-Differential-Expression-Analysis.md)
- [6. 功能富集分析 (Enrichment)](06-Functional-Enrichment-Analysis.md)

### 进阶专题 (Advanced Topics)

- [7. WGCNA 共表达网络分析](07-WGCNA.md)
- [8. 时间序列分析](08-Time-series-Analysis.md)
- [9. 多因素与交互作用设计](09-Multi-factor-Design.md)

---

## 推荐学习路径

1. 第一次接触 bulk RNA-seq：按 `0 -> 6` 顺序完整走一遍。
2. 已有基础分析经验：重点阅读 `4 -> 9`。
3. 准备论文与复现材料：优先检查每章的“翻车点”和“检查清单”。

---

## 使用建议

- 使用 `samplesheet.csv` 驱动所有批处理步骤。
- 每一步都保存日志与版本信息。
- 任何样本剔除或阈值调整都要记录理由。

这套规范能显著提升结果可信度和团队协作效率。

