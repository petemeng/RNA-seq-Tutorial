# RNA-seq Best Practices

本教程采用 `scRNA-best-practice` 风格的信息组织：

1. 先给可执行步骤（Runbook）
2. 再给完整代码（Scripts）
3. 最后给可核验结果（Artifacts + Checksums）

---

## 快速入口

- 端到端实跑步骤：`10-Validated-Case-Study-PRJDB11848.md`
- 代码与结果索引：`13-Code-and-Results-Manifest.md`
- 完整脚本源码：`12-Validated-Case-Study-Full-Scripts.md`
- 故障排查与复现保障：`11-Troubleshooting-and-Reproducibility.md`

---

## 教程结构

### 基础流程

- [0. 数据获取与元数据管理](00-Data-Acquisition.md)
- [1. 环境与参考基因组准备](01-Environment-and-References.md)
- [2. 质量控制 (QC) 与数据清洗](02-Data-Acquisition-and-QC.md)
- [3. 序列比对与定量](03-Alignment-and-Quantification.md)
- [4. 标准化与探索性分析 (EDA)](04-Normalization-and-Exploratory-Analysis.md)
- [5. 差异表达分析 (DEA)](05-Differential-Expression-Analysis.md)
- [6. 功能富集分析 (Enrichment)](06-Functional-Enrichment-Analysis.md)

### 进阶专题

- [7. WGCNA 共表达网络分析](07-WGCNA.md)
- [8. 时间序列分析](08-Time-series-Analysis.md)
- [9. 多因素与交互作用设计](09-Multi-factor-Design.md)

### 实跑与复现

- [10. 可完整复现实战教程（PRJDB11848）](10-Validated-Case-Study-PRJDB11848.md)
- [11. 代码与结果索引](13-Code-and-Results-Manifest.md)
- [12. 实测案例完整脚本源码](12-Validated-Case-Study-Full-Scripts.md)
- [13. 故障排查与复现保障](11-Troubleshooting-and-Reproducibility.md)

---

## 本仓库已上传产物

- 实跑代码：`scripts/`
- 实跑结果快照：`artifacts/prjdb11848/`
- 结果图：`docs/assets/validated_case/`
