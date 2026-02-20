# RNA-seq Best Practices — AI 写作与技术指南

> **用途**：本文件定义了本教程的编写标准和技术深度，确保内容对标 "Single-cell Best Practices"。

---

## 一、 核心写作哲学

1.  **Why over How**：不仅给代码，更要解释代码背后的生物统计学原理（如：为什么要用负二项分布？为什么要选 STAR 2-pass？）。
2.  **避坑驱动 (Pitfall-Driven)**：每一章必须包含至少一个真实的“生信巨坑”及其解决方案。
3.  **可重复性 (Reproducibility)**：所有环境必须基于 Conda/Mamba，所有软件必须带版本号。
4.  **工程化思维**：推崇 Sample Sheet 驱动的自动化，反对硬编码路径。

## 二、 技术栈标准 (Gold Standard)

| 步骤 | 推荐工具 (Best Practice) | 理由 |
| :--- | :--- | :--- |
| **环境管理** | Mamba / Conda | 版本控制，环境隔离 |
| **质控/清洗** | fastp + MultiQC | 速度快，一步到位，报告精美 |
| **基因组比对** | STAR (2-pass) | 准确度极高，对可变剪切支持最好 |
| **转录本定量** | Salmon | 极速，纠正序列偏好性 (GC/Positional) |
| **差异表达分析** | DESeq2 | 统计模型稳健，收缩估算法 (Shrinkage) 权威 |
| **富集分析** | clusterProfiler | 社区支持最强，数据更新快 |

## 三、 文档排版规范

*   **📘 核心概念**：用于解释底层数学或生物学原理。
*   **⚠️ 避坑指南**：用于标记新手最容易犯错的地方。
*   **💡 Best Practice**：推荐的行业标准操作。
*   **🛠️ 实战演练**：包含具体的可执行代码块。

## 四、 目录结构规范
```text
RNA-seq-Tutorial/
├── docs/                   # 教程 Markdown 文件
├── assets/                 # 图片与资源
└── RNA-seq-Best-Practices-SKILL.md
```
