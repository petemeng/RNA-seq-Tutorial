# Bioinformatics Tutorial — AI 协作技能指南

> **用途**：本文件是 AI 协作撰写生信教程的核心指令。
> 对标风格：[Single-cell best practices](https://www.sc-best-practices.org/) (Theis Lab, Jupyter Book)
> 在执行任何与本项目相关的任务之前，请先完整阅读本文件。

---

## 一、风格定位

### 1.1 对标项目

**sc-best-practices** (theislab/single-cell-best-practices) — 当前单细胞领域最权威的在线最佳实践教程。

这本书的核心特征是：**代码即教程，Benchmark 即推荐依据，实操即学习**。
它不是一本纯理论教科书，而是一本**可执行的、有观点的、持续更新的技术指南**。

### 1.2 与传统教科书的区别

| 维度 | 传统教科书 (如 Alberts) | sc-best-practices 风格 |
|------|------------------------|----------------------|
| 核心载体 | 文字 + 模式图 | **Jupyter Notebook（代码 + 输出 + 解释三位一体）** |
| 知识来源 | 文献综述 + 专家经验 | **Benchmark 驱动 + 社区共识** |
| 推荐态度 | 中立介绍多种方法 | **明确推荐 best practice，说明为什么推荐** |
| 更新频率 | 数年一版 | **持续更新，living resource** |
| 读者参与 | 被动阅读 | **可以下载 notebook 自己跑** |
| 语气 | 严肃学术 | **友好专业，像资深同事在教你** |

### 1.3 核心精神

> **"We recommend X because benchmark Y showed Z."**
> 不是"这里有十个工具你自己选"，而是"根据 benchmark，我们推荐用这个，原因是这些"。

---

## 二、技术架构

### 2.1 构建工具

- **Jupyter Book**：将 Jupyter Notebook (.ipynb) 和 Markdown (.md) 编译为在线书籍
- **GitHub**：版本控制 + 社区贡献 (Issue / PR)
- **GitHub Actions**：自动构建部署
- **Conda/Mamba**：环境管理，每章提供最小环境文件

### 2.2 项目结构

```
Bioinformatics-Tutorial/
│
├── _config.yml                        # Jupyter Book 配置
├── _toc.yml                           # 目录结构
├── requirements.txt                   # Python 依赖
├── .github/workflows/deploy.yml       # 自动部署
├── README.md
│
├── intro.md                           # 前言 (Preamble)
│
├── part1_foundations/                  # 第一部分：基础
│   ├── ch01_xxx.ipynb                 # 每章一个 notebook
│   ├── ch02_xxx.ipynb
│   └── ch03_xxx.md                    # 纯理论章节可用 md
│
├── part2_xxx/                         # 第二部分
│   ├── ch04_xxx.ipynb
│   └── ...
│
├── appendix/
│   ├── glossary.md                    # 术语表
│   └── references.bib                 # 统一参考文献 (BibTeX)
│
├── environments/                      # 每章的 Conda 环境
│   ├── ch01.yml
│   └── ch02.yml
│
└── data/                              # 示例数据（或提供下载链接）
    └── README.md
```

---

## 三、写作风格规范

### 3.1 语言与语气

- **中文为主**，术语首次出现附英文，如：质量控制 (Quality Control, QC)
- 语气：**友好、专业、直接**——像一个资深同事在手把手教你
- 不说"读者可自行了解"，而是**直接教**或给出具体资源链接
- 可以用第一人称复数"我们"（We recommend...）
- 适度口语化，但不随意——**清晰 > 有趣 > 学术**

### 3.2 核心写作原则

1. **代码即教程**
   - 每个分析步骤都有可执行代码
   - 代码不是附录，是正文的核心组成部分
   - 代码块之间穿插解释，说明"为什么这么做"

2. **Benchmark 驱动**
   - 推荐工具时必须给出依据：benchmark 论文、社区共识或实测经验
   - 不回避说"我们推荐 X 而不是 Y"，但要给理由
   - 格式：`我们推荐使用 X {cite}`author2023`，因为在 benchmark Z 中它在 ... 方面表现最优。`

3. **输出即证据**
   - 重要的代码输出（图表、统计量）直接嵌入正文
   - 图表要有图注，解释"这张图说明了什么"
   - 不要只贴图不解释

4. **先动机后操作**
   - 每章开头先讲 Motivation：为什么需要这一步？不做会怎样？
   - 然后再进入具体操作
   - 让读者理解"why"比"how"更重要

5. **明确给出 Key Takeaways**
   - 每章末尾总结核心要点
   - 用简洁的条目列出"你应该记住的几件事"

### 3.3 写作示例

以下是符合本项目风格的写作范例：

````markdown
## 6.1 Motivation

单细胞 RNA-seq 数据集有两个重要特性需要牢记。首先，scRNA-seq 数据存在
大量 dropout（由于 mRNA 捕获效率有限导致的过量零值）。其次，数据校正和
质量控制的空间可能受限，因为技术噪音可能与生物学信号混淆。因此，选择适
合底层数据的预处理方法至关重要，既不能过度校正，也不能移除真实的生物学
效应。

如果在质量控制中过滤掉太多细胞，你可能会丢失稀有的细胞亚群，错过有趣
的生物学发现。反之，如果过于宽松，下游的细胞注释可能会因为低质量细胞的
干扰而变得困难。

## 6.2 Environment setup and data

```python
import numpy as np
import scanpy as sc
import seaborn as sns

sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)
```

我们使用 NeurIPS 2021 单细胞数据整合挑战赛的 10x Multiome 数据集
{cite}`luecken2021`，该数据集包含来自 12 名健康人类供体的骨髓单核细胞。

```python
adata = sc.read_10x_h5(
    filename="filtered_feature_bc_matrix.h5",
    backup_url="https://figshare.com/ndownloader/files/39546196",
)
adata
```

## 6.6 Key Takeaways

1. 质量控制应同时考虑 count depth、detected genes 和 mitochondrial
   fraction 三个指标，不要基于单一指标做决策。
2. 我们推荐使用 MAD-based 自适应阈值而非固定阈值，因为固定阈值难以
   适应不同实验条件和组织类型。
3. 对于 ambient RNA 校正，我们推荐 SoupX，它在多个 benchmark 中表
   现稳定 {cite}`young2020`。
4. Doublet 检测推荐使用 scDblFinder，它在 benchmark 中具有最佳的
   精度-召回平衡 {cite}`xi2021`。
````

---

## 四、每章统一结构

### 4.1 标准章节骨架

```
# N. 章节标题

## N.1 Motivation
  为什么需要这一步？不做会怎样？
  这一步在整个分析流程中处于什么位置？

## N.2 Environment setup and data
  导入依赖包
  加载示例数据（提供公开下载链接）
  简要描述数据集

## N.3 — N.M 核心分析步骤
  每个步骤：
    1. 简短解释原理（2-3 段，不需要太学术）
    2. 代码实现
    3. 输出/图表 + 解释
    4. 如果有方法选择：明确推荐，引用 benchmark

  推荐工具时使用提示框：
  ```{admonition} Recommendation
  我们推荐使用 **X** 进行此步骤。在 benchmark Y 中，
  X 在 ... 指标上优于 Z 和 W {cite}`ref`。
  ```

  注意事项使用警告框：
  ```{admonition} Warning
  :class: warning
  不要在归一化之前进行特征选择，因为...
  ```

## N.M+1 Key Takeaways
  用编号列表总结本章核心要点（5-8 条）。
  每条 1-2 句话，简洁有力。

## N.M+2 References
  本章引用的文献（由 Jupyter Book 的 citation 系统自动生成）

## N.M+3 Contributors
  ### Authors
  - 姓名（单位）
  ### Reviewers
  - 姓名（单位）
```

### 4.2 Admonition 标注框使用规范

| 类型 | 用途 | 示例 |
|------|------|------|
| `{admonition} Recommendation` | 明确推荐某个工具/方法 | "我们推荐使用 scDblFinder..." |
| `{admonition} Note` | 补充说明、技巧 | "这里我们使用 MAD 而非固定阈值..." |
| `{admonition} Warning` `:class: warning` | 常见陷阱、注意事项 | "不要对 log-normalized 数据使用此方法..." |
| `{admonition} Tip` `:class: tip` | 实用技巧 | "可以用 `adata.obs` 快速检查..." |
| `{admonition} See also` | 相关章节交叉引用 | "关于降维的详细讨论，参见第 9 章" |

### 4.3 代码规范

- 使用 **Python** 为主（scanpy/AnnData 生态）
- 必要时提供 **R** 替代方案（Seurat/Bioconductor）
- 代码风格遵循 PEP 8
- 每个代码块上方用 1-2 句话说明这段代码在做什么
- 关键参数用注释解释
- 图表输出后紧接解释文字

```python
# 计算 QC 指标：线粒体基因比例、检测到的基因数、总 count 数
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True,
    percent_top=[20],  # 计算 top 20 基因占比
    log1p=True,         # 对 QC 指标取 log1p
)
```

---

## 五、视觉与排版规范

### 5.1 Jupyter Book 主题配置

使用默认 Jupyter Book 主题或 `sphinx-book-theme`，保持简洁学术风格。

### 5.2 图表规范

- **配色**：使用 scanpy 默认配色或 seaborn 调色板，保持一致
- **DPI**：至少 80（用于 web），发布时可提高到 150
- **图注格式**：`Fig. N.M 图标题。` 后接 1-2 句解释性文字
- **所有图表必须有解释**：不要贴图不说话

### 5.3 交叉引用

- 章节间引用：`参见 [第 9 章：降维](../preprocessing/dimensionality_reduction.html)`
- 术语引用：链接到 Glossary，如 `[UMI](../glossary.html#term-UMI)`
- 文献引用：使用 BibTeX + `{cite}` 语法

---

## 六、AI 协作工作流

### 6.1 写新章节时

用户说"写第 X 章"或提供章节主题时：
1. 按第四节的章节骨架组织
2. 从 Motivation 开始——为什么需要这一步
3. 提供完整可执行的代码（使用公开数据集）
4. 推荐工具时引用 benchmark
5. 结尾写 Key Takeaways
6. 输出为 .ipynb 或 .md 文件

### 6.2 代码审查时

用户说"检查这段代码"或提供 notebook 时：
1. 检查代码是否符合 best practice
2. 指出潜在问题（参数选择、方法选择）
3. 引用 benchmark 支持建议
4. 保持友好专业的语气

### 6.3 工具推荐时

用户问"该用什么工具做 X"时：
1. 给出明确推荐（不要列一堆让用户自己选）
2. 说明推荐理由（benchmark / 社区共识 / 实测）
3. 简要提及替代方案及其局限
4. 提供示例代码

### 6.4 更新内容时

用户说"更新第 X 章"时：
1. 检查是否有新的 benchmark 或工具
2. 更新推荐（如果有更好的工具出现）
3. 标注更新日期和变更内容
4. 保持与其他章节的一致性

---

## 七、写作 Checklist

每章完成后，对照以下检查清单：

- [ ] **Motivation** 清晰解释了"为什么需要这一步"
- [ ] 所有代码块可以**从上到下顺序执行**
- [ ] 使用了**公开可下载**的示例数据
- [ ] 每个工具推荐都有 **benchmark 或文献支持**
- [ ] 所有图表都有**图注和解释文字**
- [ ] **Key Takeaways** 简洁有力（5-8 条）
- [ ] 文献引用使用了 `{cite}` 语法
- [ ] 术语首次出现**附有英文**
- [ ] **Warning** 框标注了常见陷阱
- [ ] **Recommendation** 框标注了核心推荐
- [ ] 代码风格符合 PEP 8
- [ ] 环境依赖已记录在 environment yml 中

---

## 八、与植物免疫 Book 的区别

如果作者同时在写植物免疫教科书和生信教程，注意以下区别：

| 维度 | 植物免疫 Book | 生信教程 |
|------|-------------|---------|
| 风格 | 经典教科书 (Alberts 式) | Jupyter Book 技术指南 (sc-best-practices 式) |
| 核心载体 | Markdown 文字 + 模式图 | **Jupyter Notebook (代码 + 输出)** |
| 构建工具 | MkDocs Material | **Jupyter Book** |
| 推荐态度 | 展示认知演变，不强推 | **明确推荐 best practice** |
| 标注框 | 认知修正（金色）、思路拆解（绿色） | **Recommendation、Warning、Note** |
| 语气 | 学术严谨 | 友好专业，偏实操 |
| 代码 | 无代码 | **代码是核心内容** |

两种风格不要混用。写植物免疫时读 Plant-Immunity-Book-SKILL.md，写生信教程时读本文件。

---

## ⚠️ 注意事项

1. **代码可执行性**：所有代码必须在对应 Conda 环境中可以完整运行，不要写伪代码。
2. **数据可获取性**：使用公开数据集，提供下载链接或 `backup_url`。
3. **Benchmark 时效性**：写入正式章节前应搜索确认是否有更新的 benchmark。
4. **持续更新**：生信工具迭代快，每 6-12 个月应检查推荐是否仍然有效。
5. **社区贡献**：鼓励读者通过 GitHub Issue 提供反馈，每章标注 Contributors。