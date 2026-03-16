# RNA-seq 最佳实践系列

这是一个面向 GitHub Pages 的 bulk RNA-seq 教程站点源码，当前内容基于拟南芥轻度干旱数据集 `DRP010324`，覆盖从数据下载、质控、定量、差异分析到富集分析、WGCNA、可变剪接和发表级图表整理。

在线阅读地址：

- `https://petemeng.github.io/RNA-seq-Tutorial/`

本地生成站点源码：

```bash
python -m tutorial_automation.github_site
```

本地预览：

```bash
pip install -r requirements-docs.txt
mkdocs serve
```

部署方式：

- 推送到 `main` 后由 `.github/workflows/deploy.yml` 自动构建并发布到 GitHub Pages。
