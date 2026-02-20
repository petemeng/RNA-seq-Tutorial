# Chapter 11: 故障排查与复现保障

---

## 11.1 高优先级排查：物种与参考是否匹配

在本项目实测中，最先遇到的严重问题是：

- 下载到的 SRR（`SRR8694017-4020`）实际物种是 `Gossypium`（棉花）
- 却使用了拟南芥 `TAIR10` 参考
- 直接导致超低比对率与“无 DEG”

建议在下载前先核验 ENA 元数据：

```bash
for r in SRR8694017 SRR8694018 SRR8694019 SRR8694020; do
  wget -qO- "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${r}&result=read_run&fields=run_accession,scientific_name,study_accession,sample_accession"
done
```

---

## 11.2 设计矩阵不满秩（DESeq2）

报错特征：

```text
the model matrix is not full rank
```

常见原因：

- `batch` 与 `condition` 完全共线（例如一个批次只含一个条件）。

处理策略：

1. 优先修正实验设计（最佳）。
2. 若当前数据无法分离批次与处理，降阶到可估计模型（如 `~ condition`）。
3. 在结果报告中显式写明该限制。

---

## 11.3 `lfcShrink(type="apeglm")` 版本兼容

在新版本 DESeq2 中，`apeglm` 推荐配 `coef`，不是 `contrast`：

```r
coef_name <- grep("^condition_", resultsNames(dds), value = TRUE)[1]
res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
```

---

## 11.4 GSEA 并行端口受限

在受限环境中，`BiocParallel` 可能无法开端口，导致 `gseGO()` 失败。  
可强制串行：

```r
gse <- gseGO(
  geneList = gene_rank_entrez,
  OrgDb = org.At.tair.db,
  ont = "BP",
  keyType = "ENTREZID",
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
```

---

## 11.5 WGCNA 运行时间过长

大规模基因直接构网可能耗时非常久。  
教学与实战折中策略：

1. 先过滤低表达基因。
2. 只保留高变基因子集（例如 5k）。
3. 在报告中声明筛选规则，保证可审计性。

---

## 11.6 复现最小清单

建议每次分析都固定输出：

1. `samplesheet.csv`
2. 软件版本清单
3. 关键日志（下载、定量、R）
4. `sessionInfo()`
5. 中间结果矩阵与最终统计表

这样才能做到“可复现”而不是“可运行一次”。
