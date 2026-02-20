# 12. 实测案例完整脚本源码（可直接复现）

本章收录 `PRJDB11848` 实测流程的完整脚本源码。你可以按第 10 章步骤直接运行，或在此复制单个脚本使用。

## 12.1 `scripts/01_prepare_prjdb11848_samplesheet.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${1:-validation_run_downstream}"
META_DIR="${BASE_DIR}/metadata"
mkdir -p "${META_DIR}"

echo "[INFO] Fetching PRJDB11848 metadata from ENA ..."
wget -qO "${META_DIR}/prjdb11848.tsv" \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJDB11848&result=read_run&fields=run_accession,sample_title,scientific_name,library_layout,fastq_ftp"

echo "[INFO] Building samplesheet.csv ..."
awk -F"\t" '
BEGIN{
  OFS=",";
  print "sample_id,run,genotype,condition,time,replicate,layout,fastq_1,fastq_2,batch";
}
NR>1{
  run=$1; title=$2; layout=$4; ftp=$5;
  n=split(title,a,"_");
  genotype=a[1]; condition=a[2]; time=a[3]; rep=a[4];
  split(ftp,f,";");
  gsub(/^ftp\.sra\.ebi\.ac\.uk\//,"",f[1]);
  gsub(/^ftp\.sra\.ebi\.ac\.uk\//,"",f[2]);
  fq1="http://ftp.sra.ebi.ac.uk/" f[1];
  fq2="http://ftp.sra.ebi.ac.uk/" f[2];
  sample_id=run;
  print sample_id,run,genotype,condition,time,rep,layout,fq1,fq2,"B1";
}
' "${META_DIR}/prjdb11848.tsv" > "${META_DIR}/samplesheet.csv"

echo "[INFO] samplesheet rows: $(wc -l < "${META_DIR}/samplesheet.csv")"
echo "[INFO] Preview:"
head -n 8 "${META_DIR}/samplesheet.csv"

```

## 12.2 `scripts/02_download_fastq.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${1:-validation_run_downstream}"
SHEET="${BASE_DIR}/metadata/samplesheet.csv"
RAW_DIR="${BASE_DIR}/data/raw_data"
LOG_DIR="${BASE_DIR}/logs"
LOG_FILE="${LOG_DIR}/download_fastq.log"

mkdir -p "${RAW_DIR}" "${LOG_DIR}"

if [[ ! -f "${SHEET}" ]]; then
  echo "[ERROR] Missing samplesheet: ${SHEET}" >&2
  exit 1
fi

: > "${LOG_FILE}"
echo "[INFO] Download start: $(date -Iseconds)" | tee -a "${LOG_FILE}"

awk -F',' 'NR>1 {print $8"\n"$9}' "${SHEET}" | while read -r url; do
  [[ -z "${url}" ]] && continue
  echo "[$(date -Iseconds)] ${url}" | tee -a "${LOG_FILE}"
  ok=0
  for i in $(seq 1 10); do
    if wget -c -nv --tries=1 --timeout=30 "${url}" -P "${RAW_DIR}/" >> "${LOG_FILE}" 2>&1; then
      ok=1
      break
    fi
    echo "retry ${i} failed: ${url}" | tee -a "${LOG_FILE}"
    sleep 2
  done
  if [[ "${ok}" -ne 1 ]]; then
    echo "[ERROR] failed: ${url}" | tee -a "${LOG_FILE}"
    exit 1
  fi
done

echo "[INFO] Download done: $(date -Iseconds)" | tee -a "${LOG_FILE}"
echo "[INFO] FASTQ file count: $(ls "${RAW_DIR}"/*.fastq.gz | wc -l)"

```

## 12.3 `scripts/03_quantify_salmon.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="${1:-validation_run_downstream}"
SHEET="${BASE_DIR}/metadata/samplesheet.csv"
RAW_DIR="${BASE_DIR}/data/raw_data"
QUANT_DIR="${BASE_DIR}/data/quant"
REF_DIR="${BASE_DIR}/data/reference"
LOG_DIR="${BASE_DIR}/logs"

SALMON_BIN="${SALMON_BIN:-salmon}"
THREADS="${THREADS:-8}"
SALMON_INDEX="${SALMON_INDEX:-${REF_DIR}/salmon_index}"

mkdir -p "${QUANT_DIR}" "${LOG_DIR}"

if [[ ! -f "${SHEET}" ]]; then
  echo "[ERROR] Missing samplesheet: ${SHEET}" >&2
  exit 1
fi
if [[ ! -d "${SALMON_INDEX}" ]]; then
  echo "[ERROR] Missing salmon index: ${SALMON_INDEX}" >&2
  exit 1
fi

while IFS=',' read -r sample_id run genotype condition time replicate layout fastq1 fastq2 batch; do
  [[ "${sample_id}" == "sample_id" ]] && continue
  f1="${RAW_DIR}/$(basename "${fastq1}")"
  f2="${RAW_DIR}/$(basename "${fastq2}")"

  "${SALMON_BIN}" quant \
    -i "${SALMON_INDEX}" \
    -l A \
    -1 "${f1}" \
    -2 "${f2}" \
    --validateMappings \
    --gcBias \
    --seqBias \
    -p "${THREADS}" \
    -o "${QUANT_DIR}/${sample_id}" \
    > "${LOG_DIR}/salmon_quant_${sample_id}.log" 2>&1

  echo "salmon done: ${sample_id}"
done < "${SHEET}"

echo "[INFO] Quantification complete. Samples: $(ls "${QUANT_DIR}" | wc -l)"

```

## 12.4 `scripts/04_downstream_ch4_to_ch9.R`

```r
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(DESeq2)
  library(tximport)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(BiocParallel)
  library(pheatmap)
  library(WGCNA)
})

base_dir <- Sys.getenv("BASE_DIR", "validation_run_downstream")
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "ch4"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "ch5"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "ch6"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "ch7"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "ch8"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_dir, "ch9"), recursive = TRUE, showWarnings = FALSE)

samples <- read_csv(file.path(base_dir, "metadata", "samplesheet.csv"), show_col_types = FALSE)
samples$genotype <- factor(samples$genotype, levels = c("WT", "clf"))
samples$condition <- factor(samples$condition, levels = c("mock", "AvrRpm1"))
samples$time <- factor(samples$time, levels = c("0h", "0.5h", "3h"))
samples$batch <- factor(samples$batch)

files <- file.path(base_dir, "data", "quant", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id
tx2gene <- read_tsv(
  file.path(base_dir, "metadata", "tx2gene.tsv"),
  col_names = c("TXNAME", "GENEID"),
  show_col_types = FALSE
)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")

dds <- DESeqDataSetFromTximport(
  txi,
  colData = as.data.frame(samples),
  design = ~ genotype + time + condition + genotype:condition
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)

# Chapter 4
vsd <- vst(dds, blind = FALSE)
pca_df <- plotPCA(vsd, intgroup = c("genotype", "condition", "time"), returnData = TRUE)
write.csv(pca_df, file.path(results_dir, "ch4", "PCA_data.csv"), row.names = FALSE)
sample_dists <- dist(t(assay(vsd)))
sample_mat <- as.matrix(sample_dists)
pdf(file.path(results_dir, "ch4", "sample_distance_heatmap.pdf"), width = 8, height = 6)
pheatmap(
  sample_mat,
  annotation_col = as.data.frame(colData(dds)[, c("genotype", "condition", "time")])
)
dev.off()

# Chapter 5
results_names <- resultsNames(dds)
writeLines(results_names, file.path(results_dir, "ch5", "results_names.txt"))
coef_wt <- grep("^condition_", results_names, value = TRUE)[1]
res_wt <- lfcShrink(dds, coef = coef_wt, type = "apeglm")
res_wt_df <- as.data.frame(res_wt)
res_wt_df <- res_wt_df[order(res_wt_df$padj), ]
write.csv(res_wt_df, file.path(results_dir, "ch5", "DEG_WT_AvrRpm1_vs_mock_shrunk.csv"))
deg_wt <- subset(res_wt_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg_wt, file.path(results_dir, "ch5", "DEG_WT_AvrRpm1_vs_mock_sig.csv"))

res_inter <- results(dds, name = "genotypeclf.conditionAvrRpm1")
res_inter_df <- as.data.frame(res_inter)
res_inter_df <- res_inter_df[order(res_inter_df$padj), ]
write.csv(res_inter_df, file.path(results_dir, "ch5", "DEG_interaction_clf_vs_WT.csv"))
deg_inter <- subset(res_inter_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(deg_inter, file.path(results_dir, "ch5", "DEG_interaction_clf_vs_WT_sig.csv"))

# Chapter 6 ORA + GSEA
sig_gene_tair <- rownames(deg_wt)
universe_tair <- rownames(res_wt_df)

if (length(sig_gene_tair) > 0) {
  id_map <- bitr(sig_gene_tair, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db)
  universe_map <- bitr(universe_tair, fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db)
  gene_entrez <- unique(id_map$ENTREZID)
  universe_entrez <- unique(universe_map$ENTREZID)
  if (length(gene_entrez) >= 10 && length(universe_entrez) > 0) {
    ego <- enrichGO(
      gene = gene_entrez,
      universe = universe_entrez,
      OrgDb = org.At.tair.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    write.csv(as.data.frame(ego), file.path(results_dir, "ch6", "GO_ORA_WT.csv"), row.names = FALSE)
  }
}

gene_rank <- res_wt_df$log2FoldChange
names(gene_rank) <- rownames(res_wt_df)
gene_rank <- gene_rank[!is.na(gene_rank)]
gene_rank <- sort(gene_rank, decreasing = TRUE)
rank_map <- bitr(names(gene_rank), fromType = "TAIR", toType = "ENTREZID", OrgDb = org.At.tair.db)
rank_map <- rank_map[!duplicated(rank_map$ENTREZID), ]
if (nrow(rank_map) > 50) {
  gene_rank_entrez <- gene_rank[rank_map$TAIR]
  names(gene_rank_entrez) <- rank_map$ENTREZID
  gene_rank_entrez <- sort(gene_rank_entrez, decreasing = TRUE)
  gse <- gseGO(
    geneList = gene_rank_entrez,
    OrgDb = org.At.tair.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
  )
  write.csv(as.data.frame(gse), file.path(results_dir, "ch6", "GO_GSEA_WT.csv"), row.names = FALSE)
}

# Chapter 7 WGCNA
allowWGCNAThreads(nThreads = 16)
expr <- assay(vsd)
expr <- expr[rowMeans(counts(dds, normalized = TRUE)) > 10, ]
if (nrow(expr) > 5000) {
  vars <- apply(expr, 1, var)
  keep_genes <- names(sort(vars, decreasing = TRUE))[1:5000]
  expr <- expr[keep_genes, , drop = FALSE]
}
datExpr <- t(expr)
gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
powers <- 1:12
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", corFnc = "bicor", verbose = 5)
write.csv(sft$fitIndices, file.path(results_dir, "ch7", "wgcna_soft_threshold.csv"), row.names = FALSE)
softPower <- 6
net <- blockwiseModules(
  datExpr,
  power = softPower,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE,
  numericLabels = TRUE,
  verbose = 2
)
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
write.csv(data.frame(gene = colnames(datExpr), module = moduleColors), file.path(results_dir, "ch7", "wgcna_modules.csv"), row.names = FALSE)
trait <- as.numeric(samples$condition[match(rownames(datExpr), samples$sample_id)] == "AvrRpm1")
moduleTraitCor <- bicor(MEs, trait, use = "pairwise.complete.obs")
moduleTraitP <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))
write.csv(moduleTraitCor, file.path(results_dir, "ch7", "module_trait_cor.csv"))
write.csv(moduleTraitP, file.path(results_dir, "ch7", "module_trait_p.csv"))

# Chapter 8 time-series LRT (WT only)
wt_samples <- samples[samples$genotype == "WT", ]
wt_files <- file.path(base_dir, "data", "quant", wt_samples$sample_id, "quant.sf")
names(wt_files) <- wt_samples$sample_id
txi_wt <- tximport(wt_files, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM")
dds_wt <- DESeqDataSetFromTximport(
  txi_wt,
  colData = as.data.frame(wt_samples),
  design = ~ condition + time + condition:time
)
keep_wt <- rowSums(counts(dds_wt) >= 10) >= 3
dds_wt <- dds_wt[keep_wt, ]
dds_wt_lrt <- DESeq(dds_wt, test = "LRT", reduced = ~ condition + time)
res_lrt <- results(dds_wt_lrt)
res_lrt_df <- as.data.frame(res_lrt)
res_lrt_df <- res_lrt_df[order(res_lrt_df$padj), ]
write.csv(res_lrt_df, file.path(results_dir, "ch8", "time_series_LRT_WT.csv"))
sig_dynamic <- subset(res_lrt_df, !is.na(padj) & padj < 0.05)
write.csv(sig_dynamic, file.path(results_dir, "ch8", "time_series_LRT_WT_sig.csv"))

# Chapter 9 multifactor
res_treat_in_wt <- results(dds, name = "condition_AvrRpm1_vs_mock")
res_interaction <- results(dds, name = "genotypeclf.conditionAvrRpm1")
res_treat_in_clf <- results(dds, list(c("condition_AvrRpm1_vs_mock", "genotypeclf.conditionAvrRpm1")))
write.csv(as.data.frame(res_treat_in_wt), file.path(results_dir, "ch9", "res_treat_in_wt.csv"))
write.csv(as.data.frame(res_interaction), file.path(results_dir, "ch9", "res_interaction.csv"))
write.csv(as.data.frame(res_treat_in_clf), file.path(results_dir, "ch9", "res_treat_in_clf.csv"))

session <- capture.output(sessionInfo())
writeLines(session, file.path(results_dir, "R_sessionInfo.txt"))

```

## 12.5 `scripts/05_generate_case_figures.R`

```r
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
})

base_dir <- Sys.getenv("BASE_DIR", "validation_run_downstream")
fig_dir <- Sys.getenv("FIG_DIR", "docs/assets/validated_case")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_png <- function(plot_obj, file, width, height, dpi = 150) {
  grDevices::png(filename = file, width = width, height = height, units = "in", res = dpi)
  print(plot_obj)
  grDevices::dev.off()
}

# Chapter 4: PCA
pca <- read_csv(file.path(base_dir, "results/ch4/PCA_data.csv"), show_col_types = FALSE)
if (all(c("PC1", "PC2", "condition", "genotype") %in% colnames(pca))) {
  p <- ggplot(pca, aes(PC1, PC2, color = condition, shape = genotype)) +
    geom_point(size = 3, alpha = 0.9) +
    theme_bw(base_size = 12) +
    labs(title = "PCA (PRJDB11848)", subtitle = "condition + genotype")
  save_png(p, file.path(fig_dir, "ch4_pca.png"), width = 6.5, height = 5)
}

# Chapter 5: WT DEG volcano + MA
deg <- read_csv(file.path(base_dir, "results/ch5/DEG_WT_AvrRpm1_vs_mock_shrunk.csv"), show_col_types = FALSE)
if (all(c("log2FoldChange", "padj", "baseMean") %in% colnames(deg))) {
  deg$neglog10padj <- -log10(pmax(deg$padj, 1e-300))
  deg$sig <- ifelse(!is.na(deg$padj) & deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, "DEG", "NS")

  p_vol <- ggplot(deg, aes(log2FoldChange, neglog10padj, color = sig)) +
    geom_point(alpha = 0.6, size = 1.1) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
    scale_color_manual(values = c("DEG" = "#d73027", "NS" = "grey70")) +
    theme_bw(base_size = 12) +
    labs(
      title = "Volcano Plot: WT AvrRpm1 vs mock",
      x = "log2FoldChange",
      y = "-log10(adjusted p-value)"
    )
  save_png(p_vol, file.path(fig_dir, "ch5_volcano_wt.png"), width = 6.8, height = 5.2)

  p_ma <- ggplot(deg, aes(baseMean, log2FoldChange, color = sig)) +
    geom_point(alpha = 0.6, size = 1.0) +
    scale_x_log10() +
    scale_color_manual(values = c("DEG" = "#d73027", "NS" = "grey70")) +
    theme_bw(base_size = 12) +
    labs(title = "MA Plot: WT AvrRpm1 vs mock", x = "baseMean (log10)", y = "log2FoldChange")
  save_png(p_ma, file.path(fig_dir, "ch5_ma_wt.png"), width = 6.8, height = 5.2)
}

# Chapter 6: ORA / GSEA top terms
ora_file <- file.path(base_dir, "results/ch6/GO_ORA_WT.csv")
if (file.exists(ora_file)) {
  ora <- read_csv(ora_file, show_col_types = FALSE)
  if (nrow(ora) > 0) {
    top_ora <- ora[order(ora$p.adjust), ][1:min(15, nrow(ora)), ]
    top_ora$Description <- factor(top_ora$Description, levels = rev(top_ora$Description))
    p_ora <- ggplot(top_ora, aes(Description, -log10(p.adjust))) +
      geom_col(fill = "#1f78b4") +
      coord_flip() +
      theme_bw(base_size = 11) +
      labs(title = "GO ORA Top Terms", x = "", y = "-log10(adjusted p-value)")
    save_png(p_ora, file.path(fig_dir, "ch6_go_ora_top15.png"), width = 8, height = 6)
  }
}

gsea_file <- file.path(base_dir, "results/ch6/GO_GSEA_WT.csv")
if (file.exists(gsea_file)) {
  gsea <- read_csv(gsea_file, show_col_types = FALSE)
  if (nrow(gsea) > 0) {
    top_gsea <- gsea[order(gsea$p.adjust), ][1:min(15, nrow(gsea)), ]
    top_gsea$Description <- factor(top_gsea$Description, levels = rev(top_gsea$Description))
    p_gsea <- ggplot(top_gsea, aes(Description, NES, fill = NES > 0)) +
      geom_col() +
      coord_flip() +
      scale_fill_manual(values = c("TRUE" = "#e31a1c", "FALSE" = "#1f78b4"), guide = "none") +
      theme_bw(base_size = 11) +
      labs(title = "GO GSEA Top Terms", x = "", y = "NES")
    save_png(p_gsea, file.path(fig_dir, "ch6_go_gsea_top15.png"), width = 8, height = 6)
  }
}

# Chapter 8: LRT p-value distribution
lrt <- read_csv(file.path(base_dir, "results/ch8/time_series_LRT_WT.csv"), show_col_types = FALSE)
if ("padj" %in% colnames(lrt)) {
  p_lrt <- ggplot(lrt, aes(padj)) +
    geom_histogram(bins = 50, fill = "#33a02c", color = "white") +
    theme_bw(base_size = 12) +
    labs(title = "LRT adjusted p-value distribution (WT)", x = "padj", y = "Gene count")
  save_png(p_lrt, file.path(fig_dir, "ch8_lrt_padj_hist.png"), width = 6.5, height = 4.8)
}

# Chapter 9: interaction volcano
inter <- read_csv(file.path(base_dir, "results/ch9/res_interaction.csv"), show_col_types = FALSE)
if (all(c("log2FoldChange", "padj") %in% colnames(inter))) {
  inter$neglog10padj <- -log10(pmax(inter$padj, 1e-300))
  inter$sig <- ifelse(!is.na(inter$padj) & inter$padj < 0.05 & abs(inter$log2FoldChange) > 1, "Sig", "NS")
  p_inter <- ggplot(inter, aes(log2FoldChange, neglog10padj, color = sig)) +
    geom_point(alpha = 0.6, size = 1.0) +
    scale_color_manual(values = c("Sig" = "#6a3d9a", "NS" = "grey70")) +
    theme_bw(base_size = 12) +
    labs(title = "Interaction Volcano (clf:AvrRpm1)", x = "log2FoldChange", y = "-log10(padj)")
  save_png(p_inter, file.path(fig_dir, "ch9_interaction_volcano.png"), width = 6.8, height = 5.2)
}

cat("Figures written to:", normalizePath(fig_dir), "\n")

```
