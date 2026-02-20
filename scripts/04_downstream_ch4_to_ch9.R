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
