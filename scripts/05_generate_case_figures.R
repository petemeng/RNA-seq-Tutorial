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
