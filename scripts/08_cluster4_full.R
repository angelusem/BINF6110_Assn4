# cluster4 full

# 0. DE run
library(Seurat)
library(dplyr)
library(DESeq2)
library(tibble)
library(ggplot2)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")
Idents(obj) <- "seurat_clusters"

# subset cluster 4 endothelial cells, RM only, Naive vs D02
cl4_rm <- subset(
  obj,
  subset = seurat_clusters == "4" &
    organ_custom == "RM" &
    time %in% c("Naive", "D02")
)

# confirm replicate sizes
cl4_rm_counts <- cl4_rm[[]] %>%
  dplyr::count(mouse_id, time, name = "n_cells") %>%
  dplyr::arrange(time, mouse_id)

print(cl4_rm_counts)

write.csv(
  cl4_rm_counts,
  "results/tables/cluster4_RM_Naive_vs_D02_cell_counts.csv",
  row.names = FALSE
)

# pseudobulk by mouse
cl4_rm$pb_id <- cl4_rm$mouse_id

pb <- AggregateExpression(
  cl4_rm,
  assays = "RNA",
  group.by = "pb_id",
  return.seurat = TRUE
)

pb_counts <- GetAssayData(pb, assay = "RNA", layer = "counts")

# metadata for DESeq2
col_data <- cl4_rm[[]] %>%
  as.data.frame() %>%
  dplyr::select(mouse_id, time, organ_custom, pb_id) %>%
  dplyr::distinct()

# Seurat changed underscores to dashes in pseudobulk column names
col_data$pb_id <- gsub("_", "-", col_data$pb_id)

# order metadata to match pseudobulk columns
col_data <- col_data[match(colnames(pb_counts), col_data$pb_id), , drop = FALSE]

stopifnot(all(colnames(pb_counts) == col_data$pb_id))

rownames(col_data) <- col_data$pb_id
col_data$time <- factor(col_data$time, levels = c("Naive", "D02"))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pb_counts)),
  colData = col_data,
  design = ~ time
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("time", "D02", "Naive")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(padj)

write.csv(
  res,
  "results/tables/cluster4_RM_D02_vs_Naive_DESeq2.csv",
  row.names = FALSE
)

saveRDS(
  dds,
  "results/tables/cluster4_RM_D02_vs_Naive_dds.rds"
)

head(res, 20)
summary(res)

# 1. Making interpretable summary tables:
# full DE table, significant DE table, interpretable significant DE table

# full table is already saved above as:
# results/tables/cluster4_RM_D02_vs_Naive_DESeq2.csv

# significant with effect size
res_sig <- res %>%
  dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

write.csv(
  res_sig,
  "results/tables/cluster4_RM_D02_vs_Naive_DESeq2_sig.csv",
  row.names = FALSE
)

# interpretable = remove ribosomal + mitochondrial genes
res_sig_interpretable <- res_sig %>%
  dplyr::filter(!grepl("^Rpl|^Rps|^mt-", gene))

write.csv(
  res_sig_interpretable,
  "results/tables/cluster4_RM_D02_vs_Naive_DESeq2_sig_interpretable.csv",
  row.names = FALSE
)

# summary counts
sig_all <- res %>% dplyr::filter(!is.na(padj), padj < 0.05)
sig_up <- res %>% dplyr::filter(!is.na(padj), padj < 0.05, log2FoldChange > 1)
sig_down <- res %>% dplyr::filter(!is.na(padj), padj < 0.05, log2FoldChange < -1)

cat("padj < 0.05:", nrow(sig_all), "\n")
cat("padj < 0.05 & |log2FC| > 1:", nrow(res_sig), "\n")
cat("interpretable strong DE genes:", nrow(res_sig_interpretable), "\n")
cat("up in D02:", nrow(sig_up), "\n")
cat("up in Naive:", nrow(sig_down), "\n")

head(res_sig_interpretable, 30)

# 2. making my plots: MA and volcano
library(apeglm)

# shrink log2 fold changes
res_lfc <- lfcShrink(dds, coef = "time_D02_vs_Naive", type = "apeglm")

# save shrunk results
res_shrunk <- as.data.frame(res_lfc) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(padj)

write.csv(
  res_shrunk,
  "results/tables/cluster4_RM_D02_vs_Naive_DESeq2_shrunk.csv",
  row.names = FALSE
)

# MA plot
png("results/figures/cluster4_RM_D02_vs_Naive_MA.png", width = 1800, height = 1400, res = 200)
plotMA(res_lfc, ylim = c(-4, 4), main = "Cluster 4 RM endothelial cells: D02 vs Naive")
dev.off()

# volcano
volcano_df <- res_shrunk %>%
  dplyr::mutate(
    category = dplyr::case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 1 ~ "Up in D02",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Up in Naive",
      TRUE ~ "Not significant"
    )
  )

p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = -log10(padj), color = category)) +
  geom_point(alpha = 0.7, size = 1.2) +
  theme_bw() +
  labs(
    title = "Cluster 4 RM endothelial cells: D02 vs Naive",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value"
  )

print(p_volcano)

ggsave(
  "results/figures/cluster4_RM_D02_vs_Naive_volcano.png",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

# 3. Making my feature plots
# to show why cluster 4 is endothelial and which DE genes changed in the comparison
# A. canonical endothelial markers
FeaturePlot(
  obj,
  features = c("Gpihbp1", "Rbp7", "Sox18", "Ptprb", "Emcn", "Robo4", "Tie1")
)

DotPlot(
  obj,
  features = c("Gpihbp1", "Rbp7", "Sox18", "Ptprb", "Emcn", "Robo4", "Tie1")
) + RotatedAxis()

VlnPlot(
  obj,
  features = c("Gpihbp1", "Rbp7", "Sox18", "Ptprb", "Emcn", "Robo4", "Tie1"),
  group.by = "seurat_clusters",
  pt.size = 0
)
# B. Cluster 4 genes of interest
FeaturePlot(
  obj,
  features = c("Hes1", "Cd200", "Nmt1", "Slc38a2", "Cped1")
)

DotPlot(
  obj,
  features = c("Hes1", "Cd200", "Nmt1", "Slc38a2", "Cped1")
) + RotatedAxis()

VlnPlot(
  obj,
  features = c("Hes1", "Cd200", "Nmt1", "Slc38a2", "Cped1"),
  group.by = "seurat_clusters",
  pt.size = 0
)
# saving 
p_feat1 <- FeaturePlot(
  obj,
  features = c("Gpihbp1", "Rbp7", "Sox18", "Ptprb", "Emcn", "Robo4", "Tie1")
)

ggsave(
  "results/figures/cluster4_endothelial_markers_featureplot.png",
  p_feat1,
  width = 12,
  height = 8,
  dpi = 300
)

p_feat2 <- FeaturePlot(
  obj,
  features = c("Hes1", "Cd200", "Nmt1", "Slc38a2", "Cped1")
)

ggsave(
  "results/figures/cluster4_DE_genes_featureplot.png",
  p_feat2,
  width = 10,
  height = 8,
  dpi = 300
)

# 4. Enrichment analysis
# A. ORA
library(clusterProfiler)
library(org.Mm.eg.db)

sig_up_interpretable <- res_shrunk %>%
  dplyr::filter(!is.na(padj), padj < 0.05, log2FoldChange > 1) %>%
  dplyr::filter(!grepl("^Rpl|^Rps|^mt-", gene))

gene_map <- bitr(
  unique(sig_up_interpretable$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

all_map <- bitr(
  unique(res_shrunk$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ego_bp <- enrichGO(
  gene = unique(gene_map$ENTREZID),
  universe = unique(all_map$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_cc <- enrichGO(
  gene = unique(gene_map$ENTREZID),
  universe = unique(all_map$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_mf <- enrichGO(
  gene = unique(gene_map$ENTREZID),
  universe = unique(all_map$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

write.csv(as.data.frame(ego_bp), "results/tables/cluster4_RM_D02_vs_Naive_ORA_BP.csv", row.names = FALSE)
write.csv(as.data.frame(ego_cc), "results/tables/cluster4_RM_D02_vs_Naive_ORA_CC.csv", row.names = FALSE)
write.csv(as.data.frame(ego_mf), "results/tables/cluster4_RM_D02_vs_Naive_ORA_MF.csv", row.names = FALSE)

nrow(as.data.frame(ego_bp))
nrow(as.data.frame(ego_cc))
nrow(as.data.frame(ego_mf))
head(as.data.frame(ego_bp), 20)
head(as.data.frame(ego_cc), 20)
head(as.data.frame(ego_mf), 20)
# plotting when there are terms in ego
if (nrow(as.data.frame(ego_bp)) > 0) {
  p_ora_bp <- dotplot(ego_bp, showCategory = 20, title = "Cluster 4 ORA BP")
  print(p_ora_bp)
  ggsave("results/figures/cluster4_ORA_BP.png", p_ora_bp, width = 9, height = 7, dpi = 300)
}

if (nrow(as.data.frame(ego_cc)) > 0) {
  p_ora_cc <- dotplot(ego_cc, showCategory = 20, title = "Cluster 4 ORA CC")
  print(p_ora_cc)
  ggsave("results/figures/cluster4_ORA_CC.png", p_ora_cc, width = 9, height = 7, dpi = 300)
}

if (nrow(as.data.frame(ego_mf)) > 0) {
  p_ora_mf <- dotplot(ego_mf, showCategory = 20, title = "Cluster 4 ORA MF")
  print(p_ora_mf)
  ggsave("results/figures/cluster4_ORA_MF.png", p_ora_mf, width = 9, height = 7, dpi = 300)
}
# B. GSEA (using the full-ranked shrunk log2fc list)
rank_df <- res_shrunk %>%
  dplyr::filter(!is.na(log2FoldChange)) %>%
  dplyr::select(gene, log2FoldChange)

gene_map_rank <- bitr(
  unique(rank_df$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

rank_mapped <- rank_df %>%
  dplyr::inner_join(gene_map_rank, by = c("gene" = "SYMBOL")) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::summarise(log2FoldChange = mean(log2FoldChange), .groups = "drop")

gene_list <- rank_mapped$log2FoldChange
names(gene_list) <- rank_mapped$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_bp <- gseGO(
  geneList = gene_list,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 3,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

write.csv(as.data.frame(gsea_bp), "results/tables/cluster4_RM_D02_vs_Naive_GSEA_BP.csv", row.names = FALSE)

head(as.data.frame(gsea_bp), 20)