# scripts/08_cluster1_macrophage_RM_D02_vs_Naive_DE.R

# determining which comparison to make in cluster 1
library(Seurat)
library(dplyr)

obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")
Idents(obj) <- "seurat_clusters"

mac1 <- subset(obj, subset = seurat_clusters == "1")

# broad checks
table(mac1$organ_custom)
table(mac1$time)
table(mac1$organ_custom, mac1$time)

# sample-level counts for pseudobulk planning
mac1_sample_counts <- mac1[[]] %>%
  count(organ_custom, time, biosample_id, mouse_id, name = "n_cells") %>%
  arrange(organ_custom, time, biosample_id)

mac1_sample_counts

write.csv(
  mac1_sample_counts,
  "results/tables/cluster1_macrophage_sample_counts.csv",
  row.names = FALSE
)

.rs.restartR()
# clean environment after 
# ============================================================= #

library(Seurat)
library(dplyr)
library(DESeq2)
library(tibble)

obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")
Idents(obj) <- "seurat_clusters"

# cluster 1 macrophages, RM only, Naive vs D02
mac_rm <- subset(
  obj,
  subset = seurat_clusters == "1" &
    organ_custom == "RM" &
    time %in% c("Naive", "D02")
)

# confirm replicate sizes
mac_rm_counts <- mac_rm[[]] %>%
  dplyr::count(mouse_id, time, name = "n_cells") %>%
  dplyr::arrange(time, mouse_id)

mac_rm_counts
write.csv(mac_rm_counts, "results/tables/cluster1_RM_Naive_vs_D02_cell_counts.csv", row.names = FALSE)

# use mouse_id as pseudobulk sample ID
mac_rm$pb_id <- mac_rm$mouse_id

# aggregate raw counts by mouse
pb <- AggregateExpression(
  mac_rm,
  assays = "RNA",
  group.by = "pb_id",
  return.seurat = TRUE
)

# extract pseudobulk count matrix
pb_counts <- GetAssayData(pb, assay = "RNA", layer = "counts")

# build sample metadata
col_data <- mac_rm[[]] %>%
  as_tibble() %>%
  select(pb_id, mouse_id, time, organ_custom) %>%
  distinct() %>%
  arrange(match(pb_id, colnames(pb_counts)))

col_data$time <- factor(col_data$time, levels = c("Naive", "D02"))
rownames(col_data) <- col_data$pb_id

# column names' underscores replaced with dashes: fixing this 
col_data <- mac_rm[[]] %>%
  as_tibble() %>%
  dplyr::select(pb_id, mouse_id, time, organ_custom) %>%
  dplyr::distinct()

# match Seurat's underscore-to-dash conversion
col_data <- col_data %>%
  dplyr::mutate(pb_id = gsub("_", "-", pb_id))

# order metadata to match pseudobulk columns
col_data <- col_data %>%
  dplyr::slice(match(colnames(pb_counts), pb_id))

# convert to plain data.frame before rownames
col_data <- as.data.frame(col_data)
rownames(col_data) <- col_data$pb_id

# then

# reorder count matrix to match metadata
pb_counts <- pb_counts[, rownames(col_data)]

# ensuring that time is a factor before building dds to keep the reference level explicit
col_data$time <- factor(col_data$time, levels = c("Naive", "D02"))

# DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pb_counts)),
  colData = col_data,
  design = ~ time
)

# prefilter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# run DE
dds <- DESeq(dds)

res <- results(dds, contrast = c("time", "D02", "Naive")) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(padj)

write.csv(res, "results/tables/cluster1_RM_D02_vs_Naive_DESeq2.csv", row.names = FALSE)

# save full DESeq object too
saveRDS(dds, "results/tables/cluster1_RM_D02_vs_Naive_dds.rds")

# inspect top hits
head(res, 20)
summary(res)
# one quick DE plot-- before making the interpretable table
library(ggplot2)
library(dplyr)

res <- read.csv("results/tables/cluster1_RM_D02_vs_Naive_DESeq2.csv")

res$significance <- "NS"
res$significance[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 1] <- "Up in D02"
res$significance[!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -1] <- "Up in Naive"

p <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.7, size = 1.5) +
  theme_minimal() +
  labs(
    title = "Cluster 1 macrophages, RM: D02 vs Naive",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  )
p

ggsave("results/figures/cluster1_RM_D02_vs_Naive_volcano.png", p, width = 8, height = 6, dpi = 300)


# Making a cleaner DE table

# shrinking LFCs for presentation
#BiocManager::install("apeglm")
res_shrunk <- lfcShrink(dds, coef = "time_D02_vs_Naive", type = "apeglm") %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(padj)

write.csv(
  res_shrunk,
  "results/tables/cluster1_RM_D02_vs_Naive_DESeq2_shrunk.csv",
  row.names = FALSE
)

# significance + effect size filter
res_sig <- res_shrunk %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

write.csv(
  res_sig,
  "results/tables/cluster1_RM_D02_vs_Naive_DESeq2_sig.csv",
  row.names = FALSE
)

# for interpretation purposes only: removing ribosomal and mitochondrial transcripts to highlight more condition-specific genes

res_sig_interpretable <- res_sig %>%
  filter(
    !grepl("^Rpl|^Rps|^mt-", gene)
  )
# what exactly is this removing? 
res_sig %>%
  filter(grepl("^Rpl|^Rps|^mt-", gene)) %>%
  head(20)
# what remains?
head(res_sig_interpretable, 20)

head(res_sig_interpretable, 30)


write.csv(
  res_sig_interpretable,
  "results/tables/cluster1_RM_D02_vs_Naive_DESeq2_sig_interpretable.csv",
  row.names = FALSE
)

# checking how many genes I actually have: getting simple counts
nrow(res_sig)
sum(res_sig$log2FoldChange > 0)
sum(res_sig$log2FoldChange < 0)

nrow(res_sig_interpretable)
sum(res_sig_interpretable$log2FoldChange > 0)
sum(res_sig_interpretable$log2FoldChange < 0)

# making my main DE figure
# MA plot
#plotMA(as.matrix = FALSE, res_shrunk, ylim = c(-4, 4)) # doesnt work

res_lfc <- lfcShrink(dds, coef = "time_D02_vs_Naive", type = "apeglm")

plotMA(res_lfc, ylim = c(-4, 4))

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

png(
  filename = "results/figures/cluster1_RM_D02_vs_Naive_MA.png",
  width = 1800,
  height = 1400,
  res = 200
)
plotMA(
  res_lfc,
  ylim = c(-4, 4),
  main = "Cluster 1 RM macrophages: D02 vs Naive"
)
dev.off()

# volcano plot
library(ggplot2)

volcano_df <- res_shrunk %>%
  mutate(
    category = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up in D02",
      padj < 0.05 & log2FoldChange < -1 ~ "Up in Naive",
      TRUE ~ "Not significant"
    )
  )

p_volcano_cl1 <- ggplot(
  volcano_df,
  aes(x = log2FoldChange, y = -log10(padj), color = category)
) +
  geom_point(alpha = 0.7, size = 1.2) +
  theme_bw() +
  labs(
    title = "Cluster 1 RM macrophages: D02 vs Naive",
    x = "Shrunken log2 fold change",
    y = "-log10 adjusted p-value"
  )

p_volcano_cl1

ggsave(
  filename = "results/figures/cluster1_RM_D02_vs_Naive_volcano_shrunk.png",
  plot = p_volcano_cl1,
  width = 8,
  height = 6,
  dpi = 300
)

# Checking whether top genes are consistently different across all 3 mice per group (so that I do not rely only on the ranked table) ==
norm_counts <- counts(dds, normalized = TRUE) %>%
  as.data.frame() %>%
  rownames_to_column("gene")

genes_to_check <- c("Plek", "Btg2", "Ier2", "Nmt1", "Cbr2", "Mup4")

plot_df <- norm_counts %>%
  dplyr::filter(gene %in% genes_to_check) %>%
  tidyr::pivot_longer(-gene, names_to = "pb_id", values_to = "norm_count") %>%
  dplyr::left_join(
    col_data[, c("pb_id", "mouse_id", "time", "organ_custom")],
    by = "pb_id"
  )

ggplot(plot_df, aes(x = time, y = norm_count, color = mouse_id)) +
  geom_point(size = 3) +
  geom_line(aes(group = mouse_id)) +
  facet_wrap(~ gene, scales = "free_y") +
  theme_bw()
# =================================================================

# validating whether cluster 1 is really clean enough to report
FeaturePlot(obj, features = c(
  "Adgre1", "Csf1r", "Lyz2", "Tyrobp", "Ctss", "C1qa",
  "Plek", "Hk2", "Gpr183",
  "Omp", "Gng13", "S100a5", "Stoml3"
))
# checking them in cluster 1 only 
VlnPlot(obj, features = c(
  "Adgre1", "Csf1r", "Lyz2", "Tyrobp", "Ctss", "C1qa",
  "Plek", "Hk2", "Gpr183",
  "Omp", "Gng13", "S100a5", "Stoml3"
), group.by = "seurat_clusters", pt.size = 0)
#
DimPlot(obj, label = TRUE)
DotPlot(
  obj,
  features = c("Adgre1","Csf1r","Lyz2","Tyrobp","Ctss","C1qa","Omp","Gng13","S100a5","Stoml3")
) + RotatedAxis()
# =======================================================================

# ORA
library(clusterProfiler)
library(org.Mm.eg.db)

sig_up <- res_sig_interpretable %>%
  filter(log2FoldChange > 1)

gene_map <- bitr(
  sig_up$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

all_map <- bitr(
  res_shrunk$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ego <- enrichGO(
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

dotplot(ego, showCategory = 20) # plot has zero enriched terms! 

# attempt 2: softer threshold and lower mingssize
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

sig_up <- res_shrunk %>%
  filter(!is.na(padj), padj < 0.05, log2FoldChange > 0.5)

gene_map <- bitr(
  unique(sig_up$gene),
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

ego <- enrichGO(
  gene = unique(gene_map$ENTREZID),
  universe = unique(all_map$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  minGSSize = 3,
  readable = TRUE
)

nrow(as.data.frame(ego))
head(as.data.frame(ego), 20)

dotplot(ego, showCategory = 20) # still zero

# attempt 3
library(clusterProfiler)
library(org.Mm.eg.db)

sig_up <- res_sig_interpretable %>%
  filter(log2FoldChange > 1)

gene_map <- bitr(
  sig_up$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

all_map <- bitr(
  res_shrunk$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ego <- enrichGO(
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

dotplot(ego, showCategory = 20) # 1 single gene

# attempt 4
library(clusterProfiler)
library(org.Mm.eg.db)

sig_up <- res_sig_interpretable %>%
  filter(log2FoldChange > 1)

gene_map <- bitr(
  sig_up$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

all_map <- bitr(
  res_shrunk$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ego <- enrichGO(
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

dotplot(ego, showCategory = 20) # looks better!
# inspecting the results
as.data.frame(ego)[, c("Description", "geneID", "Count", "p.adjust")]

# reduce GO redundancy in CC result
ego_cc_simple <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

p_ora_cl1_cc <- dotplot(
  ego_cc_simple,
  showCategory = 10,
  title = "Cluster 1 macrophages: ORA GO Cellular Component"
)

p_ora_cl1_cc

ggsave(
  filename = "results/figures/cluster1_ORA_CC.png",
  plot = p_ora_cl1_cc,
  width = 9,
  height = 7,
  dpi = 300
)

write.csv(
  as.data.frame(ego_cc_simple),
  "results/tables/cluster1_RM_D02_vs_Naive_ORA_CC.csv",
  row.names = FALSE
)

# attempt 5: using whole sig gene list 

library(clusterProfiler)
library(org.Mm.eg.db)

sig_all <- res_shrunk %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

gene_map <- bitr(
  sig_all$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

all_map <- bitr(
  res_shrunk$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

ego <- enrichGO(
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

dotplot(ego, showCategory = 20) # none again

# trying gsea
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

rank_df <- res_shrunk %>%
  filter(!is.na(log2FoldChange)) %>%
  dplyr::select(gene, log2FoldChange)

gene_map <- bitr(
  unique(rank_df$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

rank_mapped <- rank_df %>%
  inner_join(gene_map, by = c("gene" = "SYMBOL")) %>%
  group_by(ENTREZID) %>%
  summarise(log2FoldChange = mean(log2FoldChange), .groups = "drop")

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

head(as.data.frame(gsea_bp), 20)
dotplot(gsea_bp, showCategory = 20)
# inspect what is driving each GSEA term
as.data.frame(gsea_bp)[, c("Description", "NES", "p.adjust", "setSize", "core_enrichment")]

