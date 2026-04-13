#cluster 4 

library(Seurat)
library(dplyr)
library(DESeq2)
library(tibble)

obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")
Idents(obj) <- "seurat_clusters"

# subset cluster 4 endothelial cells, RM only, Naive vs D02
cl4_rm <- subset(
  obj,
  subset = seurat_clusters == "4" &
    organ_custom == "RM" &
    time %in% c("Naive", "D02")
)

# check replicate sizes
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

# match metadata order to pseudobulk count columns
col_data <- col_data[match(colnames(pb_counts), col_data$pb_id), , drop = FALSE]

# make sure nothing failed to match
stopifnot(all(colnames(pb_counts) == col_data$pb_id))

rownames(col_data) <- col_data$pb_id
col_data$time <- factor(col_data$time, levels = c("Naive", "D02"))

# DESeq2
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

# for clean summary numbers
res_sig <- res %>%
  dplyr::filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)

res_sig_interpretable <- res_sig %>%
  dplyr::filter(!grepl("^Rpl|^Rps|^mt-", gene))

nrow(res %>% dplyr::filter(!is.na(padj), padj < 0.05))
nrow(res_sig)
nrow(res_sig_interpretable)
sum(res_sig$log2FoldChange > 0)
sum(res_sig$log2FoldChange < 0)

# applying shrinkage
library(apeglm)

res_lfc <- lfcShrink(dds, coef = "time_D02_vs_Naive", type = "apeglm")
plotMA(res_lfc, ylim = c(-4, 4))