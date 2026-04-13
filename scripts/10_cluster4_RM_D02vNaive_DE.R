# cluster 4 
library(Seurat)
library(dplyr)
library(DESeq2)

obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")

cl4_rm <- subset(
  obj,
  subset = seurat_clusters == "4" &
    organ_custom == "RM" &
    time %in% c("Naive", "D02")
)

cl4_rm$pb_id <- cl4_rm$mouse_id

pb <- AggregateExpression(
  cl4_rm,
  assays = "RNA",
  group.by = "pb_id",
  return.seurat = TRUE
)

pb_counts <- GetAssayData(pb, assay = "RNA", layer = "counts")

col_data <- cl4_rm@meta.data %>%
  dplyr::select(mouse_id, time) %>%
  distinct() %>%
  mutate(pb_id = gsub("_", "-", mouse_id)) %>%
  tibble::column_to_rownames("pb_id")

col_data <- col_data[colnames(pb_counts), , drop = FALSE]
col_data$time <- factor(col_data$time, levels = c("Naive", "D02"))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pb_counts)),
  colData = col_data,
  design = ~ time
)

dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("time", "D02", "Naive"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

write.csv(res_df, "results/tables/cluster4_RM_D02_vs_Naive_DESeq2.csv", row.names = FALSE)