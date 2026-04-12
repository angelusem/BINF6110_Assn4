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
  count(mouse_id, time, name = "n_cells") %>%
  arrange(time, mouse_id)

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

# reorder count matrix to match metadata
pb_counts <- pb_counts[, rownames(col_data)]

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