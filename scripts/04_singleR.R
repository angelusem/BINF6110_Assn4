#07_singler_annotation.R
#BiocManager::install("celldex")

library(Seurat)
library(SingleR)
library(celldex)
library(dplyr)

# load your clustered object
obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")
Idents(obj) <- "seurat_clusters"

# use the log-normalized RNA data for annotation
expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")

# mouse reference
ref <- MouseRNAseqData()

# annotate clusters, not individual cells
cluster_pred <- SingleR(
  test = expr_mat,
  ref = ref,
  labels = ref$label.main,
  clusters = obj$seurat_clusters
)

# make a clean results table
cluster_labels <- data.frame(
  cluster = rownames(cluster_pred),
  singler_label = cluster_pred$labels,
  singler_pruned = cluster_pred$pruned.labels,
  stringsAsFactors = FALSE
)

cluster_labels <- cluster_labels %>%
  arrange(as.numeric(as.character(cluster)))

print(cluster_labels)

write.csv(
  cluster_labels,
  "results/tables/singler_cluster_labels.csv",
  row.names = FALSE
)

# map SingleR labels back onto the Seurat object
obj$SingleR_label <- cluster_labels$singler_label[
  match(as.character(obj$seurat_clusters), cluster_labels$cluster)
]

obj$SingleR_pruned <- cluster_labels$singler_pruned[
  match(as.character(obj$seurat_clusters), cluster_labels$cluster)
]

# save the updated object
saveRDS(obj, "data/seurat_ass4_lognorm_clustered_singler.rds")

# UMAP colored by SingleR labels
p1 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "SingleR_label",
  label = TRUE,
  repel = TRUE
)

p2 <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "SingleR_pruned",
  label = TRUE,
  repel = TRUE
)

print(p1)
print(p2)

ggsave(
  filename = "results/figures/umap_singler_label.png",
  plot = p1,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = "results/figures/umap_singler_pruned.png",
  plot = p2,
  width = 10,
  height = 8,
  dpi = 300
)

# optional: score heatmap
pdf("results/figures/singler_score_heatmap.pdf", width = 10, height = 8)
plotScoreHeatmap(cluster_pred)
dev.off()