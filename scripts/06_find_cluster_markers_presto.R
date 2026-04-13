library(Seurat)
library(dplyr)
library(presto)

obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")
Idents(obj) <- "seurat_clusters"

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10) %>%
  ungroup()

write.csv(markers, "results/tables/all_cluster_markers.csv", row.names = FALSE)
write.csv(top_markers, "results/tables/top_markers_by_cluster.csv", row.names = FALSE)

saveRDS(markers, "results/tables/all_cluster_markers.rds")
saveRDS(top_markers, "results/tables/top_markers_by_cluster.rds")

top_markers
#
#print all and look at them 
print(top_markers, n = 360)
# making a compact annotation table

annotation_table <- top_markers %>%
  group_by(cluster) %>%
  summarise(
    top10 = paste(gene, collapse = ", "),
    .groups = "drop"
  )

annotation_table

write.csv(annotation_table, "results/tables/annotation_helper_table.csv", row.names = FALSE)

# after curation of manual annotation, using feature plots to confirm a smaller subset
library(Seurat)
library(dplyr)

obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")
Idents(obj) <- "seurat_clusters"

FeaturePlot(obj, features = c("Cnga4", "Ano2", "Adcy3"))
FeaturePlot(obj, features = c("Flt1", "Pecam1", "Cdh5"))
FeaturePlot(obj, features = c("Ctss", "Lyz2", "C1qa"))
FeaturePlot(obj, features = c("Col1a1", "Col1a2", "Prrx1"))
FeaturePlot(obj, features = c("Krt8", "Krt18", "Krt19"))

# To connect that region to a numbered cluster, put the labeled cluster UMAP next to the feature plot and compare positions-- why would i do next- where did i get the above groups from??

DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE)

# dot plot for clearer view
DotPlot(
  obj,
  features = c(
    "Cnga4","Ano2","Adcy3",
    "Flt1","Pecam1","Cdh5",
    "Col1a1","Col1a2","Prrx1",
    "Ctss","Lyz2","C1qa",
    "Krt8","Krt18","Krt19"
  )
) + RotatedAxis()