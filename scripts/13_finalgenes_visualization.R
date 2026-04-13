library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# load clustered object
obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")
Idents(obj) <- "seurat_clusters"

# final discussion genes
genes_main <- c("Gpr183", "Hk2", "Ap2b1", "Hes1", "Cd200")
genes_alt  <- c("Ap2a2")  # optional swap-in if you want Ap2a2 instead of Ap2b1

# choose one final set
genes_to_plot <- genes_main

# check which genes exist in the object
genes_present <- genes_to_plot[genes_to_plot %in% rownames(obj)]
genes_missing <- setdiff(genes_to_plot, genes_present)

if (length(genes_missing) > 0) {
  message("These genes were not found in the object and will be skipped: ",
          paste(genes_missing, collapse = ", "))
}

stopifnot(length(genes_present) > 0)


# 1. FeaturePlot panel on full UMAP

p_feature <- FeaturePlot(
  obj,
  features = genes_present,
  reduction = "umap",
  ncol = 2,
  order = TRUE
)

ggsave(
  filename = "results/figures/final_gene_panel_featureplot.png",
  plot = p_feature,
  width = 12,
  height = 10,
  dpi = 300
)

# 2. DotPlot across all clusters

p_dot <- DotPlot(
  obj,
  features = genes_present,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  labs(
    title = "Selected discussion genes across clusters",
    x = "Gene",
    y = "Cluster"
  ) +
  theme_bw()

ggsave(
  filename = "results/figures/final_gene_panel_dotplot.png",
  plot = p_dot,
  width = 8,
  height = 5,
  dpi = 300
)

# 3. VlnPlot focused only on clusters 1 and 4

obj_14 <- subset(obj, idents = c("1", "4"))

# nicer labels
obj_14$cluster_label <- dplyr::case_when(
  obj_14$seurat_clusters == "1" ~ "1_macrophage",
  obj_14$seurat_clusters == "4" ~ "4_endothelial",
  TRUE ~ as.character(obj_14$seurat_clusters)
)

p_vln <- VlnPlot(
  obj_14,
  features = genes_present,
  group.by = "cluster_label",
  pt.size = 0,
  ncol = 2
)

ggsave(
  filename = "results/figures/final_gene_panel_vlnplot_clusters1_4.png",
  plot = p_vln,
  width = 10,
  height = 8,
  dpi = 300
)


# 4. one combined multi-panel figure

# For a cleaner combined figure, using a reduced feature set
genes_combined <- genes_present[1:min(length(genes_present), 4)]

p_feature_small <- FeaturePlot(
  obj,
  features = genes_combined,
  reduction = "umap",
  ncol = 2,
  order = TRUE
)

p_dot_small <- DotPlot(
  obj,
  features = genes_combined,
  group.by = "seurat_clusters"
) +
  RotatedAxis() +
  theme_bw() +
  labs(title = "Dot plot")

p_combined <- p_feature_small / p_dot_small +
  plot_annotation(title = "Final discussion gene panel")

ggsave(
  filename = "results/figures/final_gene_panel_combined.png",
  plot = p_combined,
  width = 12,
  height = 12,
  dpi = 300
)

message("Saved figures to results/figures/")