#01_load_data

library(Seurat)

obj <- readRDS("data/seurat_ass4.rds")
obj
colnames(obj[[]])

#scripts/02_metadata_and_qc_setup.R =================================================================
library(Seurat)
library(tidyverse)

obj <- readRDS("data/seurat_ass4.rds")

# 1. Inspect the metadata values we actually have
lapply(
  obj[[]][, c("biosample_id", "organ_custom", "time", "mouse_id", "disease__ontology_label")],
  function(x) sort(unique(as.character(x)))
)

# 2. Cell counts per metadata field
table(obj$organ_custom)
table(obj$time)
table(obj$disease__ontology_label)

# 3. Cross-tab tissue by time
table(obj$organ_custom, obj$time)

# 4. Sample-level summary
sample_summary <- obj[[]] %>%
  as_tibble(rownames = "cell_id") %>%
  count(biosample_id, mouse_id, organ_custom, time, disease__ontology_label, name = "n_cells") %>%
  arrange(organ_custom, time, biosample_id)

sample_summary
write.csv(sample_summary, "results/tables/sample_summary.csv", row.names = FALSE)

# 5. Check mitochondrial gene naming before QC
sum(grepl("^mt-", rownames(obj)))
sum(grepl("^MT-", rownames(obj)))
head(rownames(obj), 20)

# QC
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
)

FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
# if the above is too complex, i might have run this:
library(Seurat)
library(dplyr)

obj <- readRDS("data/seurat_ass4.rds")

unique(obj$organ_custom)
unique(obj$time)

table(obj$organ_custom, obj$time)

sample_summary <- obj[[]] %>%
  count(biosample_id, mouse_id, organ_custom, time, disease__ontology_label, name = "n_cells") %>%
  arrange(organ_custom, time)

sample_summary

write.csv(sample_summary, "results/tables/sample_summary.csv", row.names = FALSE)

# scripts/03_check_mouse_id.R: NOT SURE IF I DID THIS !!! ===================================
library(Seurat)

obj <- readRDS("data/seurat_ass4.rds")

table(obj$mouse_id == "")
table(obj$organ_custom, obj$mouse_id == "")
table(obj$time, obj$mouse_id == "")

# If the blank mouse_id count is small, we will just remove those cells before downstream group comparisons.

# scripts/04_qc_and_filter.R =================================================================

library(Seurat)

obj <- readRDS("data/seurat_ass4.rds")

# mitochondrial percentage
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# remove cells with missing mouse_id
obj <- subset(obj, subset = mouse_id != "")

# basic QC plots
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# first-pass filter
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)

# save filtered object
saveRDS(obj, "data/seurat_ass4_qc.rds")

# scripts/05_cluster_umap.R: =============================================================================
library(Seurat)

obj <- readRDS("data/seurat_ass4_qc.rds")

obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

ElbowPlot(obj)

obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)

DimPlot(obj, reduction = "umap", label = TRUE)
DimPlot(obj, reduction = "umap", group.by = "organ_custom")
DimPlot(obj, reduction = "umap", group.by = "time")

saveRDS(obj, "data/seurat_ass4_clustered.rds") 

# did this but different!!! --- used log-normalized instead?=============================== still 05

library(Seurat)

obj <- readRDS("data/seurat_ass4.rds")

# QC metrics
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# basic filtering
obj <- subset(obj, subset = mouse_id != "")
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)

# lighter normalization path from lecture workflow
obj <- NormalizeData(obj, normalization.method = "LogNormalize")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  vars.to.regress = "percent.mt"
)
obj <- RunPCA(obj, features = VariableFeatures(obj))

ElbowPlot(obj)

obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)

DimPlot(obj, reduction = "umap", label = TRUE)
DimPlot(obj, reduction = "umap", group.by = "organ_custom")
DimPlot(obj, reduction = "umap", group.by = "time")

saveRDS(obj, "data/seurat_ass4_lognorm_clustered.rds")

# scripts/06_find_cluster_markers_presto.R ========================================================

# 07_singler_annotation.R ============================ only for the purpose of validation-- why? is this not reliable or what 