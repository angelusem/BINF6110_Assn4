library(Seurat)

obj <- readRDS("data/seurat_ass4.rds")

# QC metrics
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

# basic filtering
obj <- subset(obj, subset = mouse_id != "")
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)

# lighter normalization path 
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