library(Seurat)
library(tidyverse)

obj <- readRDS("data/seurat_ass4.rds")
obj
colnames(obj[[]])

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

# using "^mt-" 
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")

VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
)

FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")

# next steps! plots do not indicate high-mito content so using a simple threshold 

obj <- subset(obj, subset = mouse_id != "")
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 20)

# 

obj <- SCTransform(obj, vars.to.regress = "percent.mt", verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.5)
obj <- RunUMAP(obj, dims = 1:30)

#
library(future)

plan(sequential)
obj <- SCTransform(
  obj,
  vars.to.regress = "percent.mt",
  method = "glmGamPoi",
  verbose = FALSE
)
