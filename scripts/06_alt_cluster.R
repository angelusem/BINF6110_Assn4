obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")
library(dplyr)
library(tidyr)
library(readr)

# sanity check
stopifnot("seurat_clusters" %in% colnames(obj@meta.data))
stopifnot("organ_custom" %in% colnames(obj@meta.data))
stopifnot("time" %in% colnames(obj@meta.data))
stopifnot("mouse_id" %in% colnames(obj@meta.data))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

meta <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    seurat_clusters = as.character(seurat_clusters),
    organ_custom = as.character(organ_custom),
    time = as.character(time),
    mouse_id = as.character(mouse_id)
  )

# 1) cluster x tissue x time x mouse counts
cluster_counts_by_mouse <- meta %>%
  dplyr::count(seurat_clusters, organ_custom, time, mouse_id, name = "n_cells") %>%
  dplyr::arrange(as.numeric(seurat_clusters), organ_custom, time, mouse_id)

write_csv(
  cluster_counts_by_mouse,
  "results/tables/cluster_counts_by_mouse.csv"
)

# 2) cluster x tissue x time totals

cluster_counts_by_tissue_time <- meta %>%
  dplyr::count(seurat_clusters, organ_custom, time, name = "n_cells") %>%
  dplyr::arrange(as.numeric(seurat_clusters), organ_custom, time)


write_csv(
  cluster_counts_by_tissue_time,
  "results/tables/cluster_counts_by_tissue_time.csv"
)

# 3) overall cluster sizes
cluster_sizes <- meta %>%
  dplyr::count(seurat_clusters, name = "total_cells") %>%
  dplyr::arrange(as.numeric(seurat_clusters))


write_csv(
  cluster_sizes,
  "results/tables/cluster_sizes.csv"
)

# 4) quick screen for pairwise comparisons with 3 mice per group
# adjust candidate clusters here if you want
candidate_clusters <- c("2", "4", "5", "6", "13", "15", "31")

pairwise_screen <- cluster_counts_by_mouse %>%
  filter(seurat_clusters %in% candidate_clusters) %>%
  group_by(seurat_clusters, organ_custom, time) %>%
  summarise(
    mice_present = n_distinct(mouse_id),
    min_cells_per_mouse = min(n_cells),
    median_cells_per_mouse = median(n_cells),
    total_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  arrange(as.numeric(seurat_clusters), organ_custom, time)

write_csv(
  pairwise_screen,
  "results/tables/candidate_cluster_screen.csv"
)

# 5) optional wide table for easy visual inspection
pairwise_screen_wide <- pairwise_screen %>%
  unite("group", organ_custom, time, sep = "_") %>%
  select(seurat_clusters, group, mice_present, min_cells_per_mouse, total_cells) %>%
  pivot_wider(
    names_from = group,
    values_from = c(mice_present, min_cells_per_mouse, total_cells)
  ) %>%
  arrange(as.numeric(seurat_clusters))

write_csv(
  pairwise_screen_wide,
  "results/tables/candidate_cluster_screen_wide.csv"
)

message("Done. Files written to results/tables/")

# then, running this for my candidate genes that popped out in my top markers by cluster table
library(dplyr)
library(readr)

markers <- read_csv("top_markers_by_cluster.csv", show_col_types = FALSE)

candidate_clusters <- c("2", "4", "5", "6", "13", "15", "31")

candidate_markers <- markers %>%
  mutate(cluster = as.character(cluster)) %>%
  filter(cluster %in% candidate_clusters)

write_csv(
  candidate_markers,
  "results/tables/candidate_cluster_top_markers.csv"
)

candidate_markers