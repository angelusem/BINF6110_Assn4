# 12_additional_analysis: to see whether macrophages and endothelial cells appear to shift in relative abundance between Naive and D02
## for each mouse, calculating the fraction of RM cells that are cluster 1 or cluster 4, then making a plot split by time


library(dplyr)
library(ggplot2)

obj <- readRDS("data/seurat_ass4_lognorm_clustered.rds")

comp_df <- obj[[]] %>%
  filter(organ_custom == "RM", time %in% c("Naive", "D02")) %>%
  mutate(
    cluster_group = case_when(
      seurat_clusters == "1" ~ "Cluster 1 macrophage",
      seurat_clusters == "4" ~ "Cluster 4 endothelial",
      TRUE ~ "Other"
    )
  ) %>%
  count(mouse_id, time, cluster_group, name = "n_cells") %>%
  group_by(mouse_id, time) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  filter(cluster_group != "Other")

# 1 : per-mouse

p_comp_points <- ggplot(comp_df, aes(x = time, y = prop, color = cluster_group)) +
  geom_point(position = position_jitter(width = 0.05, height = 0), size = 3) +
  facet_wrap(~ cluster_group) +
  theme_bw() +
  labs(
    y = "Proportion of RM cells",
    x = NULL
  )

ggsave(
  filename = "results/figures/cluster_composition_RM_clusters1_4_points_only.png",
  plot = p_comp_points,
  width = 8,
  height = 5.5,
  dpi = 300
)


# 2: with group mean marker

p_comp_mean <- ggplot(comp_df, aes(x = time, y = prop, color = cluster_group)) +
  geom_point(position = position_jitter(width = 0.05, height = 0), size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  facet_wrap(~ cluster_group) +
  theme_bw() +
  labs(
    y = "Proportion of RM cells",
    x = NULL
  )

p_comp_mean

ggsave(
  filename = "results/figures/cluster_composition_RM_clusters1_4.png",
  plot = p_comp_mean,
  width = 8,
  height = 5.5,
  dpi = 300
)