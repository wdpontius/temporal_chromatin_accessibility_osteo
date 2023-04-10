# setwd("/mnt/rds/genetics01/ScacheriLab/wdp14/ATAC/in_vivo/Filtered_bams")

library(glue)
library(wesanderson)    # color palette for figures
library(stringr)
library(patchwork)
library(cowplot)        # ggplot formatting
library(tidyverse)

files <- list.files(pattern = "rm-dup_rpkm.txt")

names(files) <- c(
  paste(rep("invitro", 3), 1:3, sep = "_"),
  paste(rep("early", 2),   1:2, sep = "_"),
  paste(rep("late", 5),    1:5, sep = "_")
)

df_rpkm_long <- map_dfr(files, ~read_tsv(.x, col_names = FALSE), .id = "sample")

df_rpkm_wide <-
  df_rpkm_long %>%
  pivot_wider(names_from = sample, values_from = X4) %>%
  unite("peak", c("X1", "X2", "X3"))

matrix_rpkm <- 
  df_rpkm_wide %>%
  column_to_rownames("peak") %>%
  as.matrix()

# principal component analysis -----------------------------------
# log2 transform
matrix_log2_rpkm <- apply(matrix_rpkm, 2, function(x) log2(x))

# qnorm
matrix_log2_rpkm_qnorm <- preprocessCore::normalize.quantiles(matrix_log2_rpkm)

# lose row and column names after quantile normalization
colnames(matrix_log2_rpkm_qnorm) <- colnames(matrix_log2_rpkm)
rownames(matrix_log2_rpkm_qnorm) <- rownames(matrix_log2_rpkm)

pc_log2_rpkm <- prcomp(matrix_log2_rpkm_qnorm)

df_pc_anno <-
  pc_log2_rpkm$rotation %>%
  # use data frame instead of tibble to retain rownames
  as.data.frame() %>% 
  rownames_to_column() %>%
  # pivot_longer(-rowname, names_to = "principal_component", values_to = "value") %>%
  mutate(type = str_replace(rowname, "_.*", ""))


# plot samples to see clustering
# ggplot(df_pc_anno) +
#   geom_point(aes(x = PC1, y = PC2, color = type)) 

# k-means clustering ---------------------------------------------
# prepare matrix
matrix_rpkm_qnorm <- preprocessCore::normalize.quantiles(matrix_rpkm)

# replace lost row and column names after quantile normalization
colnames(matrix_rpkm_qnorm) <- colnames(matrix_rpkm)
rownames(matrix_rpkm_qnorm) <- rownames(matrix_rpkm)

# make data frame of mean for each group
df_rpkm_qnorm_mean <-
  matrix_rpkm_qnorm %>%
  as.data.frame() %>%
  rownames_to_column(var = "peak") %>%
  pivot_longer(-peak, names_to = "sample", values_to = "value") %>%
  mutate(type = str_replace(sample, "_.*", "")) %>%
  group_by(type, peak) %>%
  summarise(mean = mean(value)) %>%
  ungroup()

# Step 1: filter stable/dynamic peaks---------------------
df_rpkm_qnorm_mean_var <-
  df_rpkm_qnorm_mean %>%
  pivot_wider(names_from = type,
              values_from = mean) %>%
  rowwise() %>%
  mutate(coeff_variation = sd(c(early, invitro, late))/mean(c(early, invitro, late))*100) %>%
  ungroup()

# filter dynamic and static peaks based on 
# variability between conditions
static_peaks <- df_rpkm_qnorm_mean_var %>% filter(coeff_variation < 10)
dynamic_peaks <- df_rpkm_qnorm_mean_var %>% filter(coeff_variation >= 10)

# Step 2: cluster dynamic peaks with k-means --------------
# z-score dynamic peaks

matrix_dynamic <-
  dynamic_peaks %>%
  select(-coeff_variation) %>%
  column_to_rownames("peak") %>%
  as.matrix()

# have to transpose in order to calculate
# z-score across rows
matrix_dynamic_zscore <- t(scale(t(matrix_dynamic)))

set.seed(216)

clustering <- kmeans(matrix_dynamic_zscore, 8)

df_cluster <-
  as.data.frame(clustering$cluster) %>%
  dplyr::rename("cluster" = "clustering$cluster") %>%
  rownames_to_column(var = "peak")

df_dynamic_zscore <-
  as.data.frame(matrix_dynamic_zscore) %>%
  rownames_to_column(var = "peak")

df_dynamic_zscore_clusters <-
  left_join(df_dynamic_zscore, df_cluster, by = "peak") %>%
  bind_rows()

# z score static peaks and combine with dynamic
# for visualization
static_peaks_zscore <-
  static_peaks %>%
  select(-coeff_variation) %>%
  mutate(cluster = "stable") %>%
  pivot_longer(-c(peak,cluster),
               names_to = "condition",
               values_to = "value") %>%
  group_by(peak) %>%
  mutate(value = (value - mean(value))/sd(value)) %>%
  ungroup() %>%
  pivot_wider(names_from = condition,
              values_from = value)

df_combined <-
  df_dynamic_zscore_clusters %>%
  mutate(cluster = as.character(cluster)) %>%
  bind_rows(static_peaks_zscore)

# plot dynamics of clusters
# reorder clusters for aesthetic purposes
reorder_clusters <-
  tibble(cluster = c(8,3,4,5,2,6,1,7, "stable"),
         recluster = c(1:8,"stable"))

df_plot <-
  df_combined %>%
  add_count(cluster) %>% 
  left_join(reorder_clusters, by = "cluster") %>%
  pivot_longer(c(early, invitro, late),
               names_to = "condition",
               values_to = "value") %>%
  mutate(condition = fct_relevel(condition,
                                 "invitro",
                                 "early",
                                 "late"))

p <-
  df_plot %>% 
  filter(cluster != "stable") %>% 
  ggplot(aes(condition, value)) +
  geom_line(aes(group = peak), alpha = 0.5) +
  stat_summary(fun.y = "median", geom = "line",
               aes(group = 1),
               col = "coral") +
  theme_cowplot() + 
  panel_border(color = "black") +
  background_grid(major = "y") +
  facet_grid(glue::glue("{recluster}: {n} peaks")~.) +
  labs(x = NULL, y = "z-score") +
  theme(
    # aspect.ratio = 1,
    axis.line = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.background = element_blank(),
    panel.spacing = unit(0.1, "lines")
  )

# ggsave("2020-04-27_cluster-redo.png", p, type = "cairo", width = 20, height = 10)

p2 <-
  df_combined %>%
  # filter(cluster != "stable") %>% 
  add_count(cluster) %>% 
  left_join(reorder_clusters, by = "cluster") %>%
  group_by(recluster) %>% 
  sample_frac(size = 0.1) %>% 
  ungroup() %>% 
  pivot_longer(c(early, invitro, late),
               names_to = "condition",
               values_to = "z_score") %>%
  mutate(condition = fct_relevel(condition,
                                 "invitro",
                                 "early",
                                 "late")) %>% 
  ggplot(aes(condition, peak)) +
  geom_tile(aes(fill = z_score)) +
  # annotate("segment",x=Inf,xend=Inf,y=Inf,yend=-Inf,color="black") +
  facet_grid(glue::glue("Cluster {recluster}: {n} peaks")~.,
             scales = "free", space = "free") +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "firebrick") +
  labs(x = NULL) +
  theme_cowplot() +
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "left",
    legend.text.align = 1,
    panel.spacing = unit(0.1, "lines")
  )

# ggsave("2020-04-30_heatmap.png", p2, type = "cairo")
