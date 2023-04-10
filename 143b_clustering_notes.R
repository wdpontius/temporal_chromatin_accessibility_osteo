# setwd("~/143b_atac/FASTQ/Filtered_bams_round1")

library(preprocessCore)
library(wesanderson)
library(cowplot)
library(tidyverse)

options(scipen = 999)

files <- list.files(pattern = glob2rx("143b*ATAC*day1_combined*txt"))

names(files) <- files

df <-
  map_dfr(files, ~ read_tsv(.x, col_names = FALSE), .id = "file")

# df_sample <-
#   df %>%
#   mutate(sample = case_when(
#     str_detect(file, "in_vitro") ~ "invitro",
#     str_detect(file, "Day_1") ~ "early",
#     str_detect(file, "Day_22") ~ "late"
#   ))

matrix_rpkm <-
  df %>%
  unite("peak", c("X1", "X2", "X3")) %>%
  pivot_wider(names_from = file, values_from = X4) %>%
  column_to_rownames("peak") %>%
  as.matrix()


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

p <- 
  df_pc_anno %>%
  mutate(type = case_when(
    str_detect(rowname, "in_vitro") ~ "invitro",
    str_detect(rowname, "Day_1") ~ "early",
    str_detect(rowname, "Day_22") ~ "late"
  )) %>%
  mutate(type = fct_relevel(type, "invitro", "early", "late")) %>%
  ggplot() + 
  geom_point(aes(x = PC1, y = PC2, color = type), size = 2) +
  theme_cowplot() +
  background_grid() + 
  scale_color_manual(values = wes_palette("BottleRocket2")) +
  labs(x = "PC1 - 68.6% variance explained",
       y = "PC2 - 20.2% variance explained") +
  panel_border(color = "black") +
  theme(
    aspect.ratio = 1
  )

save_plot("143b_pca.pdf", p)


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
  mutate(type = case_when(
    str_detect(sample, "in_vitro") ~ "invitro",
    str_detect(sample, "Day_1") ~ "early",
    str_detect(sample, "Day_22") ~ "late"
  )) %>%
  # mutate(type = str_replace(sample, "_.*", "")) %>%
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

# cluster into 8 clusters based on shape
clustering <- kmeans(matrix_dynamic_zscore, 7)

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

# visualize clusters

p <- 
  df_combined %>%
  pivot_longer(c(early, invitro, late),
               names_to = "condition",
               values_to = "z_score") %>%
  mutate(condition = fct_relevel(condition,
                                 "invitro",
                                 "early",
                                 "late")) %>%
  ggplot(aes(condition, z_score)) +
  geom_line(aes(group = peak)) +
  facet_grid(.~cluster)

df_combined %>%
  pivot_longer(c(early, invitro, late),
               names_to = "condition",
               values_to = "z_score") %>%
  mutate(condition = fct_relevel(condition,
                                 "invitro",
                                 "early",
                                 "late")) %>%
  ggplot(aes(condition, z_score)) +
  stat_summary(fun.y = "mean", geom = "line", aes(group = 1)) +
  facet_grid(.~cluster)

# ggsave("143b_cluster.png", p, type = "cairo-png")

reorder_clusters <-
  tibble(cluster = c(2,4,3,6,5,1,7,"stable"),
         recluster = c(1:7,"stable"))

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
  # geom_line(aes(group = peak), alpha = 0.5) +
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

ggsave("143b_cluster_reorder.png", p, type = "cairo-png")
ggsave("143b_cluster_reorder_wide.png", p, type = "cairo-png", width = 20, height = 10)

df_plot %>%
  distinct(peak, recluster) %>%
  separate(peak, c("chr", "start", "stop"), sep = "_") %>%
  mutate_at(vars(start, stop), ~as.numeric(.)) %>%
  nest(-recluster) %>%
  pwalk( ~ write_tsv(..2, paste0("143b_recluster_", ..1, ".bed"), col_names = FALSE))


# make data for gimme motifs

df_gimme <- 
  df_plot %>%
  select("loc" = peak,
         "cluster" = recluster) %>%
  mutate(loc = str_replace(loc, "_", ":"),
         loc = str_replace(loc, "_", "-")) 

write_tsv(df_gimme, "143b_gimmemotifs_input.txt", col_names = TRUE)
  

# get only largest change peaks -------------------------------------------

peaks_bigchange <-
  df_rpkm_qnorm_mean %>%
  pivot_wider(names_from = type,
              values_from = mean) %>%
  mutate(early_fc = log2(early / invitro),
         late_fc = log2(late / invitro)) %>%
  filter(abs(early_fc) > 1 | abs(late_fc) > 1) %>%
  pull(peak)

df_plot %>%
  filter(peak %in% peaks_bigchange) %>%
  distinct(peak, recluster) %>%
  separate(peak, c("chr", "start", "stop"), sep = "_") %>%
  mutate_at(vars(start, stop), ~as.numeric(.)) %>%
  nest(-recluster) %>%
  pwalk( ~ write_tsv(..2, paste0("143b_recluster_log2fc-gt-1_", ..1, ".bed"), col_names = FALSE))



df_gimme <- 
  df_plot %>% 
  select("loc" = peak,
         "cluster" = recluster) %>%
  # this step is to get rid of parts of bed file with scientific notation (2e+05)
  mutate(loc = str_replace(loc, ":", "_")) %>%
  separate("loc", c("chr", "start", "stop"), sep = "_") %>% 
  mutate_at(vars(start, stop), ~as.numeric(.)) %>% 
  unite("loc", c("chr", "start", "stop"), sep = "_") %>%
  # reformat for gimme motifs input file
  mutate(loc = str_replace(loc, "_", ":"),
         loc = str_replace(loc, "_", "-"),
         cluster = paste0("cluster", cluster)) %>%
  distinct(loc, cluster)

write_tsv(df_gimme, "143b_gimmemotifs_input.txt", col_names = TRUE)
