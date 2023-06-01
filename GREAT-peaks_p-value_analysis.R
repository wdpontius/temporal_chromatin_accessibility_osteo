# ----------------- anova for p-values
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

matrix_rpkm_qnorm <- preprocessCore::normalize.quantiles(matrix_rpkm)

# replace lost row and column names after quantile normalization
colnames(matrix_rpkm_qnorm) <- colnames(matrix_rpkm)
rownames(matrix_rpkm_qnorm) <- rownames(matrix_rpkm)

df_rpkm_qnorm_long <- 
  matrix_rpkm_qnorm %>%
  as.data.frame() %>%
  rownames_to_column(var = "peak") %>%
  pivot_longer(-peak, names_to = "sample", values_to = "value") %>%
  mutate(type = str_replace(sample, "_.*", ""))


model_fit <- function(data){aov(value ~ type, data = data)}

df_rpkm_qnorm_long_p <- 
  df_rpkm_qnorm_long %>%
  select(-sample) %>%
  group_nest(peak) %>%
  mutate(model = map(data, model_fit)) %>%
  mutate(results = map(model, summary)) %>%
  mutate(p = map(results, ~.x[[1]]$Pr[1])) %>%
  unnest(p)

# read in peaks used for GREAT previously
files <- list.files(pattern = glob2rx("2020-04-24_cluster-redo_recluster*log2fc-gt-1.bed"))

names(files) <- paste0("cluster_", 1:8)

df_log2fc <- map_dfr(files, ~read_tsv(.x, col_names = FALSE), .id = "cluster")

peaks <- 
  df_log2fc %>%
  unite("peak", c("X1", "X2", "X3")) %>%
  pull(peak)

# check how many of great peaks are not significantly different due to anova
df_rpkm_qnorm_long_p %>%
      filter(peak %in% peaks) %>%
      count(p < 0.05)


# 1494 p >= 0.05, 42074 < 0.05. ~3% are not significant.
