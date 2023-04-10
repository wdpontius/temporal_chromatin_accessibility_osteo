library(tidyverse)
library(cowplot)
library(broom)
library(valr)

numbers <- 1:7

list_633 <- map(numbers, ~ read_bed(
  paste0(
    "/Users/deanbot/Desktop/projects/tracks/mg633/atac/2020-04-24_cluster-redo_recluster-",
    .x,
    ".bed"
  )
))

list_143b <- map(numbers, ~ read_bed(
  paste0(
    "/Users/deanbot/Desktop/projects/crispr-screen/2022-03-01_143b-atac/143b_recluster_",
    .x,
    ".bed"
  )
))

x633_x143b_universe <- read_bed("/Users/deanbot/Desktop/projects/crispr-screen/2022-03-01_143b-atac/143b_633_combined_universe_autosome_noblacklist.bed")

results_list <- 
  expand_grid("x" = 1:7, "y" = 1:7) %>%
  pmap(function(...){x633_in_x143b <- nrow(bed_intersect(list_633[[..1]], list_143b[[..2]]))
  
  x633_notin_x143b <-
    nrow(bed_subtract(list_633[[..1]], list_143b[[..2]], any = TRUE))
  
  x143b_notin_x633 <-
    nrow(bed_subtract(list_143b[[..2]], list_633[[..1]], any = TRUE))
  
  # notin_other <- nrow(bed_subtract(
  #   bind_rows(list_633,
  #             list_143b),
  #   bind_rows(list_633[[..1]],
  #             list_143b[[..2]])
  # ))
  
  notin_other <-
  nrow(bed_subtract(x633_x143b_universe, bind_rows(list_633[[..1]], list_143b[[..2]])))
  
  fisher_df <-
    data.frame(
      "in_b" = c(x633_in_x143b, x143b_notin_x633),
      "not_in_b" = c(x633_notin_x143b, notin_other)
    )
  
  rownames(fisher_df) <- c("in_a", "not_in_a")
  
  tidy(fisher.test(fisher_df)) %>% mutate(x633 = ..1, x143b = ..2)
  })

df_results <- bind_rows(results_list)

# plot heatmap showing all vs. all comparison
p_results <- 
  df_results %>%
  mutate(shape = case_when(x633 == x143b ~ "identity",
                           TRUE ~ "other")) %>%
  rename("odds_ratio" = "estimate") %>%
  ggplot(aes(factor(x633), factor(x143b))) +
  geom_tile(aes(fill = -log10(p.value))) +
  geom_point(aes(size = odds_ratio, color = odds_ratio)) +
  # geom_point(aes(shape = shape), fill = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "red") +
  scale_color_gradient2(low = "yellow", mid = "white", high = "black", midpoint = 1) +
  # scale_shape_manual(values = c(21, 26)) + 
  labs(x = "mg633 clusters",
       y = "143b clusters") +
  theme_cowplot() +
  panel_border(color = "black")

ggsave("cluster_overlap_stats.pdf", p_results)

# plot barplot just comparing equivalent clusters
p_cluster_overlap_bar <- 
  df_results %>%
  filter(x633 == x143b) %>%
  ggplot(aes(x = factor(x633), y = estimate)) +
  geom_col(fill = "red") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.25) +
  geom_text(label = "***", aes(y = conf.high + 0.1)) +
  geom_hline(aes(yintercept = 1), lty = "dashed") +
  geom_text(aes(x = "1", y = 1.5), label = "OR = 1") +
  labs(x = "Cluster", y = "Odds ratio")

ggsave("cluster_overlap_barplot.pdf", height = 4, width = 6)

