library(tidyverse)
library(janitor)
library(preprocessCore)
library(patchwork)
library(tidytext)
library(simplifyEnrichment)
library(widyr)

# analysis based of this https://www.badgrammargoodsyntax.com/compbio/2017/12/16/compbio-017-is-your-overlap-significant

great_143b <-
  map_dfr(c(1:7),
          ~ read_tsv(
            paste0(
              "143b_recluster",
              .x,
              "_all-ontology.tsv"
            ),
            skip = 3
          ) %>%
            clean_names() %>%
            mutate(cell_line = "x143b"),
          .id = "cluster")

great_633 <- 
  map_dfr(c(1:7),
          ~ read_tsv(
            paste0(
              "633_recluster",
              .x,
              "_all-ontology.tsv"
            ),
            skip = 3
          ) %>%
            clean_names() %>%
            mutate(cell_line = "x633"),
          .id = "cluster")

df_combined <- bind_rows(great_143b, great_633)

df_combined %>%
  ggplot(aes(hyper_fdr_q)) +
  geom_density(aes(fill = cell_line), alpha = 0.5) +
  facet_grid(.~cluster)

# test on cluster 1 first



df_results <- map_dfr(1:7,
    function(x) {
      # significant terms in 633, per cluster
      sig_633_terms <-
        df_combined %>%
        filter(cluster == x) %>%
        filter(str_detect(id, "GO")) %>%
        filter(hyper_fdr_q < 0.005, cell_line == "x633") %>%
        pull(id)
      
      # significant terms in 143b, per cluster
      sig_143b_terms <- 
        df_combined %>%
        filter(cluster == x) %>% 
        filter(str_detect(id, "GO")) %>%
        filter(hyper_fdr_q < 0.005, cell_line == "x143b") %>%
        pull(id)
      
      # number of sig 143b terms
      sig_143b <- length(sig_143b_terms)
      
      # number of sig 633 terms
      sig_633 <- length(sig_633_terms)
      
      # number of overlapping sig terms
      overlap <- length(intersect(sig_143b_terms, sig_633_terms))
      
      # total terms in 143b
      total_143b <-
        df_combined %>%
        filter(cell_line == "x143b", cluster == x) %>%
        filter(str_detect(id, "GO")) %>%
        nrow()
      
      p <- 
        phyper(
        overlap - 1,
        sig_633,
        total_143b - sig_633,
        sig_143b,
        lower.tail = FALSE,
        log.p = FALSE
      )
      
      results <- c(x, sig_633, overlap, sig_143b, p)
      set_names(results, c("cluster", "x633", "overlap", "x143b", "p_value"))
      # names(results) <- c("x633", "overlap", "x143b", "p_value")
    })

p_results <- 
  df_results %>%
  pivot_longer(-c(cluster, p_value), names_to = "fill", values_to = "num") %>%
  mutate(fill = fct_relevel(fill, "x633", "overlap", "x143b")) %>%
  ggplot(aes(factor(cluster), num)) +
  geom_col(aes(fill = fill), position = "stack") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), aes(fill = fill, label = num)) +
  labs(x = NULL, y = "number of terms (q < 0.005)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
  # geom_text(aes(label = num), position = "fill")

p_sig <-
  df_results %>%
  pivot_longer(-c(cluster, p_value), names_to = "fill", values_to = "num") %>%
  ggplot(aes(x = factor(cluster))) +
  geom_text(aes(label = round(p_value, 5), y = "a")) +
  labs(x = "cluster", y = "p-value") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_results/p_sig + plot_layout(heights = c(5, 1)) & theme_bw()

ggsave("great_compare_barplot.pdf", last_plot(), height = 4, width = 6)
  # scale_color_gradient(low = "black", high = "red")

# # Testing all cluster combos -----------------------------------
# df_results_all <- 
#   expand_grid(1:7, 1:7) %>%
#   pmap(
#     function(...) {
#       # significant terms in 633, per cluster
#       sig_633_terms <-
#         df_combined %>%
#         filter(cluster == ..1) %>%
#         filter(str_detect(id, "GO")) %>%
#         filter(hyper_fdr_q < 0.005, cell_line == "x633") %>%
#         pull(id)
#       
#       # significant terms in 143b, per cluster
#       sig_143b_terms <- 
#         df_combined %>%
#         filter(cluster == ..2) %>% 
#         filter(str_detect(id, "GO")) %>%
#         filter(hyper_fdr_q < 0.005, cell_line == "x143b") %>%
#         pull(id)
#       
#       # number of sig 143b terms
#       sig_143b <- length(sig_143b_terms)
#       
#       # number of sig 633 terms
#       sig_633 <- length(sig_633_terms)
#       
#       # number of overlapping sig terms
#       overlap <- length(intersect(sig_143b_terms, sig_633_terms))
#       
#       # total terms in 143b
#       total_143b <-
#         df_combined %>%
#         filter(cell_line == "x143b", cluster == ..2) %>%
#         filter(str_detect(id, "GO")) %>%
#         nrow()
#       
#       p <- 
#         phyper(
#           overlap - 1,
#           sig_633,
#           total_143b - sig_633,
#           sig_143b,
#           lower.tail = FALSE,
#           log.p = FALSE
#         )
#       
#       data.frame(
#         x633_cluster = ..1,
#         x143b_cluster = ..2,
#         p_value = p
#       )
#       
#       # results <- c(x, sig_633, overlap, sig_143b, p)
#       # set_names(results, c("cluster", "x633", "overlap", "x143b", "p_value"))
#       # names(results) <- c("x633", "overlap", "x143b", "p_value")
#     }) %>% bind_rows()
# 
# ggplot(df_results_all, aes(factor(x633_cluster), factor(x143b_cluster))) +
#   geom_tile(aes(fill = -log10(p_value))) +
#   geom_text(aes(label = round(p_value, 3))) +
#   scale_fill_gradient(
#     low = "white",
#     high = "red",
#     # set upper limit of fill as 0.005 i.e. everything more sig is red
#     limits = c(-log10(1),-log10(0.005)),
#     oob = scales::squish
#   )
# 
# # Testing simplify enrichment -----------------------------------
# test_vec_143b <- 
#  map(
#   1:7,
#   ~great_143b %>%
#     filter(cluster == .x) %>%
#     arrange(hyper_fdr_q) %>%
#     filter(str_detect(id, "GO")) %>%
#     # dplyr::slice(1:50) %>%
#     pull(hyper_fdr_q, name = "id")
# )
#   
# 
# test_vec_633 <- 
#   map(
#     1:7,
#     ~great_633 %>%
#       filter(cluster == .x) %>%
#       arrange(hyper_fdr_q) %>%
#       filter(str_detect(id, "GO")) %>%
#       # dplyr::slice(1:50) %>%
#       pull(hyper_fdr_q, name = "id")
#   )
# 
# list_p <- list(test_vec_633[[1]], test_vec_143b[[1]])
# 
# pdf("633_143b_compare_great_cluster_2_3_combined.pdf", height = 5, width = 10)
# simplifyGOFromMultipleLists(list_p, padj_cutoff = 0.001)
# dev.off()
# 
# simplifyGO(test_vec_633[[1]])