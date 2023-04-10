library(tidyverse)

motif_output_633 <- read_tsv("nat-comm_revisions/recluster_maelstrom_output/final.out.txt")
motif_output_143b <- read_tsv("nat-comm_revisions/143b_gimme-maelstrom_nostable_reorder_output/final.out.txt")

colnames(motif_output_633) <- c("factor", paste0("cluster", 1:8))
colnames(motif_output_143b) <- c("factor", paste0("cluster", 1:7))

motif_output_633_long <- 
  motif_output_633 %>%
  pivot_longer(-factor, names_to = "cluster", values_to = "enrichment_633")

motif_output_143b_long <-
  motif_output_143b %>%
  pivot_longer(-factor, names_to = "cluster", values_to = "enrichment_143b")

df_combined <- 
  motif_output_633_long %>%
  left_join(motif_output_143b_long, by = c("factor", "cluster")) 
  # filter(cluster != "cluster8")

df_correlate_clusters <- 
  df_combined %>%
  filter(cluster != "cluster8") %>%
  pivot_longer(-c(factor, cluster), names_pattern = "enrichment_(.*)") %>%
  unite("id", c("cluster", "name")) %>%
  filter(abs(value) > 2) %>%
  pivot_wider(names_from = id, values_from = value) %>%
  select(-factor) %>%
  correlate()

df_correlate_clusters %>%
  pivot_longer(-term, names_to = "term2", values_to = "r") %>%
  separate("term", c("cluster1", "cellline1")) %>%
  separate("term2", c("cluster2", "cellline2")) %>%
  filter(cellline1 == "633", cellline2 == "143b") %>%
  unite("term1", c("cluster1", "cellline1")) %>%
  unite("term2", c("cluster2", "cellline2")) %>%
  ggplot(aes(term1, term2)) +
  geom_tile(aes(fill = r)) +
  geom_text(aes(label = round(r, 3))) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  coord_fixed() +
  # theme_bw() +
  theme(
    panel.grid = element_blank()
  )

df_summary_stats <- 
  df_combined %>%
  mutate(sig = case_when((abs(enrichment_633) > 2 &
         abs(enrichment_143b) > 2) ~ "sig",
         TRUE ~ "non_sig")) %>%
  group_by(cluster, sig) %>%
  summarize(cor = round(cor(enrichment_633, enrichment_143b), digits = 3)) %>%
  ungroup() %>%
  pivot_wider(names_from = sig, values_from = cor)

p_motif_comparison <- 
  df_combined %>%
  filter(cluster != "cluster8") %>%
  mutate(sig = case_when((abs(enrichment_633) > 2 &
                            abs(enrichment_143b) > 2) ~ "yes",
                         TRUE ~ "no")) %>%
  ggplot(aes(enrichment_633, enrichment_143b)) +
  geom_point(aes(color = sig)) +
  geom_abline(slope = 1, color = "red") +
  geom_text(data = df_summary_stats %>%
              filter(cluster != "cluster8"),
            aes(x = 8, y = -5, label = paste0("r = ", non_sig))) +
  geom_text(data = df_summary_stats %>%
              filter(cluster != "cluster8"),
            aes(x = 8, y = -7, label = paste0("r = ", sig)), color = "red") +
  scale_color_manual(values = c("black", "red")) +
  coord_fixed() +
  labs(x = "633 motif enrichment",
       y = "143b motif enrichment", 
       color = "Threshold > 2") +
  facet_grid(.~cluster) 

ggsave("143b_633_maelstrom_compare.pdf", p_motif_comparison, width = 15)
