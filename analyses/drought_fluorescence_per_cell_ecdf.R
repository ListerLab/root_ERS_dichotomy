library(dplyr)
library(ggplot2)


df <- read.csv("/home/moliva/root_sc_paper/analyses/drought/per_cell_segmentation_fluorecence_quants.csv")

df <- df %>%
  mutate(
    batch = factor(batch),
    water = factor(water, levels = c(100, 20)),
    root_id = interaction(gene, batch, root, drop = TRUE)
  )
  
control_scale <- df %>%
  filter(water == "100") %>%
  group_by(gene, batch, root_id) %>%
  summarise(
    control_root_mean = mean(mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(gene, batch) %>%
  summarise(
    control_scale = mean(control_root_mean, na.rm = TRUE),
    .groups = "drop"
  )

df_norm <- df %>%
  left_join(control_scale, by = c("gene", "batch")) %>%
  mutate(norm = mean / control_scale)

p_ecdf <- ggplot(
  df_norm,
  aes(
    x = norm,
    colour = batch,
    linetype = water,
    group = interaction(batch, root_id, water)
  )
) +
  stat_ecdf(linewidth = 0.8) +
  facet_wrap(~ gene, scales = "free_x") +
  scale_linetype_manual(
    values = c("100" = "solid", "20" = "dashed"),
    labels = c("100%", "20%")
  ) +
  labs(
    x = "Normalised fluorescence intensity",
    y = "Cumulative fraction of segmented cells",
    colour = "Batch",
    linetype = "Water content"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "right"
  )

pdf("/home/moliva/root_sc_paper/analyses/drought/percell_quant_ecdf.pdf", height=4, width=10)
p_ecdf
dev.off()


