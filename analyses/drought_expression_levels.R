library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

drought_100 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_100percentH2O/integrated_drought_100percentH2O_predicted_id.rds")
drought_20 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/integrated_drought_20percentH2O_predicted_id.rds")

drought_100@meta.data$water <- 100
drought_20@meta.data$water <- 20
df <- rbind(drought_100@meta.data[,c("nCount_RNA","nFeature_RNA","predicted.identity","water","orig.ident")],drought_20@meta.data[,c("nCount_RNA","nFeature_RNA","predicted.identity","water","orig.ident")])
df$state <- df$predicted.identity
df$state <- gsub(".*(ERS[+-])$", "\\1", df$state)
df$state[!grepl("ERS[+-]$", df$state)] <- "non-ERS"
df$state <- factor(df$state)
df$water <- factor(df$water)

plot_resp <- function(df, variable, outfile) {

  rep_summary <- df %>%
    group_by(state, water, orig.ident) %>%
    summarise(median_val = median(.data[[variable]]), .groups = "drop")

  cond_summary <- rep_summary %>%
    group_by(state, water) %>%
    summarise(cond_median = median(median_val), .groups="drop")

  effect_df <- cond_summary %>%
    pivot_wider(names_from = water, values_from = cond_median) %>%
    mutate(FC = `20` / `100`)

  label_pos <- df %>%
    group_by(state) %>%
    summarise(y = max(.data[[variable]]) * 1.1)

  effect_df <- left_join(effect_df, label_pos, by = "state")

  dodge_width <- 0.9

  pdf(outfile)
  print(
    ggplot(df, aes(x = state, y = .data[[variable]], fill = water)) +
      geom_violin(position = position_dodge(dodge_width), scale = "width", alpha = 0.5) +
      geom_boxplot(width = 0.2, position = position_dodge(dodge_width), outlier.shape = NA) +
      geom_point(data = rep_summary,
                 aes(x = state, y = median_val),
                 position = position_dodge(dodge_width),
                 size = 3, colour = "black") +
      geom_text(data = effect_df,
                aes(x = state, y = y, label = paste0("FC=", round(FC,2))),
                inherit.aes = FALSE) +
      theme_classic() +
      labs(x = "State", y = variable, fill = "Water")
  )
  dev.off()
}

plot_resp(df,
        "nCount_RNA",
        "/home/moliva/root_sc_paper/analyses/drought/RNA_per_state_water.pdf")

plot_resp(df,
        "nFeature_RNA",
        "/home/moliva/root_sc_paper/analyses/drought/Features_per_state_water.pdf")