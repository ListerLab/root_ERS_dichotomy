library(Seurat)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(splines)




drought_100 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_100percentH2O/integrated_drought_100percentH2O_predicted_id.rds")
drought_20 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/integrated_drought_20percentH2O_predicted_id.rds")

modules <- read.csv("/home/moliva/root_sc_paper/analyses/hdWGCNA/modules.tsv", sep='\t')
mods <- list(subset(modules, module == "Module 20")$gene_name, subset(modules, module == "Module 21")$gene_name)
drought_100 <- AddModuleScore(drought_100, features=mods)
drought_20 <- AddModuleScore(drought_20, features=mods)

drought_100@meta.data$water <- 100
drought_20@meta.data$water <- 20
df <- rbind(drought_100@meta.data[,c("nCount_RNA","predicted.identity","water","orig.ident","Cluster1","Cluster2")],drought_20@meta.data[,c("nCount_RNA","predicted.identity","water","orig.ident","Cluster1","Cluster2")])
df$water <- relevel(factor(df$water), ref = "100")
df$logUMI <- log10(df$nCount_RNA +1)

df <- df[grep("ERS",df$predicted.identity),]
df$state <- df$predicted.identity
df$state <- gsub(".*(ERS[+-])$", "\\1", df$state)
df$state <- factor(df$state)


m20_score <- lmer(
  Cluster1 ~ water + ns(logUMI, df = 4) + (1 | orig.ident),
  data = df
)

m21_score <- lmer(
  Cluster2 ~ water + ns(logUMI, df = 4) + (1 | orig.ident),
  data = df
)

make_binned_score_plot <- function(df, module_col, module_label, model, n_bins = 10) {
  
  stat <- broom.mixed::tidy(model, effects = "fixed") %>%
    filter(term == "water20")
  
  pval <- stat$p.value
  beta <- stat$estimate
  
  label <- paste0(
    "Adjusted for logUMI\n",
    "beta = ", signif(beta, 3),
    ", P = ", format.pval(pval, digits = 2, eps = 1e-4)
  )
  
  df_bins <- df %>%
    mutate(UMI_bin = ntile(logUMI, n_bins)) %>%
    group_by(water, UMI_bin) %>%
    summarise(
      median_logUMI = median(logUMI, na.rm = TRUE),
      mean_score = mean(.data[[module_col]], na.rm = TRUE),
      se_score = sd(.data[[module_col]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  ggplot(df_bins, aes(x = median_logUMI, y = mean_score,
                      colour = water, group = water)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = mean_score - se_score,
                      ymax = mean_score + se_score),
                  width = 0) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = label,
      hjust = 1.05,
      vjust = 1.3,
      size = 3
    ) +
    labs(
      x = "Median log10 UMI count",
      y = paste0("Mean ERS ", module_label, " score"),
      colour = "Water content"
    ) +
    theme_classic()
}

p20 <- make_binned_score_plot(df, "Cluster1", "module 20", m20_score)
p21 <- make_binned_score_plot(df, "Cluster2", "module 21", m21_score)

pdf("/home/moliva/root_sc_paper/analyses/drought/binned_score_module20.pdf")
p20
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/drought/binned_score_module21.pdf")
p21
dev.off()

m_state_unadj <- glmer(
  state ~ water + (1 | orig.ident),
  data = df,
  family = binomial
)

m_state_adj <- glmer(
  state ~ water + logUMI + (1 | orig.ident),
  data = df,
  family = binomial
)

state_stats <- bind_rows(
  broom.mixed::tidy(
    m_state_unadj,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(term == "water20") %>%
    mutate(model = "Unadjusted"),

  broom.mixed::tidy(
    m_state_adj,
    effects = "fixed",
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    filter(term == "water20") %>%
    mutate(model = "Adjusted for logUMI")
) %>%
  mutate(
    model = factor(model, levels = c("Unadjusted", "Adjusted for logUMI")),
    label = paste0(
      "OR = ", signif(estimate, 2),
      " [", signif(conf.low, 2), ", ", signif(conf.high, 2), "]\n",
      "P = ", format.pval(p.value, digits = 2, eps = 1e-4)
    )
  )

pdf("/home/moliva/root_sc_paper/analyses/drought/state_odds.pdf")
ggplot(state_stats, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), linewidth = 0.6) +
  geom_text(
    aes(label = label),
    nudge_y = 0.18,
    hjust = 0,
    size = 3.2
  ) +
  scale_x_continuous(
    limits = c(0.8, 5),
    breaks = pretty(c(1, state_stats$conf.high))
  ) +
  labs(
    x = "Odds ratio for ERS+ identity\n20% water versus 100% water",
    y = NULL
  ) +
  theme_classic()
dev.off()
