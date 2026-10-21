library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(mclust)

options(future.globals.maxSize = 2500 * 1024^2)

obj <- readRDS("/home/moliva/root_sc_paper/published_datasets/shahan/GSE152766_Root_Atlas.rds")

ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

DefaultAssay(obj) <- "RNA"
obj <- UpdateSeuratObject(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f= obj$orig.ident)
obj <- NormalizeData(obj)
obj <- ScaleData(obj)
obj <- JoinLayers(obj)


label_transfer <- function(ref, query, dims, k.anchor) {

	anchors <- FindTransferAnchors(reference=ref, query=query, dims=1:dims, k.anchor=k.anchor, reference.reduction="integrated.rpca")
 	preds <- TransferData(anchorset=anchors, refdata=ref$identity, dims=1:dims)
	out <- data.frame(cell=rownames(preds), pred=preds$predicted.id, score= preds$prediction.score.max)
	out
	
}

param_grid <- expand.grid(dims=c(30,40,50), k.anchor=c(5, 10, 20))

runs <- lapply(seq_len(nrow(param_grid)), function(i) {
  p <- param_grid[i, ]
  df <- label_transfer(ref=ref, query=obj, dims=p$dims, k.anchor=p$k.anchor)
  df$run_id <- i
  df
})

all_runs <- do.call(rbind, runs)


eligible_cells <- all_runs %>%
  filter(run_id == 3) %>%
  filter(str_detect(pred, "ERS")) %>% 
  pull(cell) %>%
  unique()

df <- all_runs %>%
  filter(cell %in% eligible_cells)
  
df <- df %>%
  mutate(
    state = str_extract(pred, "ERS[+-]"),
    celltype = str_trim(str_remove(pred, "\\s*ERS[+-]\\s*$"))
  )

baseline_type <- df %>%
  filter(run_id == 3) %>%
  select(cell, baseline_celltype = celltype, baseline_state = state)
  
celltype_stability <- df %>%
  inner_join(baseline_type %>% select(cell, baseline_celltype), by = "cell") %>%
  group_by(cell) %>%
  summarise(
    baseline_celltype = first(baseline_celltype),
    celltype_stability = mean(celltype == baseline_celltype),
    .groups = "drop"
  )
  
df2 <- df %>%
  inner_join(baseline_type, by = "cell") %>%
  mutate(in_baseline_type = (celltype == baseline_celltype))

ers_stability_cell <- df2 %>%
  group_by(cell, baseline_celltype) %>%
  summarise(
    n_runs = n(),
    n_in_type = sum(in_baseline_type),
    ers_stability_in_type = ifelse(n_in_type > 0,
                                   mean(state[in_baseline_type] == "ERS+"),
                                   NA_real_),
    .groups = "drop"
  ) %>%
  inner_join(celltype_stability %>% select(cell, celltype_stability), by = "cell")
  
  ers_stability_cell <- ers_stability_cell %>%
  mutate(ers_state_stability = pmax(ers_stability_in_type, 1 - ers_stability_in_type))
  
ers_frac_by_type_run <- df2 %>%
  filter(in_baseline_type) %>%
  group_by(baseline_celltype, run_id) %>%
  summarise(
    n_cells = n(),
    frac_ers_pos = mean(state == "ERS+"),
    .groups = "drop"
  )

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/ERS_fraction_stability.pdf")
ggplot(ers_frac_by_type_run,
             aes(x = baseline_celltype, y = frac_ers_pos)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.3) +
  coord_flip() +
  xlab(NULL) +
  ylab("Proportion") +
  theme_classic()
dev.off()


heat_cells <- df2 %>%
  filter(run_id != 3) %>%
  distinct(cell, baseline_celltype, celltype) %>% 
  count(baseline_celltype, celltype, name = "n_cells")

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/heatmap_celltype_flips.pdf")
ggplot(heat_cells,
       aes(x = celltype, y = baseline_celltype, fill = n_cells)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10") +
  theme_classic() +
  xlab("Assigned cell type (other runs)") +
  ylab("Cell type in paper run") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/stability_score_distribution.pdf")
ggplot(ers_stability_cell, aes(x = ers_state_stability)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ baseline_celltype, scales = "free_y") +
  xlab("Per-cell ERS state stability within baseline cell type") +
  ylab("Number of cells") +
  theme_classic()
dev.off()


pred_wide <- df %>%
  select(cell, run_id, pred) %>%
  distinct() %>%
  pivot_wider(names_from = run_id, values_from = pred)

run_ids <- sort(unique(df$run_id))
ari_mat <- matrix(NA_real_, length(run_ids), length(run_ids),
                  dimnames = list(run_ids, run_ids))

for (i in seq_along(run_ids)) {
  for (j in seq_along(run_ids)) {
    ri <- as.character(run_ids[i])
    rj <- as.character(run_ids[j])
    ari_mat[i, j] <- adjustedRandIndex(pred_wide[[ri]], pred_wide[[rj]])
  }
}

ari_df <- as.data.frame(as.table(ari_mat))
colnames(ari_df) <- c("run_i", "run_j", "ARI")

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/ARI_heatmap.pdf")
ggplot(ari_df, aes(x = run_i, y = run_j, fill = ARI)) +
  geom_tile() +
  theme_classic() +
  xlab("Run") + ylab("Run")
dev.off()


celltype_stability_run <- df2 %>%
  group_by(run_id) %>%
  summarise(
    frac_same_celltype = mean(celltype == baseline_celltype),
    .groups = "drop"
  )
  
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/celltype_stability_per_run.pdf")
ggplot(celltype_stability_run,
       aes(x = factor(run_id), y = frac_same_celltype)) +
  geom_point(size = 3) +
  geom_line(group = 1) +
  geom_vline(xintercept = 3, linetype = 2) +
  ylim(0,1) +
  ylab("Fraction of cells retaining paper-run cell type") +
  xlab("Parameter set (run)") +
  theme_classic()
dev.off()

ers_switch_summary <- df2 %>%
  filter(run_id != 3) %>%
  mutate(
    switch_category = case_when(
      is.na(state) ~ "non-ERS cell type",
      state == baseline_state ~ "same ERS state",
      TRUE ~ "ERS state flipped"
    )
  ) %>%
  count(run_id, switch_category) %>%
  group_by(run_id) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()
  
ers_switch_summary$switch_category <- factor(
  ers_switch_summary$switch_category,
  levels = c("same ERS state",
             "ERS state flipped",
             "non-ERS cell type")
)

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/ERS_switch_per_run.pdf")
ggplot(ers_switch_summary,
       aes(x = switch_category, y = frac,
           group = run_id, colour = factor(run_id))) +
  geom_point(size = 3) +
  geom_line(alpha = 0.4) +
  ylim(0,1) +
  theme_classic()
dev.off()