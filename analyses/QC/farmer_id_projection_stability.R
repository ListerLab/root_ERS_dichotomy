library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(mclust)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep1/filtered_feature_bc_matrix/")
rna1 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep2/filtered_feature_bc_matrix/")
rna2 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep3/filtered_feature_bc_matrix/")
rna3 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep4/filtered_feature_bc_matrix/")
rna4 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep5/filtered_feature_bc_matrix/")
rna5 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

obj.rna <- merge(rna1, list(rna2,rna3,rna4,rna5))

obj.rna <- NormalizeData(obj.rna)
obj.rna <- ScaleData(obj.rna)

obj.rna <- FindVariableFeatures(obj.rna, nfeatures=25000,verbose=F)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj.rna), value=TRUE)
chloro <- grep("ATC", row.names(obj.rna), value=TRUE)
VariableFeatures(obj.rna) <- setdiff(VariableFeatures(obj.rna),c(proto$V1,mito,chloro))

obj.rna <- RunPCA(obj.rna, npcs=50,verbose=F)

obj.rna <- IntegrateLayers(obj.rna, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50, verbose=FALSE)
obj.rna <- JoinLayers(obj.rna)

ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")



label_transfer <- function(ref, query, dims, k.anchor) {

	anchors <- FindTransferAnchors(reference=ref, query=query, dims=1:dims, k.anchor=k.anchor, reference.reduction="integrated.rpca")
 	preds <- TransferData(anchorset=anchors, refdata=ref$identity, dims=1:dims)
	out <- data.frame(cell=rownames(preds), pred=preds$predicted.id, score= preds$prediction.score.max)
	out
	
}

param_grid <- expand.grid(dims=c(30,40,50), k.anchor=c(5, 10, 20))

runs <- lapply(seq_len(nrow(param_grid)), function(i) {
  p <- param_grid[i, ]
  df <- label_transfer(ref=ref, query=obj.rna, dims=p$dims, k.anchor=p$k.anchor)
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
  
df2 <- df %>%
  inner_join(baseline_type, by = "cell") %>%
  mutate(in_baseline_type = (celltype == baseline_celltype))

ers_frac_by_type_run <- df2 %>%
  filter(in_baseline_type) %>%
  group_by(baseline_celltype, run_id) %>%
  summarise(
    n_cells = n(),
    frac_ers_pos = mean(state == "ERS+"),
    .groups = "drop"
  )

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/ERS_fraction_stability.pdf")
ggplot(ers_frac_by_type_run,
             aes(x = baseline_celltype, y = frac_ers_pos)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.3) +
  coord_flip() +
  xlab(NULL) +
  ylab("Proportion") +
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

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/ARI_heatmap.pdf")
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
  
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/celltype_stability_per_run.pdf")
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

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/ERS_switch_per_run.pdf")
ggplot(ers_switch_summary,
       aes(x = switch_category, y = frac,
           group = run_id, colour = factor(run_id))) +
  geom_point(size = 3) +
  geom_line(alpha = 0.4) +
  ylim(0,1) +
  theme_classic()
dev.off()