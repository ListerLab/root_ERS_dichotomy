library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)


obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")

obj@reductions[["umap"]] <- obj@reductions[["umap.rpca"]]
cds <- as.cell_data_set(obj)
cds <- cluster_cells(cds, resolution=0.0002)
cds <- learn_graph(cds, learn_graph_control=list(minimal_branch_len=20))

pdf("/home/moliva/root_sc_paper/analyses/QC/trajectories/monocle/monocle_trajectories_partitions.pdf")
plot_cells(cds, color_cells_by="partition",label_leaves=FALSE,label_branch_points=FALSE,label_groups_by_cluster=FALSE)
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/trajectories/monocle/monocle_trajectories_identities.pdf")
plot_cells(cds, color_cells_by="identity", label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
dev.off()

root_cells <- c()

for(i in unique(obj$lineage)){
	meta <- subset(obj@meta.data, lineage == i)
	cell <- row.names(meta[meta$early_epi_genes_score==max(meta$early_epi_genes_score),])
	root_cells <- c(root_cells,cell)
}

cds <- order_cells(cds, root_cells=root_cells)

pdf("/home/moliva/root_sc_paper/analyses/QC/trajectories/monocle/monocle_trajectories_pseudotime.pdf")
plot_cells(cds, color_cells_by="pseudotime",label_roots=FALSE,label_groups_by_cluster=FALSE,label_leaves=FALSE,label_branch_points=FALSE)
dev.off()


mono_pt <- pseudotime(cds)
obj$monocle_pt <- mono_pt[colnames(obj)]

rescale <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (any(!is.finite(rng)) || diff(rng) == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

obj$monocle_pt_scaled <- ave(obj$monocle_pt, obj$lineage, FUN = rescale)

cor_df <- obj@meta.data %>%
  group_by(lineage) %>%
  summarise(
    rho = cor(pseudotime, monocle_pt_scaled, method = "spearman", use = "complete.obs"),
    n = n(),
    x = 0.05,
    y = 0.95,
    label = paste0("rho = ", round(rho, 2), "\n", "n = ", n),
    .groups = "drop"
  )

pdf("/home/moliva/root_sc_paper/analyses/QC/trajectories/monocle/scatter_corr_pseudotime.pdf")
ggplot(obj@meta.data, aes(x = pseudotime, y = monocle_pt_scaled)) +
  geom_point(size = 0.3, alpha = 0.4) +
  facet_wrap(~ lineage) +
  geom_text(
    data = cor_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  labs(
    x = "PAGA pseudotime",
    y = "Monocle pseudotime (scaled 0–1)",
    title = "Pseudotime concordance between PAGA and Monocle"
  ) +
  theme_classic()
dev.off()

