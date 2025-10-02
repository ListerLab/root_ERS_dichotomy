library(Seurat)
library(SeuratWrappers)
library(monocle3)

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
