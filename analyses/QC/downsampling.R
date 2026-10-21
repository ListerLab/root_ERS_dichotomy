library(Seurat)
library(ggplot2)
library(Matrix)
library(DropletUtils)

obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj), value=TRUE)
chloro <- grep("ATC", row.names(obj), value=TRUE)


obj$cell_state <- ifelse(grepl("ERS\\+$", obj$identity), "ERS+",
                  ifelse(grepl("ERS\\-$", obj$identity), "ERS-", NA))

obj$cell_type <- sub(" ERS[+-]$", "", obj$identity)


downsample_ERSpos_by_celltype_quantile <- function(obj, seed) {
  set.seed(seed)

  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  meta <- obj@meta.data
  meta <- meta[colnames(counts), , drop = FALSE]

  all_cells <- colnames(counts)
  all_cell_types <- unique(na.omit(meta$cell_type))

  new_counts <- counts
  original_depths <- Matrix::colSums(counts)
  target_depths_all <- original_depths
  matched_group <- rep(NA_character_, length(all_cells))
  names(matched_group) <- all_cells

  summary_list <- list()

  for (i in seq_along(all_cell_types)) {
  ct <- all_cell_types[i]
  t0 <- Sys.time()

  cat(sprintf("[%d/%d] %s ... ", i, length(all_cell_types), ct))
  flush.console()

  cells_ct <- all_cells[which(meta[all_cells, "cell_type"] == ct)]

  pos_cells <- cells_ct[which(meta[cells_ct, "cell_state"] == "ERS+")]
  neg_cells <- cells_ct[which(meta[cells_ct, "cell_state"] == "ERS-")]

  if (length(pos_cells) == 0 || length(neg_cells) == 0) {
    matched_group[cells_ct] <- paste0(ct, "_unchanged_missing_one_state")

    summary_list[[length(summary_list) + 1]] <- data.frame(
      cell_type = ct,
      n_pos = length(pos_cells),
      n_neg = length(neg_cells),
      median_pos_before = if (length(pos_cells) > 0) median(original_depths[pos_cells]) else NA,
      median_neg_before = if (length(neg_cells) > 0) median(original_depths[neg_cells]) else NA,
      median_pos_target = NA,
      status = "unchanged_missing_group"
    )

    cat(sprintf("skipped (%.1f s)\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    flush.console()
    next
  }

  pos_depths <- original_depths[pos_cells]
  neg_depths <- original_depths[neg_cells]

  ranks <- rank(pos_depths, ties.method = "random") / length(pos_depths)
  target_depths <- quantile(neg_depths, probs = ranks, type = 8)
  target_depths <- as.integer(round(target_depths))
  target_depths <- pmin(pos_depths, target_depths)
  names(target_depths) <- pos_cells

  props <- target_depths / pos_depths
  props[props > 1] <- 1
  props[is.na(props)] <- 1

  pos_mat <- counts[, pos_cells, drop = FALSE]
  pos_mat_ds <- DropletUtils::downsampleMatrix(pos_mat, prop = props, bycol = TRUE)

  new_counts[, pos_cells] <- pos_mat_ds

  target_depths_all[pos_cells] <- target_depths
  matched_group[pos_cells] <- paste0(ct, "_ERS+_matched_to_ERS-")
  matched_group[neg_cells] <- paste0(ct, "_ERS-_reference")

  summary_list[[length(summary_list) + 1]] <- data.frame(
    cell_type = ct,
    n_pos = length(pos_cells),
    n_neg = length(neg_cells),
    median_pos_before = median(pos_depths),
    median_neg_before = median(neg_depths),
    median_pos_target = median(target_depths),
    status = "processed"
  )

  cat(sprintf("done (%.1f s)\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  flush.console()
  }

  new_obj <- CreateSeuratObject(counts = new_counts, meta.data = meta)
  new_obj$cell_type <- meta$cell_type
  new_obj$cell_state <- meta$cell_state
  new_obj$identity <- meta$identity

  new_obj$original_nCount_RNA <- original_depths[colnames(new_obj)]
  new_obj$target_nCount_RNA <- target_depths_all[colnames(new_obj)]
  new_obj$matched_nCount_RNA <- Matrix::colSums(new_counts)
  new_obj$downsampling_group <- matched_group[colnames(new_obj)]

  summary_df <- do.call(rbind, summary_list)

  list(obj = new_obj, summary = summary_df)
}

down1 <- downsample_ERSpos_by_celltype_quantile(obj=obj, seed=1)


process_obj <- function(x){
	x <- NormalizeData(x, verbose=FALSE)
	x <- FindVariableFeatures(x, nfeatures=25000, verbose=FALSE)
	VariableFeatures(x) <- setdiff(VariableFeatures(x),c(proto$V1,mito,chloro))
	x <- ScaleData(x, verbose = FALSE)
	x <- RunPCA(x, npcs=50, verbose=FALSE)
	x <- RunUMAP(x, dims=1:50, verbose=FALSE)
	x$identity <- obj@meta.data[row.names(x@meta.data),"identity"]
	return(x)
}

down1$obj <- process_obj(down1$obj)

pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/umap_down1.pdf",width=9, height=6)
DimPlot(down1$obj,group.by="identity")
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/UMI_per_identity_down1.pdf")
ggplot(down1$obj@meta.data,aes(x=identity,y=nCount_RNA))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
early_genes <- readRDS("/home/moliva/root_sc_paper/analyses/trajectories/early_genes.rds")
down1$obj <- AddModuleScore(down1$obj, features=list(early_genes), name="score")
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/early_genes_score_down1.pdf")
FeaturePlot(down1$obj,"score1", min.cutoff=0, max.cutoff="q50",order=TRUE)
dev.off()



down2 <- downsample_ERSpos_by_celltype_quantile(obj=obj, seed=2)
down2$obj <- process_obj(down2$obj)
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/umap_down2.pdf",width=9, height=6)
DimPlot(down2$obj,group.by="identity")
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/UMI_per_identity_down2.pdf")
ggplot(down2$obj@meta.data,aes(x=identity,y=nCount_RNA))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
early_genes <- readRDS("/home/moliva/root_sc_paper/analyses/trajectories/early_genes.rds")
down2$obj <- AddModuleScore(down2$obj, features=list(early_genes), name="score")
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/early_genes_score_down2.pdf")
FeaturePlot(down2$obj,"score1", min.cutoff=0, max.cutoff="q50",order=TRUE)
dev.off()

down3 <- downsample_ERSpos_by_celltype_quantile(obj=obj, seed=3)
down3$obj <- process_obj(down3$obj)
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/umap_down3.pdf",width=9, height=6)
DimPlot(down3$obj,group.by="identity")
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/UMI_per_identity_down3.pdf")
ggplot(down3$obj@meta.data,aes(x=identity,y=nCount_RNA))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
early_genes <- readRDS("/home/moliva/root_sc_paper/analyses/trajectories/early_genes.rds")
down3$obj <- AddModuleScore(down3$obj, features=list(early_genes), name="score")
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/early_genes_score_down3.pdf")
FeaturePlot(down3$obj,"score1", min.cutoff=0, max.cutoff="q50",order=TRUE)
dev.off()

saveRDS(down1,"/home/moliva/root_sc_paper/analyses/QC/downsampling/down1.rds")
saveRDS(down2,"/home/moliva/root_sc_paper/analyses/QC/downsampling/down2.rds")
saveRDS(down3,"/home/moliva/root_sc_paper/analyses/QC/downsampling/down3.rds")
