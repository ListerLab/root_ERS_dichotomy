library(reticulate)
library(hdf5r)
library(velocyto.R)
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(velociraptor)
library(ggplot2)
library(grid)
library(dplyr)



sample_names <- c("RL650_col_cvi_r3","RL970_col_ler_r4","RL971_col_c24_r5","RL2757_col_fullMS_r1","RL2758_col_fullMS_r2","RL2759_col_halfMS_r1","RL2760_col_halfMS_r2","RL4207_col_fullMS_r3","RL4208_col_halfMS_r3")

spliced_list   <- list()
unspliced_list <- list()

for (x in seq_along(sample_names)) {

  s <- sample_names[x]

  loom_path <- file.path("/home/moliva/root_sc_paper/data/raw_data_per_library", s, paste0("cellranger_", s), "velocyto", paste0("cellranger_",s, ".loom"))

  ldat <- read.loom.matrices(loom_path)

  sp <- ldat$spliced
  un <- ldat$unspliced

  f <- H5File$new(loom_path)
  gene_ids <- f[["row_attrs"]][["Accession"]][]
  f$close_all()
  rownames(sp) <- gene_ids
  rownames(un) <- gene_ids

  colnames(sp) <- gsub(".*:|x$", "",colnames(sp))
  colnames(un) <- gsub(".*:|x$", "",colnames(un))
  colnames(sp) <- paste0(colnames(sp),"-1_",x)
  colnames(un) <- paste0(colnames(un),"-1_",x)

  spliced_list[[s]]   <- sp
  unspliced_list[[s]] <- un

  rm(ldat, sp, un)
  gc()
}


common_genes <- Reduce(intersect, lapply(spliced_list, rownames))

spliced_list   <- lapply(spliced_list,   function(m) m[common_genes, , drop=FALSE])
unspliced_list <- lapply(unspliced_list, function(m) m[common_genes, , drop=FALSE])

spliced_all   <- do.call(cbind, spliced_list)
unspliced_all <- do.call(cbind, unspliced_list)

obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")

obj$latent_time <- NA_real_

run_scvelo_lineage <- function(sce, genes, root_cells) {
  sce$my_root <- colnames(sce) %in% root_cells

  velociraptor::scvelo(
    sce,
    subset.row = genes,
    mode = "dynamical",
    use.dimred="RPCA",
    scvelo.params = list(
      neighbors = list(n_neighbors=30L, n_pcs=50L),
      recover_dynamics = list(),
      latent_time = list(root_key="my_root")
    )
  )
}

Idents(obj) <- "lineage"

stats <- data.frame(lineage=character(), n_genes=integer())

min_cells_un <- 50
min_total_un <- 300
keep_un_global <- (rowSums(unspliced_all > 0) >= min_cells_un) & (rowSums(unspliced_all) >= min_total_un)

min_cells_sp <- 50
min_total_sp <- 1000
keep_sp_global <- (rowSums(spliced_all > 0) >= min_cells_sp) & (rowSums(spliced_all) >= min_total_sp)

genes_ok_global <- intersect(rownames(unspliced_all)[keep_un_global], rownames(spliced_all)[keep_sp_global])

arrows_all <- list()

for (lin in unique(obj$lineage)) {

  obj_sub <- subset(obj, idents=lin)
  
  obj_sub <- FindVariableFeatures(obj_sub, selection.method="vst", nfeatures=4000)
  proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
  mito <- grep("ATM", row.names(obj_sub), value=TRUE)
  chloro <- grep("ATC", row.names(obj_sub), value=TRUE)
  VariableFeatures(obj_sub) <- setdiff(VariableFeatures(obj_sub),c(proto$V1,mito,chloro))
  genes <- VariableFeatures(obj_sub)
  
  genes <- intersect(genes, genes_ok_global)
  
  stats <- rbind(stats, data.frame(lineage=lin, n_genes=length(genes)))

  sce <- SingleCellExperiment(
  	assays=list(
  		counts=GetAssayData(obj_sub,assay="RNA",slot="counts"), 
  		spliced=spliced_all[row.names(obj_sub),colnames(obj_sub)], 
  		unspliced=unspliced_all[row.names(obj_sub),colnames(obj_sub)]
  	)
  )
  
  rpca <- Embeddings(obj_sub, "integrated.rpca")
  rpca <- rpca[colnames(obj_sub), , drop=FALSE]
  reducedDim(sce, "RPCA") <- rpca
  
  scores <- obj_sub$early_epi_genes_score
  n <- ncol(obj_sub)
  root_cells <- names(sort(scores, decreasing=TRUE))[1:as.integer(round(0.01 * n))]

  sce <- run_scvelo_lineage(sce = sce, genes = genes, root_cells = root_cells)
  
  reducedDim(sce, "UMAP") <- Embeddings(obj_sub, "umap.rpca")[colnames(sce), 1:2]
  vel_umap <- embedVelocity(sce, sce, use.dimred="UMAP")
  arrows_df <- gridVectors(sce, embedded=vel_umap, use.dimred="UMAP", resolution=40, scale=TRUE, as.data.frame=TRUE)
  arrows_df$lineage <- lin
  arrows_all[[lin]] <- arrows_df
  
  obj$latent_time[colnames(sce)] <- sce$latent_time
  
}

arrows_all_df <- do.call(rbind, arrows_all)

pdf("/home/moliva/root_sc_paper/analyses/QC/trajectories/scvelo/umap_latent_time.pdf")
FeaturePlot(obj, features = "latent_time")
dev.off()

saveRDS(arrows_all_df,"/home/moliva/root_sc_paper/analyses/QC/trajectories/scvelo/arrows_all_df.rds")
saveRDS(obj@meta.data,"/home/moliva/root_sc_paper/analyses/QC/trajectories/scvelo/metadata.rds")

cor_df <- obj@meta.data %>%
  group_by(lineage) %>%
  summarise(
    rho = cor(pseudotime, latent_time, method = "spearman", use = "complete.obs"),
    n = n(),
    x = 0.05,
    y = 0.95,
    label = paste0("rho = ", round(rho, 2), "\n", "n = ", n),
    .groups = "drop"
  )

pdf("/home/moliva/root_sc_paper/analyses/QC/trajectories/scvelo/scatter_corr_pseudotime.pdf")
ggplot(obj@meta.data, aes(x = pseudotime, y = latent_time)) +
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
    x = "Diffusion pseudotime",
    y = "scvelo latent time",
    title = "Pseudotime concordance between dpt and velocity"
  ) +
  theme_classic()
dev.off()

