library(Seurat)
options(future.globals.maxSize = 2500 * 1024^2)

RL5314 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/seurat_object/RL5314_drought_20percentH2O_r1_singlets.rds")
RL5314$orig.ident <- "RL5314_drought_20percentH2O_r1"

RL5537 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5537_drought_20percentH2O_r2/seurat_object/RL5537_drought_20percentH2O_r2_singlets.rds")
RL5537$orig.ident <- "RL5537_drought_20percentH2O_r2"

obj <- merge(RL5314, RL5537)

obj <- FindVariableFeatures(obj, nfeatures=25000,verbose=F)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj), value=TRUE)
chloro <- grep("ATC", row.names(obj), value=TRUE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))

obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=50,verbose=F)

obj <- IntegrateLayers(obj, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50, verbose=FALSE)

obj <- RunUMAP(obj, reduction="integrated.rpca", dims=1:50, reduction.name="umap.rpca",verbose=F)
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:50, k.param=50, verbose=F)
obj <- FindClusters(obj, resolution = 0.3, cluster.name = "clusters", verbose=F)

pdf("/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/obj_clusters.pdf", width=6, height=6)
DimPlot(obj, group.by="clusters", label=TRUE)+NoLegend()
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/obj_umi.pdf", width=6, height=6)
FeaturePlot(obj, "nCount_RNA", min.cutoff="q10", max.cutoff="q90")
dev.off()

#One clusters of low-UMI "cells" that is labelled as epidermis and clusters with root cap will be removed

obj <- subset(obj, idents=11, invert=TRUE)
obj <- FindVariableFeatures(obj, nfeatures=25000,verbose=F)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=50,verbose=F)
obj <- IntegrateLayers(obj, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50, verbose=FALSE)
obj <- RunUMAP(obj, reduction="integrated.rpca", dims=1:50, reduction.name="umap.rpca",verbose=F)


ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")
obj <- JoinLayers(obj)
anchors <- FindTransferAnchors(reference=ref, query=obj, dims=1:50, reference.reduction="integrated.rpca")
predictions_lineage <- TransferData(anchorset=anchors, refdata=ref$lineage, dims=1:50)
obj <-AddMetaData(obj, metadata=predictions_lineage[,1,drop=FALSE], col.name="predicted.lineage")
predictions_identity <- TransferData(anchorset=anchors, refdata=ref$identity, dims=1:50)
obj <-AddMetaData(obj, metadata=predictions_identity[,1,drop=FALSE], col.name="predicted.identity")

pdf("/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/umap_drought_20percentH2O_lineage.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.lineage")
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/umap_drought_20percentH2O_identity.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.identity")
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/umap_drought_20percentH2O_orig.ident.pdf", width=9, height=6)
DimPlot(obj, group.by="orig.ident")
dev.off()

obj$clusters <- NULL

saveRDS(obj,"/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/integrated_drought_20percentH2O_predicted_id.rds")
