library(Seurat)
options(future.globals.maxSize = 2500 * 1024^2)

RL971 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL971_col_c24_r5/seurat_object/RL971_col_c24_r5_singlets_genotype.rds")
Idents(RL971) <- "genotype"
RL971 <- subset(RL971,idents="c24")
RL971$orig.ident <- "RL971_col_c24"

RL1465 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL1465_c24_N681600/seurat_object/RL1465_c24_N681600_singlets_genotype.rds")
Idents(RL1465) <- "genotype"
RL1465 <- subset(RL1465,idents="c24")
RL1465$orig.ident <- "RL1465_c24"

RL4459 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4459_c24_pyl11458/seurat_object/RL4459_c24_pyl11458_singlets_genotype.rds")
Idents(RL4459) <- "genotype"
RL4459 <- subset(RL4459,idents="c24")
RL4459$orig.ident <- "RL4459_col_c24"

obj <- merge(RL971, list(RL1465, RL4459))

obj <- FindVariableFeatures(obj, nfeatures=25000,verbose=F)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj), value=TRUE)
chloro <- grep("ATC", row.names(obj), value=TRUE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))

obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=50,verbose=F)

obj <- IntegrateLayers(obj, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50, verbose=FALSE)

obj <- RunUMAP(obj, reduction="integrated.rpca", dims=1:50, reduction.name="umap.rpca",verbose=F)
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:50, k.param=40, verbose=F)
obj <- FindClusters(obj, resolution = 3, cluster.name = "clusters", verbose=F)

pdf("/home/moliva/root_sc_paper/data/integrated/c24/obj_clusters.pdf", width=6, height=6)
DimPlot(obj, group.by="clusters", label=TRUE)+NoLegend()
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/c24/obj_umi.pdf", width=6, height=6)
FeaturePlot(obj, "nCount_RNA", min.cutoff="q10", max.cutoff="q90")
dev.off()

#One clusters of low quality "cells"

obj <- subset(obj, idents=42, invert=TRUE)
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

pdf("/home/moliva/root_sc_paper/data/integrated/c24/umap_c24_lineage.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.lineage")
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/c24/umap_c24_identity.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.identity")
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/c24/umap_c24_orig.ident.pdf", width=9, height=6)
DimPlot(obj, group.by="orig.ident")
dev.off()

obj$clusters <- NULL

saveRDS(obj,"/home/moliva/root_sc_paper/data/integrated/c24/integrated_c24_predicted_id.rds")

