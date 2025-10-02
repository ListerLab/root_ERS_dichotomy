library(Seurat)
options(future.globals.maxSize = 2500 * 1024^2)

RL650 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL650_col_cvi_r3/seurat_object/RL650_col_cvi_r3_singlets_genotype.rds")
Idents(RL650) <- "genotype"
RL650 <- subset(RL650,idents="cvi")
RL650$orig.ident <- "RL650_col_cvi"

RL860 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL860_scr3_cvi/seurat_object/RL860_scr3_cvi_singlets_genotype.rds")
Idents(RL860) <- "genotype"
RL860 <- subset(RL860, idents="cvi")
RL860$orig.ident <- "RL860_scr3_cvi"

RL4460 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4460_cvi_CDFs/seurat_object/RL4460_cvi_CDFs_singlets_genotype.rds")
Idents(RL4460) <- "genotype"
RL4460 <- subset(RL4460, idents="cvi")
RL4460$orig.ident <- "RL4460_cvi_CDFs"

obj <- merge(RL650, list(RL860, RL4460))

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
obj <- FindClusters(obj, resolution = 0.4, cluster.name = "clusters", verbose=F)

pdf("/home/moliva/root_sc_paper/data/integrated/cvi/obj_clusters.pdf", width=6, height=6)
DimPlot(obj, group.by="clusters", label=TRUE)+NoLegend()
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/cvi/obj_umi.pdf", width=6, height=6)
FeaturePlot(obj, "nCount_RNA", min.cutoff="q10", max.cutoff="q90")
dev.off()

#One diffuse clusters of low-UMI "cells" to be removed because likely low quality "cells"

obj <- subset(obj, idents=20, invert=TRUE)
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

pdf("/home/moliva/root_sc_paper/data/integrated/cvi/umap_cvi_lineage.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.lineage")
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/cvi/umap_cvi_identity.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.identity")
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/cvi/umap_cvi_orig.ident.pdf", width=9, height=6)
DimPlot(obj, group.by="orig.ident")
dev.off()

obj$clusters <- NULL

saveRDS(obj,"/home/moliva/root_sc_paper/data/integrated/cvi/integrated_cvi_predicted_id.rds")



