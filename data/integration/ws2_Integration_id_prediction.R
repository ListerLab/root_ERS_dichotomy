library(Seurat)
options(future.globals.maxSize = 2500 * 1024^2)

RL1399 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL1399_ws2_N868100/seurat_object/RL1399_ws2_N868100_singlets_genotype.rds")
Idents(RL1399) <- "genotype"
RL1399 <- subset(RL1399, idents="ws")
RL1399$orig.ident <- "RL1399_ws2"

RL1463 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL1463_ws2_N654328/seurat_object/RL1463_ws2_N654328_singlets_genotype.rds")
Idents(RL1463) <- "genotype"
RL1463 <- subset(RL1463, idents="ws")
RL1463$orig.ident <- "RL1463_ws2"

RL1464 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL1464_ws2_N655192/seurat_object/RL1464_ws2_N655192_singlets_genotype.rds")
Idents(RL1464) <- "genotype"
RL1464 <- subset(RL1464, idents="ws")
RL1464$orig.ident <- "RL1464_ws2"

obj <- merge(RL1399, list(RL1463, RL1464))

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
obj <- FindClusters(obj, resolution = 8, cluster.name = "clusters", verbose=F)

pdf("/home/moliva/root_sc_paper/data/integrated/ws2/obj_clusters.pdf", width=6, height=6)
DimPlot(obj, group.by="clusters", label=TRUE)+NoLegend()
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/ws2/obj_umi.pdf", width=6, height=6)
FeaturePlot(obj, "nCount_RNA", min.cutoff="q10", max.cutoff="q90")
dev.off()

#Two clusters of low-UMI "cells" to be removed because likely low quality "cells"

obj <- subset(obj, idents=c(17,74), invert=TRUE)
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

pdf("/home/moliva/root_sc_paper/data/integrated/ws2/umap_ws2_lineage.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.lineage")
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/ws2/umap_ws2_identity.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.identity")
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/ws2/umap_ws2_orig.ident.pdf", width=9, height=6)
DimPlot(obj, group.by="orig.ident")
dev.off()

obj$clusters <- NULL

saveRDS(obj,"/home/moliva/root_sc_paper/data/integrated/ws2/integrated_ws2_predicted_id.rds")


 
 
 
 
 