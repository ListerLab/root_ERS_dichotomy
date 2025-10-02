library(Seurat)
library(ggplot2)

obj <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4209_col_through_agar/seurat_object/RL4209_col_through_agar_singlets.rds")

ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

obj <- NormalizeData(obj)
obj <- ScaleData(obj)
obj <- FindVariableFeatures(obj, nfeatures=25000,verbose=F)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj), value=TRUE)
chloro <- grep("ATC", row.names(obj), value=TRUE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))
obj <- RunPCA(obj, npcs=50,verbose=F)
obj <- RunUMAP(obj, dims=1:50, verbose=F)

anchors <- FindTransferAnchors(reference=ref, query=obj, dims=1:50, reference.reduction="integrated.rpca")
predictions_lineage <- TransferData(anchorset=anchors, refdata=ref$lineage, dims=1:50)
obj <-AddMetaData(obj, metadata=predictions_lineage[,1,drop=FALSE], col.name="predicted.lineage")
predictions_identity <- TransferData(anchorset=anchors, refdata=ref$identity, dims=1:50)
obj <-AddMetaData(obj, metadata=predictions_identity[,1,drop=FALSE], col.name="predicted.identity")

obj <- RunUMAP(obj, dims=1:50, verbose=F)

pdf("/home/moliva/root_sc_paper/analyses/QC/through_agar/through_agar_umap.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.identity")
dev.off()

ref <- RunUMAP(ref, reduction="integrated.rpca", dims=1:45, reduction.name="umap.rpca",verbose=F, return.model=TRUE)
obj <- MapQuery(anchorset=anchors, reference=ref, query=obj, refdata=list(identity="identity"), reference.reduction="integrated.rpca", reduction.model="umap.rpca")


p1 <- DimPlot(ref, reduction="umap.rpca",group.by="identity")+NoLegend()+ggtitle("Reference Atlas")
p2 <- DimPlot(obj, reduction="ref.umap", group.by="predicted.identity")+ ggtitle("through agar projection")
pdf("/home/moliva/root_sc_paper/analyses/QC/through_agar/through_agar_umap_ref_projection.pdf", width=16, height=6)
p1+p2
dev.off()
