library(Seurat)
library(ggplot2)

obj <- readRDS("/group/llshared/rambo/processed_snapatac2_signac/processedATAC.afterIntegration.RemoveFirstComponent.rds")
ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

anchors <- FindTransferAnchors(reference=ref, query=obj, reduction='cca', dims=1:30, query.assay='ACTIVITY')
predictions.identity <- TransferData(anchorset=anchors, refdata=ref$identity, weight.reduction=obj[["lsi"]], dims=2:30)
predictions.lineage <- TransferData(anchorset=anchors, refdata=ref$lineage, weight.reduction=obj[["lsi"]], dims=2:30)

genes.use <- VariableFeatures(ref)
refdata <- GetAssayData(ref, assay="RNA",slot="data")[genes.use,]
imputation <- TransferData(anchorset=anchors, refdata=refdata, weight.reduction=obj[["lsi"]], dims=2:30)
obj[["RNA"]] <- imputation

obj$tech <- "ATAC"
ref$tech <- "RNA"
obj$lineage <- predictions.lineage$predict
obj$lineage <- predictions.lineage$predicted.id
obj$lineage_max <- predictions.lineage$prediction.score.max
obj$identity <- predictions.identity$predicted.id
obj$identity_max <- predictions.identity$prediction.score.max

coembed <- merge(ref,obj)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, npcs=50, verbose=FALSE)
coembed <- RunUMAP(coembed, dims = 1:50)

pdf("coembedding.pdf", width=15, height=6)
DimPlot(coembed, group.by = c("tech","identity"))
dev.off()














anchors <- FindTransferAnchors(reference=ref, query=obj, dims=1:50, query.assay="ACTIVITY_2kb",reduction='cca')
obj<- MapQuery(anchorset=anchors, reference=ref, query=obj, refdata=list(identity="identity"), reduction.model="umap.rpca")





pdf("/home/moliva/root_sc_paper/analyses/snATAC/snATAC_projection.pdf", width=16, height=6)
p1+p2
dev.off()