library(Seurat)
library(dplyr)
library(ggplot2)

obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")
Idents(obj) <- "identity"

pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/UMI_per_identity_atlas.pdf")
ggplot(obj@meta.data,aes(x=identity,y=nCount_RNA))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

counts <- as.matrix(GetAssayData(obj, assay="RNA", slot="counts"))

down_count_1 <- SampleUMI(counts)
down_count_2 <- SampleUMI(counts)
down_count_3 <- SampleUMI(counts)

down_1 <- CreateSeuratObject(counts= down_count_1, min.cells = 3, min.features = 0)
down_2 <- CreateSeuratObject(counts= down_count_2, min.cells = 3, min.features = 0)
down_3 <- CreateSeuratObject(counts= down_count_3, min.cells = 3, min.features = 0)

proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj), value=TRUE)
chloro <- grep("ATC", row.names(obj), value=TRUE)

process_obj <- function(x){
	x <- NormalizeData(x, verbose=FALSE)
	x <- ScaleData(x, verbose = FALSE)
	x <- FindVariableFeatures(x, nfeatures=25000, verbose=FALSE)
	VariableFeatures(x) <- setdiff(VariableFeatures(x),c(proto$V1,mito,chloro))
	x <- RunPCA(x, npcs=50, verbose=FALSE)
	x <- RunUMAP(x, dims=1:50, verbose=FALSE)
	x$identity <- obj@meta.data[row.names(x@meta.data),"identity"]
	return(x)
}

down.list <- list(down_1,down_2,down_3)
down.list <- lapply(down.list,process_obj)

down_1 <- down.list[[1]]
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/umap_down_1.pdf",width=9, height=6)
DimPlot(down_1,group.by="identity")
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/UMI_per_identity_down_1.pdf")
ggplot(down_1@meta.data,aes(x=identity,y=nCount_RNA))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
early_genes <- readRDS("/home/moliva/root_sc_paper/analyses/trajectories/early_genes.rds")
down_1 <- AddModuleScore(down_1, features=list(early_genes), name="score")
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/early_genes_score_down_1.pdf")
FeaturePlot(down_1,"score1", min.cutoff=0, max.cutoff="q50",order=TRUE)
dev.off()

down_2 <- down.list[[2]]
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/umap_down_2.pdf",width=9, height=6)
DimPlot(down_2,group.by="identity")
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/UMI_per_identity_down_2.pdf")
ggplot(down_2@meta.data,aes(x=identity,y=nCount_RNA))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
down_2 <- AddModuleScore(down_2, features=list(early_genes), name="score")
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/early_genes_score_down_2.pdf")
FeaturePlot(down_2,"score1", min.cutoff=0, max.cutoff="q50",order=TRUE)
dev.off()

down_3 <- down.list[[3]]
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/umap_down_3.pdf",width=9, height=6)
DimPlot(down_3,group.by="identity")
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/UMI_per_identity_down_3.pdf")
ggplot(down_3@meta.data,aes(x=identity,y=nCount_RNA))+geom_boxplot()+theme_classic()+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
down_3 <- AddModuleScore(down_3, features=list(early_genes), name="score")
pdf("/home/moliva/root_sc_paper/analyses/QC/downsampling/early_genes_score_down_3.pdf")
FeaturePlot(down_3,"score1", min.cutoff=0, max.cutoff="q50",order=TRUE)
dev.off()

saveRDS(down_1,"/home/moliva/root_sc_paper/analyses/QC/downsampling/down_1.rds")
saveRDS(down_2,"/home/moliva/root_sc_paper/analyses/QC/downsampling/down_2.rds")
saveRDS(down_3,"/home/moliva/root_sc_paper/analyses/QC/downsampling/down_3.rds")


