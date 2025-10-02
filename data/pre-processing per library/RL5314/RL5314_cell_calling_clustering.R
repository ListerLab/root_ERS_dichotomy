library(Seurat)
library(dplyr)
library(pheatmap)

data <- Read10X("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cellranger_RL5314_drought_20percentH2O_r1/outs/raw_feature_bc_matrix/")
cr_calling <- read.csv("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cellranger_RL5314_drought_20percentH2O_r1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sep='\t', header = FALSE)$V1

bc_order <- colSums(data)
bc_order <- names(bc_order[order(bc_order, decreasing = TRUE)])

background <- bc_order[!(bc_order %in% cr_calling)]
background <- background[1:length(cr_calling)]

obj <- CreateSeuratObject(counts= data[, cr_calling], min.cells = 3, min.features = 0)
obj <- NormalizeData(obj, verbose=FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures=15000, verbose=FALSE)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj), value=TRUE)
chloro <- grep("ATC", row.names(obj), value=TRUE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))
obj <- RunPCA(obj, npcs=50, verbose=FALSE)
obj <- RunUMAP(obj, dims=1:40, verbose=FALSE)
obj <- FindNeighbors(obj, dims=1:40, verbose=FALSE)
obj <- FindClusters(obj, resolution=0.8, verbose=FALSE)

dir.create("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/")
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_cr_calling_clusters.pdf", width =6, height=6)
DimPlot(obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

genes <- c("AT1G05010","AT1G44760","AT1G26820","AT2G04025","AT1G66470","AT2G39530","AT1G12090","AT2G14900","AT2G18800","AT1G27030","AT4G23410","AT1G68810","AT3G15353","AT4G14020","AT4G14650","AT2G28660")
avg <- AverageExpression(obj)
avg <- avg$RNA
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/pheatmap_markers_clusters.pdf", width=10, height=5)
pheatmap(avg[genes,], scale = "row", cluster_rows = FALSE)
dev.off()

root_cap <- WhichCells(obj, idents=c(3,4,8,10))
epidermis <- WhichCells(obj, idents=c(0,1,2,5,6,7,11,16,17,23))
cortex <- WhichCells(obj, idents=c(9,15,24))
endodermis <- WhichCells(obj, idents=c(12,18,22))
pericycle <- WhichCells(obj, idents=c(13,14))
phloem <- WhichCells(obj, idents=20)
xylem <- WhichCells(obj, idents=21)
procambium <- WhichCells(obj, idents=19)

process_obj <- function(obj){
	obj <- ScaleData(obj, verbose = FALSE)
	obj <- FindVariableFeatures(obj, nfeatures=5000, verbose=FALSE)
	obj <- RunPCA(obj, npcs=50, verbose=FALSE)
	VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))
	obj <- RunPCA(obj, npcs=50, verbose=FALSE)
	obj <- RunUMAP(obj, dims=1:20, verbose=FALSE)
	obj <- FindNeighbors(obj, dims=1:20, verbose=FALSE)
	return(obj)
}

#Root cap

root_cap_obj <- CreateSeuratObject(counts= data[, c(root_cap,background)], min.cells = 3, min.features = 0)
root_cap_obj <- NormalizeData(root_cap_obj, verbose=FALSE)
root_cap_obj <- process_obj(root_cap_obj)

root_cap_obj@meta.data[background,"type"] <- "background"
root_cap_obj@meta.data[root_cap,"type"] <- "root_cap"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_root_cap_1.pdf", width =6, height=6)
DimPlot(root_cap_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

root_cap_obj <- FindClusters(root_cap_obj, resolution=0.1, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_root_cap_clusters_1.pdf", width =6, height=6)
DimPlot(root_cap_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

root_cap_obj <- subset(root_cap_obj, idents=2)
root_cap_obj <- process_obj(root_cap_obj)

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_root_cap_2.pdf", width =6, height=6)
DimPlot(root_cap_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

root_cap_obj <- FindClusters(root_cap_obj, resolution=0.4, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_root_cap_clusters_2.pdf", width =6, height=6)
DimPlot(root_cap_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

root_cap_obj@meta.data[WhichCells(root_cap_obj, idents=5),"type"] <- "background"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_root_cap_3.pdf", width =6, height=6)
DimPlot(root_cap_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

Idents(root_cap_obj) <- "type"
root_cap <- WhichCells(root_cap_obj, idents="root_cap")


#Epidermis

epidermis_obj <- CreateSeuratObject(counts= data[, c(epidermis,background)], min.cells = 3, min.features = 0)
epidermis_obj <- NormalizeData(epidermis_obj, verbose=FALSE)
epidermis_obj <- process_obj(epidermis_obj)


epidermis_obj@meta.data[background,"type"] <- "background"
epidermis_obj@meta.data[epidermis,"type"] <- "epidermis"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_epidermis_1.pdf", width =6, height=6)
DimPlot(epidermis_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

epidermis_obj <- FindClusters(epidermis_obj, resolution=0.2, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_epidermis_clusters_1.pdf", width =6, height=6)
DimPlot(epidermis_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

epidermis_obj <- subset(epidermis_obj, idents=c(0,2,3,4,5,6,7))
epidermis_obj <- process_obj(epidermis_obj)

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_epidermis_2.pdf", width =6, height=6)
DimPlot(epidermis_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

epidermis_obj <- FindClusters(epidermis_obj, resolution=1, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_epidermis_clusters_2.pdf", width =6, height=6)
DimPlot(epidermis_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

epidermis_obj@meta.data[WhichCells(epidermis_obj, idents=c(0,1,2,4,5,7,8,10,15,17,19)),"type"] <- "background"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_epidermis_3.pdf", width =6, height=6)
DimPlot(epidermis_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

Idents(epidermis_obj) <- "type"
epidermis <- WhichCells(epidermis_obj, idents="epidermis")


#Cortex

cortex_obj <- CreateSeuratObject(counts= data[, c(cortex,background)], min.cells = 3, min.features = 0)
cortex_obj <- NormalizeData(cortex_obj, verbose=FALSE)
cortex_obj <- process_obj(cortex_obj)


cortex_obj@meta.data[background,"type"] <- "background"
cortex_obj@meta.data[cortex,"type"] <- "cortex"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_cortex_1.pdf", width =6, height=6)
DimPlot(cortex_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

cortex_obj <- FindClusters(cortex_obj, resolution=0.05, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_cortex_clusters_1.pdf", width =6, height=6)
DimPlot(cortex_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

cortex_obj <- subset(cortex_obj, idents=c(1,2))
cortex_obj <- process_obj(cortex_obj)

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_cortex_2.pdf", width =6, height=6)
DimPlot(cortex_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

cortex_obj <- FindClusters(cortex_obj, resolution=0.05, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_cortex_clusters_2.pdf", width =6, height=6)
DimPlot(cortex_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

cortex_obj@meta.data[WhichCells(cortex_obj, idents=1),"type"] <- "background"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_cortex_3.pdf", width =6, height=6)
DimPlot(cortex_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

Idents(cortex_obj) <- "type"
cortex <- WhichCells(cortex_obj, idents="cortex")


#Endodermis

endodermis_obj <- CreateSeuratObject(counts= data[, c(endodermis,background)], min.cells = 3, min.features = 0)
endodermis_obj <- NormalizeData(endodermis_obj, verbose=FALSE)
endodermis_obj <- process_obj(endodermis_obj)

endodermis_obj@meta.data[background,"type"] <- "background"
endodermis_obj@meta.data[endodermis,"type"] <- "endodermis"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_endodermis_1.pdf", width =6, height=6)
DimPlot(endodermis_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

endodermis_obj <- FindClusters(endodermis_obj, resolution=0.05, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_endodermis_clusters_1.pdf", width =6, height=6)
DimPlot(endodermis_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

endodermis_obj <- subset(endodermis_obj, idents=1)
endodermis_obj <- process_obj(endodermis_obj)

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_endodermis_2.pdf", width =6, height=6)
DimPlot(endodermis_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

endodermis_obj <- FindClusters(endodermis_obj, resolution=0.3, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_endodermis_clusters_2.pdf", width =6, height=6)
DimPlot(endodermis_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

endodermis_obj@meta.data[WhichCells(endodermis_obj, idents=4),"type"] <- "background"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_endodermis_3.pdf", width =6, height=6)
DimPlot(endodermis_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

Idents(endodermis_obj) <- "type"
endodermis <- WhichCells(endodermis_obj, idents="endodermis")


#Pericycle

pericycle_obj <- CreateSeuratObject(counts= data[, c(pericycle,background)], min.cells = 3, min.features = 0)
pericycle_obj <- NormalizeData(pericycle_obj, verbose=FALSE)
pericycle_obj <- process_obj(pericycle_obj)

pericycle_obj@meta.data[background,"type"] <- "background"
pericycle_obj@meta.data[pericycle,"type"] <- "pericycle"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_pericycle_1.pdf", width =6, height=6)
DimPlot(pericycle_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

pericycle_obj <- FindClusters(pericycle_obj, resolution=0.03, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_pericycle_clusters_1.pdf", width =6, height=6)
DimPlot(pericycle_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

pericycle_obj <- subset(pericycle_obj, idents=1)
pericycle_obj <- process_obj(pericycle_obj)

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_pericycle_2.pdf", width =6, height=6)
DimPlot(pericycle_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

pericycle_obj <- FindClusters(pericycle_obj, resolution=0.5, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_pericycle_clusters_2.pdf", width =6, height=6)
DimPlot(pericycle_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

pericycle_obj@meta.data[WhichCells(pericycle_obj, idents=1),"type"] <- "background"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_pericycle_3.pdf", width =6, height=6)
DimPlot(pericycle_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

Idents(pericycle_obj) <- "type"
pericycle <- WhichCells(pericycle_obj, idents="pericycle")


#Phloem

phloem_obj <- CreateSeuratObject(counts= data[, c(phloem,background)], min.cells = 3, min.features = 0)
phloem_obj <- NormalizeData(phloem_obj, verbose=FALSE)
phloem_obj <- process_obj(phloem_obj)

phloem_obj@meta.data[background,"type"] <- "background"
phloem_obj@meta.data[phloem,"type"] <- "phloem"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_phloem_1.pdf", width =6, height=6)
DimPlot(phloem_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

phloem_obj <- FindClusters(phloem_obj, resolution=0.02, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_phloem_clusters_1.pdf", width =6, height=6)
DimPlot(phloem_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

phloem_obj <- subset(phloem_obj, idents=1)
phloem_obj <- process_obj(phloem_obj)

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_phloem_2.pdf", width =6, height=6)
DimPlot(phloem_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

Idents(phloem_obj) <- "type"
phloem <- WhichCells(phloem_obj, idents="phloem")


#Xylem

xylem_obj <- CreateSeuratObject(counts= data[, c(xylem,background)], min.cells = 3, min.features = 0)
xylem_obj <- NormalizeData(xylem_obj, verbose=FALSE)
xylem_obj <- process_obj(xylem_obj)

xylem_obj@meta.data[background,"type"] <- "background"
xylem_obj@meta.data[xylem,"type"] <- "xylem"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_xylem_1.pdf", width =6, height=6)
DimPlot(xylem_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

xylem_obj <- FindClusters(xylem_obj, resolution=0.02, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_xylem_clusters_1.pdf", width =6, height=6)
DimPlot(xylem_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

xylem_obj <- subset(xylem_obj, idents=1)
xylem_obj <- process_obj(xylem_obj)

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_xylem_2.pdf", width =6, height=6)
DimPlot(xylem_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

Idents(xylem_obj) <- "type"
xylem <- WhichCells(xylem_obj, idents="xylem")


#Procambium

procambium_obj <- CreateSeuratObject(counts= data[, c(procambium,background)], min.cells = 3, min.features = 0)
procambium_obj <- NormalizeData(procambium_obj, verbose=FALSE)
procambium_obj <- process_obj(procambium_obj)

procambium_obj@meta.data[background,"type"] <- "background"
procambium_obj@meta.data[procambium,"type"] <- "procambium"

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_procambium_1.pdf", width =6, height=6)
DimPlot(procambium_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

procambium_obj <- FindClusters(procambium_obj, resolution=0.05, verbose=FALSE)
pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_procambium_clusters_1.pdf", width =6, height=6)
DimPlot(procambium_obj, group.by="seurat_clusters", label=TRUE)+NoLegend()
dev.off()

procambium_obj <- subset(procambium_obj, idents=c(1,2))
procambium_obj <- process_obj(procambium_obj)

pdf("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/cell_calling/umap_procambium_2.pdf", width =6, height=6)
DimPlot(procambium_obj, group.by="type", label=TRUE)+NoLegend()
dev.off()

Idents(procambium_obj) <- "type"
procambium <- WhichCells(procambium_obj, idents="procambium")


#Final object
cells <- c(root_cap,epidermis,cortex,endodermis,pericycle,phloem,xylem,procambium)

obj <- CreateSeuratObject(counts= data[, cells], min.cells = 3, min.features = 0)
obj <- NormalizeData(obj, verbose=FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures=15000, verbose=FALSE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))
obj <- RunPCA(obj, npcs=50, verbose=FALSE)
obj <- RunUMAP(obj, dims=1:40, verbose=FALSE)

dir.create("/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/seurat_object/")
saveRDS(obj,"/home/moliva/root_sc_paper/data/raw_data_per_library/RL5314_drought_20percentH2O_r1/seurat_object/RL5314_drought_20percentH2O_r1.rds")