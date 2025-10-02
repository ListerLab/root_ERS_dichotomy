library(Seurat)
options(future.globals.maxSize = 2500 * 1024^2)

RL650 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL650_col_cvi_r3/seurat_object/RL650_col_cvi_r3_singlets_genotype.rds")
Idents(RL650) <- "genotype"
RL650 <- subset(RL650, idents="col")
RL650$orig.ident <- "RL650_col_cvi_r3"

RL970 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL970_col_ler_r4/seurat_object/RL970_col_ler_r4_singlets_genotype.rds")
Idents(RL970) <- "genotype"
RL970 <- subset(RL970, idents="col")
RL970$orig.ident <- "RL970_col_ler_r4"

RL971 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL971_col_c24_r5/seurat_object/RL971_col_c24_r5_singlets_genotype.rds")
Idents(RL971) <- "genotype"
RL971 <- subset(RL971, idents="col")
RL971$orig.ident <- "RL971_col_c24_r5"

RL2757 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL2757_col_fullMS_r1/seurat_object/RL2757_col_fullMS_r1_singlets.rds")
RL2757$orig.ident <- "RL2757_col_fullMS_r1"
 
RL2758 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL2758_col_fullMS_r2/seurat_object/RL2758_col_fullMS_r2_singlets.rds")
RL2758$orig.ident <- "RL2758_col_fullMS_r2"
 
RL2759 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL2759_col_halfMS_r1/seurat_object/RL2759_col_halfMS_r1_singlets.rds")
RL2759$orig.ident <- "RL2759_col_halfMS_r1"
 
RL2760 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL2760_col_halfMS_r2/seurat_object/RL2760_col_halfMS_r2_singlets.rds")
RL2760$orig.ident <- "RL2760_col_halfMS_r2"
 
RL4207 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4207_col_fullMS_r3/seurat_object/RL4207_col_fullMS_r3_singlets.rds")
RL4207$orig.ident <- "RL4207_col_fullMS_r3"
 
RL4208 <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4208_col_halfMS_r3/seurat_object/RL4208_col_halfMS_r3_singlets.rds")
RL4208$orig.ident <- "RL4208_col_halfMS_r3"

obj <- merge(RL650, list(RL970, RL971, RL2757, RL2758, RL2759, RL2760, RL4207, RL4208))

obj <- FindVariableFeatures(obj, nfeatures=25000,verbose=F)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj), value=TRUE)
chloro <- grep("ATC", row.names(obj), value=TRUE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))

obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=50,verbose=F)

obj <- IntegrateLayers(obj, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50, reference=1:3,verbose=FALSE)

obj <- RunUMAP(obj, reduction="integrated.rpca", dims=1:50, reduction.name="umap.rpca",verbose=F)
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:50, k.param=50, verbose=F)
obj <- FindClusters(obj, resolution = 0.2, cluster.name = "clusters", verbose=F)

pdf("/home/moliva/root_sc_paper/data/integrated/col/obj_clusters.pdf", width=6, height=6)
DimPlot(obj, group.by="clusters", label=TRUE)+NoLegend()
dev.off()

Idents(obj) <- "clusters"
obj@meta.data[WhichCells(obj, idents=c(0,1,2,3,11)), "lineage"] <- "epidermis_rootcap"
obj@meta.data[WhichCells(obj, idents=7), "lineage"] <- "cortex"
obj@meta.data[WhichCells(obj, idents=c(9,16)), "lineage"] <- "endodermis"
obj@meta.data[WhichCells(obj, idents=c(4,6)), "lineage"] <- "xpp"
obj@meta.data[WhichCells(obj, idents=c(12,13)), "lineage"] <- "ppp"
obj@meta.data[WhichCells(obj, idents=c(8,10)), "lineage"] <- "procambium"
obj@meta.data[WhichCells(obj, idents=c(5,17)), "lineage"] <- "phloem"
obj@meta.data[WhichCells(obj, idents=c(14,15)), "lineage"] <- "xylem"

Idents(obj) <- "lineage"

#There are a few low UMI cells/low quality cells remaining in epidermis and cortex that form diffuse clusters need to be removed. We will process both epidermis and cortex independently to identify and remove them.

##Epidermis - root cap

epirc <- subset(obj, idents="epidermis_rootcap")
epirc <- FindVariableFeatures(epirc, nfeatures=10000,verbose=F)
VariableFeatures(epirc) <- setdiff(VariableFeatures(epirc),c(proto$V1,mito,chloro))
epirc <- ScaleData(epirc)
epirc <- RunPCA(epirc, npcs=50,verbose=F)
epirc <- IntegrateLayers(epirc, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50,verbose=FALSE, reference=1:3)
epirc <- RunUMAP(epirc, reduction="integrated.rpca", dims=1:15, reduction.name="umap.rpca",verbose=F)
epirc <- FindNeighbors(epirc, reduction = "integrated.rpca", dims = 1:15, k.param=30, verbose=F)
epirc <- FindClusters(epirc, resolution = 0.7, cluster.name = "clusters", verbose=F)

pdf("/home/moliva/root_sc_paper/data/integrated/col/epirc_integrated_clusters.pdf", width=6, height=6)
DimPlot(epirc, group.by="clusters", label=TRUE)+NoLegend()
dev.off()

Idents(epirc) <- "clusters"
cells_to_remove <- WhichCells(epirc, idents=15)

##Cortex

cortex <- subset(obj, idents="cortex")
cortex <- FindVariableFeatures(cortex, nfeatures=10000,verbose=F)
VariableFeatures(cortex) <- setdiff(VariableFeatures(cortex),c(proto$V1,mito,chloro))
cortex <- ScaleData(cortex)
cortex <- RunPCA(cortex, npcs=50,verbose=F)
cortex <- IntegrateLayers(cortex, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50,verbose=FALSE, reference=1:3)
cortex <- RunUMAP(cortex, reduction="integrated.rpca", dims=1:20, reduction.name="umap.rpca",verbose=F)
cortex <- FindNeighbors(cortex, reduction = "integrated.rpca", dims = 1:20, k.param=30, verbose=F)
cortex <- FindClusters(cortex, resolution = 0.15, cluster.name = "clusters", verbose=F)


pdf("/home/moliva/root_sc_paper/data/integrated/col/cortex_integrated_clusters.pdf", width=6, height=6)
DimPlot(cortex, group.by="clusters", label=TRUE)+NoLegend()
dev.off()

Idents(cortex) <- "clusters"
cells_to_remove <- c(cells_to_remove,WhichCells(cortex, idents=3))


#Processing integrated object after filtering cells

obj <- subset(obj, cells=cells_to_remove, invert=TRUE)
obj <- FindVariableFeatures(obj, nfeatures=25000,verbose=F)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs=50,verbose=F)
obj <- IntegrateLayers(obj, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50, reference=1:3, verbose=FALSE)
obj <- RunUMAP(obj, reduction="integrated.rpca", dims=1:45, reduction.name="umap.rpca",verbose=F)
obj <- FindNeighbors(obj, reduction = "integrated.rpca", dims = 1:45, k.param=50, verbose=F)
obj <- FindClusters(obj, resolution = 0.6, cluster.name = "clusters", verbose=F)

pdf("/home/moliva/root_sc_paper/data/integrated/col/umap_integrated_col_clusters.pdf", width=6, height=6)
DimPlot(obj, group.by="clusters", label=TRUE)+NoLegend()
dev.off()

pdf("/home/moliva/root_sc_paper/data/integrated/col/umap_integrated_col_orig_ident.pdf", width=8, height=6)
DimPlot(obj, group.by="orig.ident")
dev.off()

Idents(obj) <- "clusters"
obj@meta.data[WhichCells(obj, idents=c(4,5,11,21,23,26)), "lineage"] <- "epidermis"
obj@meta.data[WhichCells(obj, idents=c(0,1,2,9,19,24)), "lineage"] <- "rootcap"
obj@meta.data[WhichCells(obj, idents=c(13,25)), "lineage"] <- "cortex"
obj@meta.data[WhichCells(obj, idents=c(8,20)), "lineage"] <- "endodermis"
obj@meta.data[WhichCells(obj, idents=c(3,6)), "lineage"] <- "xpp"
obj@meta.data[WhichCells(obj, idents=c(12,15)), "lineage"] <- "ppp"
obj@meta.data[WhichCells(obj, idents=c(7,10)), "lineage"] <- "procambium"
obj@meta.data[WhichCells(obj, idents=c(14,17,22)), "lineage"] <- "phloem"
obj@meta.data[WhichCells(obj, idents=c(16,18)), "lineage"] <- "xylem"

Idents(obj) <- "lineage"

pdf("/home/moliva/root_sc_paper/data/integrated/col/umap_integrated_col_lineage.pdf", width=6, height=6)
DimPlot(obj, group.by="lineage", label=TRUE)+NoLegend()
dev.off()

obj <- JoinLayers(obj)
saveRDS(obj,"/home/moliva/root_sc_paper/data/integrated/col/integrated_col.rds")
