library(Seurat)


obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col.rds")
Idents(obj) <- "lineage"


#Epidermis

epi <- subset(obj, idents="epidermis")
epi <- FindNeighbors(epi, reduction = "integrated.rpca", dims = 1:45, k.param=200, verbose=F)
epi <- FindClusters(epi, resolution = 0.82, cluster.name = "trajectory_clusters", verbose=F)
pdf("/home/moliva/root_sc_paper/analyses/trajectories/epi_clusters.pdf", width=6, height=6)
DimPlot(epi, group.by="trajectory_clusters", label=TRUE)+NoLegend()
dev.off()

emb <- as.data.frame(epi@reductions$umap.rpca@cell.embeddings[WhichCells(epi, ident=5),])
root_cell <- row.names(emb[emb$umaprpca_2 == min(emb$umaprpca_2),])

epi_paga <- epi
epi_paga[["RNA"]] <- as(epi_paga[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)

epi_paga <- RunPAGA(epi_paga, group_by="trajectory_clusters", linear_reduction="integrated.rpca", nonlinear_reduction="umap.rpca", n_pcs=45, n_neighbors=200, infer_pseudotime=TRUE, root_cell=root_cell)

pdf("/home/moliva/root_sc_paper/analyses/trajectories/epi_paga.pdf", width=6, height=6)
PAGAPlot(epi_paga, reduction="umap.rpca", edge_threshold=0.198)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/epi_pseudotime.pdf", width=6, height=6)
FeatureDimPlot(epi_paga, features = "dpt_pseudotime", reduction="umap.rpca")
dev.off()

detach('package:SCP', unload=TRUE)
detach('package:Seurat', unload=TRUE)
library(Seurat)

obj@meta.data[row.names(epi@meta.data),"trajectory_clusters"] <- paste0("epi_",epi$trajectory_clusters)
obj@meta.data[row.names(epi_paga@meta.data),"pseudotime"] <- epi_paga$dpt_pseudotime
obj@meta.data[WhichCells(epi, idents=5), "identity"] <- "epidermis early"
obj@meta.data[WhichCells(epi, idents=c(2,3,7,10)), "identity"] <- "non-hair cells"
obj@meta.data[WhichCells(epi, idents=c(1,9)), "identity"] <- "hair cells early"
obj@meta.data[WhichCells(epi, idents=c(0,4)), "identity"] <- "hair cells ERS+"
obj@meta.data[WhichCells(epi, idents=c(6,8)), "identity"] <- "hair cells ERS-"
rm(epi_paga)

early_genes <- FindMarkers(epi, ident.1=5, only.pos=TRUE, test.use="t")
early_genes <- row.names(subset(early_genes, p_val_adj < 0.01))
saveRDS(early_genes,"/home/moliva/root_sc_paper/analyses/trajectories/early_genes.rds")

obj <- AddModuleScore(obj, features=list(early_genes))
colnames(obj@meta.data)[11] <- "early_epi_genes_score"
rm(epi)

#Root cap

rc <- subset(obj, idents="rootcap")


rc <- FindNeighbors(rc, reduction = "integrated.rpca", dims = 1:45, k.param=90, verbose=F)
rc <- FindClusters(rc, resolution = 0.28, cluster.name = "trajectory_clusters", verbose=F)
pdf("/home/moliva/root_sc_paper/analyses/trajectories/rc_clusters.pdf", width=6, height=6)
DimPlot(rc, group.by="trajectory_clusters", label=TRUE)+NoLegend()
dev.off()

emb <- as.data.frame(rc@reductions$umap.rpca@cell.embeddings[WhichCells(rc, ident=5),])
emb <- subset(emb, umaprpca_2 <5)
root_cell <- row.names(emb[emb$umaprpca_2 == max(emb$umaprpca_2),])

rc_paga <- rc
rc_paga[["RNA"]] <- as(rc[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)

rc_paga <- RunPAGA(rc_paga, group_by="trajectory_clusters", linear_reduction="integrated.rpca", nonlinear_reduction="umap.rpca", n_pcs=45, n_neighbors=90, infer_pseudotime=TRUE, root_cell=root_cell)

pdf("/home/moliva/root_sc_paper/analyses/trajectories/rc_paga.pdf", width=6, height=6)
PAGAPlot(rc_paga, reduction="umap.rpca", edge_threshold=0.06)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/rc_pseudotime.pdf", width=6, height=6)
FeatureDimPlot(rc_paga, features = "dpt_pseudotime", reduction="umap.rpca")
dev.off()

detach('package:SCP', unload=TRUE)
detach('package:Seurat', unload=TRUE)
library(Seurat)

obj@meta.data[row.names(rc@meta.data),"trajectory_clusters"] <- paste0("rc_",rc$trajectory_clusters)
obj@meta.data[row.names(rc_paga@meta.data),"pseudotime"] <- rc_paga$dpt_pseudotime
obj@meta.data[WhichCells(rc, idents=c(5,0)), "identity"] <- "root cap early"
obj@meta.data[WhichCells(rc, idents=6), "identity"] <- "columella"
obj@meta.data[WhichCells(rc, idents=4), "identity"] <- "border cells"
obj@meta.data[WhichCells(rc, idents=1), "identity"] <- "root cap tip"
obj@meta.data[WhichCells(rc, idents=c(2,3)), "identity"] <- "root cap lateral"
rm(rc_paga)
rm(rc)


#Cortex

cortex <- subset(obj, idents="cortex")

cortex <- FindNeighbors(cortex, reduction = "integrated.rpca", dims = 1:45, k.param=50, verbose=F)
cortex <- FindClusters(cortex, resolution = 0.7, cluster.name = "trajectory_clusters", verbose=F)
pdf("/home/moliva/root_sc_paper/analyses/trajectories/cortex_clusters.pdf", width=6, height=6)
DimPlot(cortex, group.by="trajectory_clusters", label=TRUE)+NoLegend()
dev.off()


root_cell <- row.names(cortex@meta.data[cortex$early_epi_genes_score==max(cortex$early_epi_genes_score),])

cortex_paga <- cortex
cortex_paga[["RNA"]] <- as(cortex[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)

cortex_paga <- RunPAGA(cortex_paga, group_by="trajectory_clusters", linear_reduction="integrated.rpca", nonlinear_reduction="umap.rpca", n_pcs=45, n_neighbors=50, infer_pseudotime=TRUE, root_cell=root_cell)

pdf("/home/moliva/root_sc_paper/analyses/trajectories/cortex_paga.pdf", width=6, height=6)
PAGAPlot(cortex_paga, reduction="umap.rpca", edge_threshold=0.06)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/cortex_pseudotime.pdf", width=6, height=6)
FeatureDimPlot(cortex_paga, features = "dpt_pseudotime", reduction="umap.rpca")
dev.off()

detach('package:SCP', unload=TRUE)
detach('package:Seurat', unload=TRUE)
library(Seurat)

obj@meta.data[row.names(cortex@meta.data),"trajectory_clusters"] <- paste0("cortex_",cortex$trajectory_clusters)
obj@meta.data[row.names(cortex_paga@meta.data),"pseudotime"] <- cortex_paga$dpt_pseudotime
obj@meta.data[WhichCells(cortex, idents=c(0,2,5)), "identity"] <- "cortex ERS+"
obj@meta.data[WhichCells(cortex, idents=c(1,3,4)), "identity"] <- "cortex ERS-"
rm(cortex_paga)
rm(cortex)

#Endodermis

endo <- subset(obj, idents="endodermis")

endo <- FindNeighbors(endo, reduction = "integrated.rpca", dims = 1:45, k.param=70, verbose=F)
endo <- FindClusters(endo, resolution = 0.3, cluster.name = "trajectory_clusters", verbose=F)
pdf("/home/moliva/root_sc_paper/analyses/trajectories/endo_clusters.pdf", width=6, height=6)
DimPlot(endo, group.by="trajectory_clusters", label=TRUE)+NoLegend()
dev.off()


emb <- as.data.frame(endo@reductions$umap.rpca@cell.embeddings)
cells <- row.names(subset(emb, umaprpca_2 >5))
root_cell <- endo@meta.data[cells,]
root_cell <- row.names(root_cell[root_cell$early_epi_genes_score==max(root_cell$early_epi_genes_score),])

endo_paga <- endo
endo_paga[["RNA"]] <- as(endo[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)

endo_paga <- RunPAGA(endo_paga, group_by="trajectory_clusters", linear_reduction="integrated.rpca", nonlinear_reduction="umap.rpca", n_pcs=45, n_neighbors=70, infer_pseudotime=TRUE, root_cell=root_cell)

pdf("/home/moliva/root_sc_paper/analyses/trajectories/endo_paga.pdf", width=6, height=6)
PAGAPlot(endo_paga, reduction="umap.rpca", edge_threshold=0.05)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/endo_pseudotime.pdf", width=6, height=6)
FeatureDimPlot(endo_paga, features = "dpt_pseudotime", reduction="umap.rpca")
dev.off()

detach('package:SCP', unload=TRUE)
detach('package:Seurat', unload=TRUE)
library(Seurat)

obj@meta.data[row.names(endo@meta.data),"trajectory_clusters"] <- paste0("endo_",endo$trajectory_clusters)
obj@meta.data[row.names(endo_paga@meta.data),"pseudotime"] <- endo_paga$dpt_pseudotime
obj@meta.data[WhichCells(endo, idents=c(0,1,4)), "identity"] <- "endodermis ERS+"
obj@meta.data[WhichCells(endo, idents=c(2,3)), "identity"] <- "endodermis ERS-"
rm(endo_paga)
rm(endo)


#XPP

xpp <- subset(obj, idents="xpp")

xpp <- FindNeighbors(xpp, reduction = "integrated.rpca", dims = 1:45, k.param=70, verbose=F)
xpp <- FindClusters(xpp, resolution = 0.3, cluster.name = "trajectory_clusters", verbose=F)
pdf("/home/moliva/root_sc_paper/analyses/trajectories/xpp_clusters.pdf", width=6, height=6)
DimPlot(xpp, group.by="trajectory_clusters", label=TRUE)+NoLegend()
dev.off()


root_cell <- row.names(xpp@meta.data[xpp$early_epi_genes_score==max(xpp$early_epi_genes_score),])

xpp_paga <- xpp
xpp_paga[["RNA"]] <- as(xpp[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)

xpp_paga <- RunPAGA(xpp_paga, group_by="trajectory_clusters", linear_reduction="integrated.rpca", nonlinear_reduction="umap.rpca", n_pcs=45, n_neighbors=70, infer_pseudotime=TRUE, root_cell=root_cell)

pdf("/home/moliva/root_sc_paper/analyses/trajectories/xpp_paga.pdf", width=6, height=6)
PAGAPlot(xpp_paga, reduction="umap.rpca", edge_threshold=0.05)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/xpp_pseudotime.pdf", width=6, height=6)
FeatureDimPlot(xpp_paga, features = "dpt_pseudotime", reduction="umap.rpca")
dev.off()

detach('package:SCP', unload=TRUE)
detach('package:Seurat', unload=TRUE)
library(Seurat)

obj@meta.data[row.names(xpp@meta.data),"trajectory_clusters"] <- paste0("xpp_",xpp$trajectory_clusters)
obj@meta.data[row.names(xpp_paga@meta.data),"pseudotime"] <- xpp_paga$dpt_pseudotime
obj@meta.data[WhichCells(xpp, idents=c(0,2)), "identity"] <- "xpp ERS+"
obj@meta.data[WhichCells(xpp, idents=c(1,3)), "identity"] <- "xpp ERS-"
rm(xpp_paga)
rm(xpp)


#PPP

ppp <- subset(obj, idents="ppp")

ppp <- FindNeighbors(ppp, reduction = "integrated.rpca", dims = 1:45, k.param=70, verbose=F)
ppp <- FindClusters(ppp, resolution = 0.3, cluster.name = "trajectory_clusters", verbose=F)
pdf("/home/moliva/root_sc_paper/analyses/trajectories/ppp_clusters.pdf", width=6, height=6)
DimPlot(ppp, group.by="trajectory_clusters", label=TRUE)+NoLegend()
dev.off()


root_cell <- row.names(ppp@meta.data[ppp$early_epi_genes_score==max(ppp$early_epi_genes_score),])

ppp_paga <- ppp
ppp_paga[["RNA"]] <- as(ppp[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)

ppp_paga <- RunPAGA(ppp_paga, group_by="trajectory_clusters", linear_reduction="integrated.rpca", nonlinear_reduction="umap.rpca", n_pcs=45, n_neighbors=70, infer_pseudotime=TRUE, root_cell=root_cell)

pdf("/home/moliva/root_sc_paper/analyses/trajectories/ppp_paga.pdf", width=6, height=6)
PAGAPlot(ppp_paga, reduction="umap.rpca", edge_threshold=0.05)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/ppp_pseudotime.pdf", width=6, height=6)
FeatureDimPlot(ppp_paga, features = "dpt_pseudotime", reduction="umap.rpca")
dev.off()

detach('package:SCP', unload=TRUE)
detach('package:Seurat', unload=TRUE)
library(Seurat)

obj@meta.data[row.names(ppp@meta.data),"trajectory_clusters"] <- paste0("ppp_",ppp$trajectory_clusters)
obj@meta.data[row.names(ppp_paga@meta.data),"pseudotime"] <- ppp_paga$dpt_pseudotime
obj@meta.data[WhichCells(ppp, idents=c(0,3)), "identity"] <- "ppp ERS+"
obj@meta.data[WhichCells(ppp, idents=c(1,2)), "identity"] <- "ppp ERS-"
rm(ppp_paga)
rm(ppp)


#procambium

procambium <- subset(obj, idents="procambium")

procambium <- FindNeighbors(procambium, reduction = "integrated.rpca", dims = 1:45, k.param=70, verbose=F)
procambium <- FindClusters(procambium, resolution = 0.3, cluster.name = "trajectory_clusters", verbose=F)
pdf("/home/moliva/root_sc_paper/analyses/trajectories/procambium_clusters.pdf", width=6, height=6)
DimPlot(procambium, group.by="trajectory_clusters", label=TRUE)+NoLegend()
dev.off()


root_cell <- row.names(procambium@meta.data[procambium$early_epi_genes_score==max(procambium$early_epi_genes_score),])

procambium_paga <- procambium
procambium_paga[["RNA"]] <- as(procambium[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)

procambium_paga <- RunPAGA(procambium_paga, group_by="trajectory_clusters", linear_reduction="integrated.rpca", nonlinear_reduction="umap.rpca", n_pcs=45, n_neighbors=70, infer_pseudotime=TRUE, root_cell=root_cell)

pdf("/home/moliva/root_sc_paper/analyses/trajectories/procambium_paga.pdf", width=6, height=6)
PAGAPlot(procambium_paga, reduction="umap.rpca", edge_threshold=0.05)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/procambium_pseudotime.pdf", width=6, height=6)
FeatureDimPlot(procambium_paga, features = "dpt_pseudotime", reduction="umap.rpca")
dev.off()

detach('package:SCP', unload=TRUE)
detach('package:Seurat', unload=TRUE)
library(Seurat)

obj@meta.data[row.names(procambium@meta.data),"trajectory_clusters"] <- paste0("procambium_",procambium$trajectory_clusters)
obj@meta.data[row.names(procambium_paga@meta.data),"pseudotime"] <- procambium_paga$dpt_pseudotime
obj@meta.data[WhichCells(procambium, idents=c(0,3)), "identity"] <- "procambium ERS+"
obj@meta.data[WhichCells(procambium, idents=c(1,2)), "identity"] <- "procambium ERS-"
rm(procambium_paga)
rm(procambium)


#phloem

phloem <- subset(obj, idents="phloem")

phloem <- FindNeighbors(phloem, reduction = "integrated.rpca", dims = 1:45, k.param=70, verbose=F)
phloem <- FindClusters(phloem, resolution = 0.3, cluster.name = "trajectory_clusters", verbose=F)
pdf("/home/moliva/root_sc_paper/analyses/trajectories/phloem_clusters.pdf", width=6, height=6)
DimPlot(phloem, group.by="trajectory_clusters", label=TRUE)+NoLegend()
dev.off()


root_cell <- row.names(phloem@meta.data[phloem$early_epi_genes_score==max(phloem$early_epi_genes_score),])

phloem_paga <- phloem
phloem_paga[["RNA"]] <- as(phloem[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)

phloem_paga <- RunPAGA(phloem_paga, group_by="trajectory_clusters", linear_reduction="integrated.rpca", nonlinear_reduction="umap.rpca", n_pcs=45, n_neighbors=70, infer_pseudotime=TRUE, root_cell=root_cell)

pdf("/home/moliva/root_sc_paper/analyses/trajectories/phloem_paga.pdf", width=6, height=6)
PAGAPlot(phloem_paga, reduction="umap.rpca", edge_threshold=0.1)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/phloem_pseudotime.pdf", width=6, height=6)
FeatureDimPlot(phloem_paga, features = "dpt_pseudotime", reduction="umap.rpca")
dev.off()

detach('package:SCP', unload=TRUE)
detach('package:Seurat', unload=TRUE)
library(Seurat)

obj@meta.data[row.names(phloem@meta.data),"trajectory_clusters"] <- paste0("phloem_",phloem$trajectory_clusters)
obj@meta.data[row.names(phloem_paga@meta.data),"pseudotime"] <- phloem_paga$dpt_pseudotime
obj@meta.data[WhichCells(phloem, idents=c(0,2)), "identity"] <- "metaphloem ERS+"
obj@meta.data[WhichCells(phloem, idents=c(1,3)), "identity"] <- "metaphloem ERS-"
obj@meta.data[WhichCells(phloem, idents=4), "identity"] <- "protophloem ERS+"
obj@meta.data[WhichCells(phloem, idents=5), "identity"] <- "protophloem ERS-"
rm(phloem_paga)
rm(phloem)


#xylem

xylem <- subset(obj, idents="xylem")

xylem <- FindNeighbors(xylem, reduction = "integrated.rpca", dims = 1:45, k.param=40, verbose=F)
xylem <- FindClusters(xylem, resolution = 0.15, cluster.name = "trajectory_clusters", verbose=F)
pdf("/home/moliva/root_sc_paper/analyses/trajectories/xylem_clusters.pdf", width=6, height=6)
DimPlot(xylem, group.by="trajectory_clusters", label=TRUE)+NoLegend()
dev.off()


root_cell <- row.names(xylem@meta.data[xylem$early_epi_genes_score==max(xylem$early_epi_genes_score),])

xylem_paga <- xylem
xylem_paga[["RNA"]] <- as(xylem[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)

xylem_paga <- RunPAGA(xylem_paga, group_by="trajectory_clusters", linear_reduction="integrated.rpca", nonlinear_reduction="umap.rpca", n_pcs=45, n_neighbors=40, infer_pseudotime=TRUE, root_cell=root_cell)

pdf("/home/moliva/root_sc_paper/analyses/trajectories/xylem_paga.pdf", width=6, height=6)
PAGAPlot(xylem_paga, reduction="umap.rpca", edge_threshold=0.05)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/xylem_pseudotime.pdf", width=6, height=6)
FeatureDimPlot(xylem_paga, features = "dpt_pseudotime", reduction="umap.rpca")
dev.off()

detach('package:SCP', unload=TRUE)
detach('package:Seurat', unload=TRUE)
library(Seurat)

obj@meta.data[row.names(xylem@meta.data),"trajectory_clusters"] <- paste0("xylem_",xylem$trajectory_clusters)
obj@meta.data[row.names(xylem_paga@meta.data),"pseudotime"] <- xylem_paga$dpt_pseudotime
obj@meta.data[WhichCells(xylem, idents=1), "identity"] <- "metaxylem ERS+"
obj@meta.data[WhichCells(xylem, idents=0), "identity"] <- "metaxylem ERS-"
obj@meta.data[WhichCells(xylem, idents=4), "identity"] <- "protoxylem ERS+"
obj@meta.data[WhichCells(xylem, idents=2), "identity"] <- "protoxylem ERS-"
obj@meta.data[WhichCells(xylem, idents=3), "identity"] <- "protoxylem early"
rm(xylem_paga)
rm(xylem)

#Final Plots

pdf("/home/moliva/root_sc_paper/analyses/trajectories/umap_col_trajectory_clusters.pdf", width=6, height=6)
DimPlot(obj, group.by="trajectory_clusters")+NoLegend()
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/umap_col_identities.pdf", width=9, height=6)
DimPlot(obj, group.by="identity")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/trajectories/umap_col_pseudotime.pdf", width=7, height=6)
FeaturePlot(obj,"pseudotime")
dev.off()

saveRDS(obj,"/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")