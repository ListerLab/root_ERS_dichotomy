library(Seurat)
library(hdWGCNA)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = 2500 * 1024^2)

obj <- readRDS("/home/moliva/root_sc_paper/published_datasets/shahan/GSE152766_Root_Atlas.rds")

ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

DefaultAssay(obj) <- "RNA"
obj <- UpdateSeuratObject(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f= obj$orig.ident)
obj <- NormalizeData(obj)
obj <- ScaleData(obj)
obj <- JoinLayers(obj)
anchors <- FindTransferAnchors(reference=ref, query=obj, dims=1:50, reference.reduction="integrated.rpca")
predictions_lineage <- TransferData(anchorset=anchors, refdata=ref$lineage, dims=1:50)
obj <-AddMetaData(obj, metadata=predictions_lineage[,1,drop=FALSE], col.name="predicted.lineage")
predictions_identity <- TransferData(anchorset=anchors, refdata=ref$identity, dims=1:50)
obj <-AddMetaData(obj, metadata=predictions_identity[,1,drop=FALSE], col.name="predicted.identity")


pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/shahan_snrna_umap_integrated.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.identity",reduction="umap",raster=FALSE)
dev.off()

obj <- ProjectModules(obj, seurat_ref= ref, wgcna_name="root_col_modules", wgcna_name_proj="projected")
obj <- ModuleExprScore(obj, n_genes="all", method="Seurat")

plot_list <- ModuleFeaturePlot(obj, features="scores", order=TRUE)
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/projected_module_scores.pdf", width=36, height=30)
wrap_plots(plot_list, ncol=6)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/shahan_snrna_umap_umi.pdf", width=7, height=6)
FeaturePlot(obj, "nCount_RNA", min.cutoff="q10", max.cutoff="q90",raster=FALSE)
dev.off()


scores <- GetModuleScores(obj)
obj$Module_20 <- scores$'Module 20'
obj$Module_21 <- scores$'Module 21'

identity_order <- c("root cap early", "columella","root cap tip","border cells","root cap lateral","epidermis early","hair cells early","hair cells ERS+","hair cells ERS-","non-hair cells","cortex ERS+","cortex ERS-","endodermis ERS+", "endodermis ERS-","ppp ERS+","ppp ERS-","xpp ERS+","xpp ERS-","procambium ERS+","procambium ERS-","metaphloem ERS+","metaphloem ERS-","protophloem ERS+","protophloem ERS-","metaxylem ERS+","metaxylem ERS-","protoxylem early","protoxylem ERS+","protoxylem ERS-")
Idents(obj) <- "predicted.identity"
Idents(obj) <- factor(Idents(obj), levels=identity_order)

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/shahan_snrna_vlnplot_module20_scores.pdf", width=10, height=6)
VlnPlot(obj,"Module_20",raster=FALSE)
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/shahan/shahan_snrna_vlnplot_module21_scores.pdf", width=10, height=6)
VlnPlot(obj,"Module_21",raster=FALSE)
dev.off()
