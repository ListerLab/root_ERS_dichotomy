library(Seurat)
library(hdWGCNA)
library(ggplot2)
library(patchwork)


rep1 <- readRDS("/home/moliva/root_sc_paper/published_datasets/wendrich/seurat_objects/wendrich_replicate1.rds")
rep1$orig.ident <- "rep1"
mat <- Read10X_h5("/home/moliva/root_sc_paper/published_datasets/wendrich/GSE141730_RAW/GSM4212550_BDR3_filtered_gene_bc_matrices.h5")
rep1$filtered <- "yes"
rep1@meta.data[intersect(colnames(mat),colnames(rep1)),"filtered"] <- "no"

rep2 <- readRDS("/home/moliva/root_sc_paper/published_datasets/wendrich/seurat_objects/wendrich_replicate2.rds")
rep2$orig.ident <- "rep2"
mat <- Read10X_h5("/home/moliva/root_sc_paper/published_datasets/wendrich/GSE141730_RAW/GSM4212551_BDR4_filtered_gene_bc_matrices.h5")
rep2$filtered <- "yes"
rep2@meta.data[intersect(colnames(mat),colnames(rep2)),"filtered"] <- "no"

rep3 <- readRDS("/home/moliva/root_sc_paper/published_datasets/wendrich/seurat_objects/wendrich_replicate3.rds")
rep3$orig.ident <- "rep3"
mat <- Read10X_h5("/home/moliva/root_sc_paper/published_datasets/wendrich/GSE141730_RAW/GSM4212552_BDR5_filtered_gene_bc_matrices.h5")
rep3$filtered <- "yes"
rep3@meta.data[intersect(colnames(mat),colnames(rep3)),"filtered"] <- "no"

obj <- merge(rep1, list(rep2,rep3))

obj <- NormalizeData(obj)
obj <- ScaleData(obj)

obj <- FindVariableFeatures(obj, nfeatures=25000,verbose=F)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj), value=TRUE)
chloro <- grep("ATC", row.names(obj), value=TRUE)
VariableFeatures(obj) <- setdiff(VariableFeatures(obj),c(proto$V1,mito,chloro))

obj <- RunPCA(obj, npcs=50,verbose=F)

obj <- IntegrateLayers(obj, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50, verbose=FALSE)

ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")
obj <- JoinLayers(obj)
anchors <- FindTransferAnchors(reference=ref, query=obj, dims=1:50, reference.reduction="integrated.rpca")
predictions_lineage <- TransferData(anchorset=anchors, refdata=ref$lineage, dims=1:50)
obj <-AddMetaData(obj, metadata=predictions_lineage[,1,drop=FALSE], col.name="predicted.lineage")
predictions_identity <- TransferData(anchorset=anchors, refdata=ref$identity, dims=1:50)
obj <-AddMetaData(obj, metadata=predictions_identity[,1,drop=FALSE], col.name="predicted.identity")

obj <- RunUMAP(obj, reduction="integrated.rpca", dims=1:50, reduction.name="umap.rpca",verbose=F)

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/wendrich/wendrich_umap_integrated.pdf", width=9, height=6)
DimPlot(obj, group.by="predicted.identity")
dev.off()

ref <- RunUMAP(ref, reduction="integrated.rpca", dims=1:45, reduction.name="umap.rpca",verbose=F, return.model=TRUE)
obj <- MapQuery(anchorset=anchors, reference=ref, query=obj, refdata=list(identity="identity"), reference.reduction="integrated.rpca", reduction.model="umap.rpca")


p1 <- DimPlot(ref, reduction="umap.rpca",group.by="identity")+NoLegend()+ggtitle("Reference Atlas")
p2 <- DimPlot(obj, reduction="ref.umap", group.by="predicted.identity")+ ggtitle("snRNA-seq projection")
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/wendrich/wendrich_umap_ref_projection.pdf", width=16, height=6)
p1+p2
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/wendrich/wendrich_umap_ref_projection_filtered_cells.pdf", width=7, height=6)
DimPlot(obj, reduction="ref.umap", group.by="filtered")
dev.off()


obj <- ProjectModules(obj, seurat_ref= ref, wgcna_name="root_col_modules", wgcna_name_proj="projected")
obj <- ModuleExprScore(obj, n_genes="all", method="Seurat")

plot_list <- ModuleFeaturePlot(obj, features="scores", order=TRUE, reduction="ref.umap")
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/wendrich/projected_module_scores.pdf", width=36, height=30)
wrap_plots(plot_list, ncol=6)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/wendrich/wendrich_umap_ref_projection_umi.pdf", width=7, height=6)
FeaturePlot(obj, "nCount_RNA", min.cutoff="q10", max.cutoff="q90", reduction="ref.umap")
dev.off()


scores <- GetModuleScores(obj)
obj$Module_20 <- scores$'Module 20'
obj$Module_21 <- scores$'Module 21'

identity_order <- c("root cap early", "columella","root cap tip","border cells","root cap lateral","epidermis early","hair cells early","hair cells ERS+","hair cells ERS-","non-hair cells","cortex ERS+","cortex ERS-","endodermis ERS+", "endodermis ERS-","ppp ERS+","ppp ERS-","xpp ERS+","xpp ERS-","procambium ERS+","procambium ERS-","metaphloem ERS+","metaphloem ERS-","protophloem ERS+","protophloem ERS-","metaxylem ERS+","metaxylem ERS-","protoxylem early","protoxylem ERS+","protoxylem ERS-")
Idents(obj) <- "predicted.identity"
Idents(obj) <- factor(Idents(obj), levels=identity_order)

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/wendrich/wendrich_vlnplot_module20_scores.pdf", width=10, height=6)
VlnPlot(obj,"Module_20")
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/wendrich/wendrich_vlnplot_module21_scores.pdf", width=10, height=6)
VlnPlot(obj,"Module_21")
dev.off()
