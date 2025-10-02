library(Seurat)
library(hdWGCNA)
library(ggplot2)
library(patchwork)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep1/filtered_feature_bc_matrix/")
rna1 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep2/filtered_feature_bc_matrix/")
rna2 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep3/filtered_feature_bc_matrix/")
rna3 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep4/filtered_feature_bc_matrix/")
rna4 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

data <- Read10X("/home/moliva/root_sc_paper/published_datasets/farmer/raw_matrices/sNucRNA-seq_rep5/filtered_feature_bc_matrix/")
rna5 <- CreateSeuratObject(counts=data, min.cells=3, min.features=0)

obj.rna <- merge(rna1, list(rna2,rna3,rna4,rna5))

obj.rna <- NormalizeData(obj.rna)
obj.rna <- ScaleData(obj.rna)

obj.rna <- FindVariableFeatures(obj.rna, nfeatures=25000,verbose=F)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj.rna), value=TRUE)
chloro <- grep("ATC", row.names(obj.rna), value=TRUE)
VariableFeatures(obj.rna) <- setdiff(VariableFeatures(obj.rna),c(proto$V1,mito,chloro))

obj.rna <- RunPCA(obj.rna, npcs=50,verbose=F)

obj.rna <- IntegrateLayers(obj.rna, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50, verbose=FALSE)

ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")
obj.rna <- JoinLayers(obj.rna)
anchors <- FindTransferAnchors(reference=ref, query=obj.rna, dims=1:50, reference.reduction="integrated.rpca")
predictions_lineage <- TransferData(anchorset=anchors, refdata=ref$lineage, dims=1:50)
obj.rna <-AddMetaData(obj.rna, metadata=predictions_lineage[,1,drop=FALSE], col.name="predicted.lineage")
predictions_identity <- TransferData(anchorset=anchors, refdata=ref$identity, dims=1:50)
obj.rna <-AddMetaData(obj.rna, metadata=predictions_identity[,1,drop=FALSE], col.name="predicted.identity")

obj.rna <- RunUMAP(obj.rna, reduction="integrated.rpca", dims=1:50, reduction.name="umap.rpca",verbose=F)

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/farmer_snrna_umap_integrated.pdf", width=9, height=6)
DimPlot(obj.rna, group.by="predicted.identity")
dev.off()

ref <- RunUMAP(ref, reduction="integrated.rpca", dims=1:45, reduction.name="umap.rpca",verbose=F, return.model=TRUE)
obj.rna <- MapQuery(anchorset=anchors, reference=ref, query=obj.rna, refdata=list(identity="identity"), reference.reduction="integrated.rpca", reduction.model="umap.rpca")


p1 <- DimPlot(ref, reduction="umap.rpca",group.by="identity")+NoLegend()+ggtitle("Reference Atlas")
p2 <- DimPlot(obj.rna, reduction="ref.umap", group.by="predicted.identity")+ ggtitle("snRNA-seq projection")
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/farmer_snrna_umap_ref_projection.pdf", width=16, height=6)
p1+p2
dev.off()

obj.rna <- ProjectModules(obj.rna, seurat_ref= ref, wgcna_name="root_col_modules", wgcna_name_proj="projected")
obj.rna <- ModuleExprScore(obj.rna, n_genes="all", method="Seurat")

plot_list <- ModuleFeaturePlot(obj.rna, features="scores", order=TRUE, reduction="ref.umap")
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/projected_module_scores.pdf", width=36, height=30)
wrap_plots(plot_list, ncol=6)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/farmer_snrna_umap_ref_projection_umi.pdf", width=7, height=6)
FeaturePlot(obj.rna, "nCount_RNA", min.cutoff="q10", max.cutoff="q90", reduction="ref.umap")
dev.off()


scores <- GetModuleScores(obj.rna)
obj.rna$Module_20 <- scores$'Module 20'
obj.rna$Module_21 <- scores$'Module 21'

identity_order <- c("root cap early", "columella","root cap tip","border cells","root cap lateral","epidermis early","hair cells early","hair cells ERS+","hair cells ERS-","non-hair cells","cortex ERS+","cortex ERS-","endodermis ERS+", "endodermis ERS-","ppp ERS+","ppp ERS-","xpp ERS+","xpp ERS-","procambium ERS+","procambium ERS-","metaphloem ERS+","metaphloem ERS-","protophloem ERS+","protophloem ERS-","metaxylem ERS+","metaxylem ERS-","protoxylem early","protoxylem ERS+","protoxylem ERS-")
Idents(obj.rna) <- "predicted.identity"
Idents(obj.rna) <- factor(Idents(obj.rna), levels=identity_order)

pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/farmer_snrna_vlnplot_module20_scores.pdf", width=10, height=6)
VlnPlot(obj.rna,"Module_20")
dev.off()
pdf("/home/moliva/root_sc_paper/analyses/QC/comparison_published_datasets/farmer/farmer_snrna_vlnplot_module21_scores.pdf", width=10, height=6)
VlnPlot(obj.rna,"Module_21")
dev.off()


