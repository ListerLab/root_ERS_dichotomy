library(Seurat)
library(ggplot2)


drought_100 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_100percentH2O/integrated_drought_100percentH2O_predicted_id.rds")
drought_50 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_50percentH2O/integrated_drought_50percentH2O_predicted_id.rds")
drought_20 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/integrated_drought_20percentH2O_predicted_id.rds")


ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")
ref <- RunUMAP(ref, reduction="integrated.rpca", dims=1:45, reduction.name="umap.rpca",verbose=F, return.model=TRUE)

anchors <- FindTransferAnchors(reference=ref, query=drought_100, dims=1:50, reference.reduction="integrated.rpca")
drought_100 <- MapQuery(anchorset=anchors, reference=ref, query=drought_100, reference.reduction="integrated.rpca", reduction.model="umap.rpca")

cols <- c("#f27457","#e02228","#7f7ebb","#43459a","#ef69a0","#db3b94","#3fb5c3","#388fc1","#426db4","#2a3684","#bf822b","#8b5422","#a4d6b5","#3dab5e","#99d4c7","#fddac1","#f9b9a0","#f2ecab","#f5e609","#f6e2b9","#dec17c","#399872","#54916e","#0f6f39","#eeb746","#cc4f27","#ef5331","#c1d8b1","#a5d06c")

p1 <- DimPlot(ref, reduction="umap.rpca",group.by="identity")+NoLegend()+ggtitle("Reference Atlas")
p2 <- DimPlot(drought_100, reduction="ref.umap", group.by="predicted.identity", cols=cols)+ ggtitle("drought 100%H2O snRNA-seq projection")
pdf("/home/moliva/root_sc_paper/analyses/drought/drought_100percentH2O_umap_ref_projection.pdf", width=16, height=6)
p1+p2
dev.off()


anchors <- FindTransferAnchors(reference=ref, query=drought_50, dims=1:50, reference.reduction="integrated.rpca")
drought_50 <- MapQuery(anchorset=anchors, reference=ref, query=drought_50, reference.reduction="integrated.rpca", reduction.model="umap.rpca")
p2 <- DimPlot(drought_50, reduction="ref.umap", group.by="predicted.identity", cols=cols)+ ggtitle("drought 50%H2O snRNA-seq projection")
pdf("/home/moliva/root_sc_paper/analyses/drought/drought_50percentH2O_umap_ref_projection.pdf", width=16, height=6)
p1+p2
dev.off()


anchors <- FindTransferAnchors(reference=ref, query=drought_20, dims=1:50, reference.reduction="integrated.rpca")
drought_20 <- MapQuery(anchorset=anchors, reference=ref, query=drought_20, reference.reduction="integrated.rpca", reduction.model="umap.rpca")
p2 <- DimPlot(drought_20, reduction="ref.umap", group.by="predicted.identity", cols=cols)+ ggtitle("drought 20%H2O snRNA-seq projection")
pdf("/home/moliva/root_sc_paper/analyses/drought/drought_20percentH2O_umap_ref_projection.pdf", width=16, height=6)
p1+p2
dev.off()
