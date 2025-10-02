library(Seurat)
library(sccomp)
library(dplyr)
library(ggplot2)

ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")
ref$genotype <- "col"
scr <- readRDS("/home/moliva/root_sc_paper/data/integrated/scr3/integrated_scr3_predicted_id.rds")
scr$genotype <- "scr3"
colnames(scr@meta.data)[7] <- "identity"

scr <- NormalizeData(scr)
scr <- ScaleData(scr)
scr <- FindVariableFeatures(scr, nfeatures=25000,verbose=F)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(scr), value=TRUE)
chloro <- grep("ATC", row.names(scr), value=TRUE)
VariableFeatures(scr) <- setdiff(VariableFeatures(scr),c(proto$V1,mito,chloro))
scr <- RunPCA(scr, npcs=50,verbose=F)
scr <- RunUMAP(scr, dims=1:50, verbose=F)

pdf("/home/moliva/root_sc_paper/analyses/scr3/umap_scr_identities.pdf",width=9, height=6)
DimPlot(scr, group.by="identity")
dev.off()

obj <- merge(ref,scr)

obj$identity <- gsub(" ERS\\+$","", obj$identity)
obj$identity <- gsub(" ERS-$","", obj$identity)
obj$identity <- gsub(" early","", obj$identity)
obj$identity <- gsub(" lateral","", obj$identity)
obj$identity <- gsub(" tip","", obj$identity)
obj$identity <- gsub("border cells","root cap", obj$identity)
obj$identity <- gsub("non-hair cells","epidermis", obj$identity)
obj$identity <- gsub("hair cells","epidermis", obj$identity)


obj$genotype <- factor(obj$genotype, levels=c("col","scr3"))
sccomp_res <- sccomp_estimate(obj, formula_composition=~genotype, .sample=orig.ident, .cell_group=identity, cores=5) %>% sccomp_test
pdf("/home/moliva/root_sc_paper/analyses/scr3/boxplot_lineage_comp.pdf")
sccomp_res %>% sccomp_boxplot(factor="genotype")
dev.off()
saveRDS(sccomp_res,"/home/moliva/root_sc_paper/analyses/scr3/sccomp_res_lineages.rds")


anchors <- FindTransferAnchors(reference=ref, query=scr, dims=1:50, reference.reduction="integrated.rpca")
ref <- RunUMAP(ref, reduction="integrated.rpca", dims=1:45, reduction.name="umap.rpca",verbose=F, return.model=TRUE)
scr <- MapQuery(anchorset=anchors, reference=ref, query=scr, refdata=list(identity="identity"), reference.reduction="integrated.rpca", reduction.model="umap.rpca")



p1 <- DimPlot(ref, reduction="umap.rpca",group.by="identity")+NoLegend()+ggtitle("Reference Atlas")
p2 <- DimPlot(scr, reduction="ref.umap", group.by="predicted.identity")+ ggtitle("scr-3 projection")
pdf("/home/moliva/root_sc_paper/analyses/scr3/scr3_umap_ref_projection.pdf", width=16, height=6)
p1+p2
dev.off()
