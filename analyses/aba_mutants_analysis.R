library(Seurat)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 5)

ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")
Idents(ref) <- "orig.ident"
ref_sub <- subset(ref, idents=c("RL650_col_cvi_r3","RL970_col_ler_r4","RL971_col_c24_r5","RL2757_col_fullMS_r1","RL2758_col_fullMS_r2","RL4207_col_fullMS_r3"))

aba <- readRDS("/home/moliva/root_sc_paper/data/integrated/aba2_3/integrated_aba2_3_predicted_id.rds")
aba$lineage <- aba$predicted.lineage
aba$identity <- aba$predicted.identity
aba$genotype <- "aba2_3"

pyl <- readRDS("/home/moliva/root_sc_paper/data/integrated/pyl11458/integrated_pyl11458_predicted_id.rds")
pyl$lineage <- pyl$predicted.lineage
pyl$identity <- pyl$predicted.identity
pyl$genotype <- "pyl11458"

ref_sub$genotype <- "col"
obj <- merge(ref_sub,list(aba,pyl))
rm(ref_sub)

obj <- JoinLayers(obj)
obj <- ScaleData(obj)
obj <- SetupForWGCNA(obj, gene_select="custom", gene_list=ref@misc$root_col_modules$wgcna_genes, wgcna_name="root_col_modules")
modules <- GetModules(ref)
obj <- SetModules(obj,modules)
obj <- ModuleExprScore(obj, n_genes="all", method="Seurat")


group1 <- obj@meta.data %>% subset(genotype=="aba2_3") %>% rownames
group2 <- obj@meta.data %>% subset(genotype=="col") %>% rownames
dscores <- FindDMEs(obj, features='ModuleScores', barcodes1=group1, barcodes2=group2, test.use='wilcox')
pdf("/home/moliva/root_sc_paper/analyses/mutants/lollipop_aba2_3.pdf",width=6,height=10)
PlotDMEsLollipop(obj, dscores, pvalue='p_val_adj', wgcna_name="root_col_modules")
dev.off()

group1 <- obj@meta.data %>% subset(genotype=="pyl11458") %>% rownames
dscores <- FindDMEs(obj, features='ModuleScores', barcodes1=group1, barcodes2=group2, test.use='wilcox')
pdf("/home/moliva/root_sc_paper/analyses/mutants/lollipop_pyl11458.pdf",width=6,height=10)
PlotDMEsLollipop(obj, dscores, pvalue='p_val_adj', wgcna_name="root_col_modules")
dev.off()


ids <- grep("ERS",obj$identity,value=TRUE)
Idents(obj) <- "identity"
obj_sub <- subset(obj, idents = ids)
obj_sub$celltype <- gsub(" ERS[+-]", "", obj_sub$identity)
obj_sub$state <- gsub(".*(ERS[+-])", "\\1", obj_sub$identity)

library(sccomp)

obj_sub$genotype <- factor(obj_sub$genotype, levels=c("col","aba2_3","pyl11458"))
sccomp_res <- sccomp_estimate(obj_sub, formula_composition=~genotype, .sample=orig.ident, .cell_group=state, cores=5) %>% sccomp_test
pdf("boxplot_ERS_total.pdf", width=3, height=6)
sccomp_res %>% sccomp_boxplot(factor="genotype")
dev.off()
saveRDS(sccomp_res,"/home/moliva/root_sc_paper/analyses/mutants/sccomp_res_total.rds")


Idents(obj_sub) <- "celltype"
for (i in unique(obj_sub$celltype)) {
	
	sub <- subset(obj_sub, idents=i)
	res <- sccomp_estimate(sub, formula_composition=~genotype, .sample=orig.ident, .cell_group=state, cores=5) %>% sccomp_test
	pdf(paste0("boxplot_ERS_",i,".pdf"), width=3, height=6)
	print(res %>% sccomp_boxplot(factor="genotype"))
	dev.off()
	saveRDS(res,paste0("/home/moliva/root_sc_paper/analyses/mutants/sccomp_res_",i,".rds"))
	
}


























regulon_scores <- GetRegulonScores(obj, target_type = "positive")
reg_assay <- Seurat::CreateAssayObject(t(regulon_scores))
cur_dregs <- FindMarkers(
            reg_assay,
            cells.1 = group1,
            cells.2 = group2,
            slot='counts', # should I make this layer at some point?
            test.use="wilcox",
            only.pos=FALSE,
            logfc.threshold=0,
            min.pct=0,
            verbose=FALSE,
            pseudocount.use=0)


obj <- ProjectModules( seurat_obj=obj, seurat_ref=ref, group.by.vars="orig.ident", wgcna_name_proj="projected")
obj <- ModuleExprScore(obj, n_genes="all", method="Seurat")









ref <- RegulonScores(ref, target_type="positive",ncores=20)
pos_regulon_scores <- GetRegulonScores(ref, target_type="positive")
pos_regulon_scores$identity <- ref@meta.data[row.names(pos_regulon_scores),"identity"]
avg_regulon_scores <- pos_regulon_scores %>% group_by(identity) %>% summarize(across(everything(),mean))
avg_regulon_scores <- as.data.frame(avg_regulon_scores)
row.names(avg_regulon_scores) <- avg_regulon_scores$identity
avg_regulon_scores$identity <- NULL

genenames <- read.csv("/home/moliva/root_sc_paper/genome/gene_descriptions_tair.txt", sep = "\t")
genenames$Primary.Gene.Symbol <- gsub(".*\\((.*)\\).*", "\\1", genenames$Primary.Gene.Symbol)
genenames[genenames$Primary.Gene.Symbol=="","Primary.Gene.Symbol"] <- as.character(genenames[genenames$Primary.Gene.Symbol=="","Locus.Identifier"])
genenames$Primary.Gene.Symbol <- as.character(genenames$Primary.Gene.Symbol)
genenames$Locus.Identifier <- as.character(genenames$Locus.Identifier)
colnames(_regulon_scores) <- genenames[match(colnames(pos_regulon_scores), genenames$Locus.Identifier),"Primary.Gene.Symbol"]

identity_order <- c("root cap early", "columella","root cap tip","border cells","root cap lateral","epidermis early","hair cells early","hair cells ERS+","hair cells ERS-","non-hair cells","cortex ERS+","cortex ERS-","endodermis ERS+", "endodermis ERS-","ppp ERS+","ppp ERS-","xpp ERS+","xpp ERS-","procambium ERS+","procambium ERS-","metaphloem ERS+","metaphloem ERS-","protophloem ERS+","protophloem ERS-","metaxylem ERS+","metaxylem ERS-","protoxylem early","protoxylem ERS+","protoxylem ERS-")
pdf("/home/moliva/root_sc_paper/analyses/TFs/heatmap_avg_regulon_scores_test.pdf",width=24, height=6)
pheatmap(avg_regulon_scores[identity_order,], cluster_rows=FALSE,cluster_cols=TRUE, scale="column",cutree_cols=8) 
dev.off()



