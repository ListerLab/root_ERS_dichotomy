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
