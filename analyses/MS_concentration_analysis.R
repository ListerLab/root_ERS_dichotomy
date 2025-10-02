library(Seurat)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 5)

obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

Idents(obj)<- "orig.ident"
obj <- subset(obj, idents=c("RL2757_col_fullMS_r1","RL2758_col_fullMS_r2","RL4207_col_fullMS_r3","RL2759_col_halfMS_r1","RL2760_col_halfMS_r2","RL4208_col_halfMS_r3"))
obj@meta.data[WhichCells(obj, idents=c("RL2757_col_fullMS_r1","RL2758_col_fullMS_r2","RL4207_col_fullMS_r3")),"condition"] <- "fullMS"
obj@meta.data[WhichCells(obj, idents=c("RL2759_col_halfMS_r1","RL2760_col_halfMS_r2","RL4208_col_halfMS_r3")),"condition"] <- "halfMS"

pdf("/home/moliva/root_sc_paper/analyses/MS/umap_MS_concentration.pdf", width=13, height=6)
DimPlot(obj, group.by="identity", split.by="condition")
dev.off()

group1 <- obj@meta.data %>% subset(condition=="halfMS") %>% rownames
group2 <- obj@meta.data %>% subset(condition=="fullMS") %>% rownames
dscores <- FindDMEs(obj, features='ModuleScores', barcodes1=group1, barcodes2=group2, test.use='wilcox')
pdf("/home/moliva/root_sc_paper/analyses/MS/lollipop_halfMSvsfullMS.pdf",width=6,height=10)
PlotDMEsLollipop(obj, subset(dscores, abs(avg_log2FC) >0.25), pvalue='p_val_adj', wgcna_name="root_col_modules")
dev.off()

ids <- grep("ERS",obj$identity,value=TRUE)
Idents(obj) <- "identity"
obj_sub <- subset(obj, idents = ids)
obj_sub$celltype <- gsub(" ERS[+-]", "", obj_sub$identity)
obj_sub$state <- gsub(".*(ERS[+-])", "\\1", obj_sub$identity)

library(sccomp)

obj_sub$condition <- factor(obj_sub$condition, levels=c("fullMS","halfMS"))
sccomp_res <- sccomp_estimate(obj_sub, formula_composition=~condition, .sample=orig.ident, .cell_group=state, cores=5) %>% sccomp_test
pdf("/home/moliva/root_sc_paper/analyses/MS/boxplot_ERS_total.pdf", width=3, height=6)
sccomp_res %>% sccomp_boxplot(factor="condition")
dev.off()
saveRDS(sccomp_res,"/home/moliva/root_sc_paper/analyses/MS/sccomp_res_total.rds")


Idents(obj_sub) <- "celltype"
for (i in unique(obj_sub$celltype)) {
	
	sub <- subset(obj_sub, idents=i)
	res <- sccomp_estimate(sub, formula_composition=~condition, .sample=orig.ident, .cell_group=state, cores=5) %>% sccomp_test
	pdf(paste0("/home/moliva/root_sc_paper/analyses/MS/boxplot_ERS_",i,".pdf"), width=3, height=6)
	print(res %>% sccomp_boxplot(factor="condition"))
	dev.off()
	saveRDS(res,paste0("/home/moliva/root_sc_paper/analyses/MS/sccomp_res_",i,".rds"))
	
}

