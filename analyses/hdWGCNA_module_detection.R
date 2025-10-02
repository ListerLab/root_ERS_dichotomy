library(Seurat)
library(WGCNA)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 20)
library(MetBrewer)

obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")
genes <- SetupForWGCNA(obj, gene_select="fraction", fraction=0.1, group.by="lineage", wgcna_name="gene_selection")
genes <- genes@misc$gene_selection$wgcna_genes
genes <- intersect(genes, VariableFeatures(obj))
obj <- SetupForWGCNA(obj, gene_select="custom", gene_list=genes, wgcna_name="root_col_modules")
obj <- MetacellsByGroups(seurat_obj=obj, group.by= c("lineage","orig.ident"), reduction="integrated.rpca",ident.group="lineage",k=20, max_shared=5)
obj <- NormalizeMetacells(obj)
obj <- SetMultiExpr(obj, assay="RNA", layer="data", multi.group.by="orig.ident", multi_groups=NULL)
obj <- ConstructNetwork(obj, soft_power=c(8,8,8,7,8,8,6,7,7), minModuleSize=30, consensus=TRUE, overwrite_tom=TRUE, deepSplit=2)

#Rename modules
 modules <- GetModules(obj)
 new_names <- c("Module 23","Module 22","Module 21","Module 11","Module 3","Module 4","Module 24","Module 2","Module 25","Module 1","Module 20","Module 5","Module 9","Module 8","Module 19","Module 7","Module 6","Module 13","Module 14","Module 10","Module 15","Module 17","Module 16","Module 12","Module 18")
 old_mods <- levels(modules$module)
 grey_ind <- which(old_mods == 'grey')
 new_names <- c(new_names[1:(grey_ind-1)], 'grey', new_names[grey_ind:length(new_names)])
 new_mod_df <- data.frame(old = old_mods ,new = new_names)
 modules$module <- factor(new_mod_df[match(modules$module, new_mod_df$old),'new'],levels = c(paste0("Module ",1:25),"grey"))
 obj <- SetModules(obj, modules)

#Change Module colours
 mod_color_df <- GetModules(obj)%>%
	dplyr::select(c(module, color)) %>%
 	distinct %>% 
 	arrange(module)
 n_mods <- nrow(mod_color_df) - 1
 new_colors <- paste0(met.brewer("Austria", n=n_mods,type="continuous"))
 obj <- ResetModuleColors(obj, new_colors)

obj <- ModuleExprScore(obj, n_genes="all", method="Seurat")
plot_list <- ModuleFeaturePlot(obj, features="scores", order=TRUE, reduction="umap.rpca")
pdf("/home/moliva/root_sc_paper/analyses/hdWGCNA/module_scores.pdf", width=36, height=30)
wrap_plots(plot_list, ncol=6)
dev.off()

obj <- AvgModuleExpr(obj, n_genes="all")
plot_list <- ModuleFeaturePlot(obj, features="average", order=TRUE, reduction="umap.rpca")
pdf("/home/moliva/root_sc_paper/analyses/hdWGCNA/module_average.pdf", width=36, height=30)
wrap_plots(plot_list, ncol=6)
dev.off()

obj <- ModuleEigengenes(obj, group.by.vars="orig.ident")
pdf("/home/moliva/root_sc_paper/analyses/hdWGCNA/hME.pdf", width=36, height=30)
wrap_plots(plot_list, ncol=6)
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/hdWGCNA/hist_genes_per_module.pdf", width=10, height=10)
modules %>% filter(module != "grey") %>% count(module) %>% ggplot(aes(x=module, y=n)) + geom_col() + theme_classic() + scale_x_discrete(labels=paste0("M",1:25))
dev.off()

modules <- GetModules(obj)
write.table(modules,"/home/moliva/root_sc_paper/analyses/hdWGCNA/modules.tsv", quote=FALSE, sep='\t', row.names=FALSE)

saveRDS(obj,"/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")