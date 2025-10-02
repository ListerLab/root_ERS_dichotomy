library(Seurat)
library(sccomp)
library(dplyr)
  
drought_100 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_100percentH2O/integrated_drought_100percentH2O_predicted_id.rds")
drought_50 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_50percentH2O/integrated_drought_50percentH2O_predicted_id.rds")
drought_20 <- readRDS("/home/moliva/root_sc_paper/data/integrated/drought_20percentH2O/integrated_drought_20percentH2O_predicted_id.rds")

drought_100$H2Ocontent <- 1
drought_50$H2Ocontent <- 0.5
drought_20$H2Ocontent <- 0.2

obj <- merge(drought_100, list(drought_20, drought_50))
colnames(obj@meta.data)[6] <- "identity"

ids <- grep("ERS",obj$identity,value=TRUE)
Idents(obj) <- "identity"
obj_sub <- subset(obj, idents = ids)
obj_sub$celltype <- gsub(" ERS[+-]", "", obj_sub$identity)
obj_sub$state <- gsub(".*(ERS[+-])", "\\1", obj_sub$identity)

sccomp_res <- sccomp_estimate(obj_sub, formula_composition=~H2Ocontent, .sample=orig.ident, .cell_group=state, cores=5) %>% sccomp_test

pdf("/home/moliva/root_sc_paper/analyses/drought/composition_plot_ERS_total.pdf", width=3, height=6)
sccomp_res %>% plot
dev.off()

saveRDS(sccomp_res,"/home/moliva/root_sc_paper/analyses/drought/sccomp_res_total.rds")


Idents(obj_sub) <- "celltype"
for (i in unique(obj_sub$celltype)) {
	
	sub <- subset(obj_sub, idents=i)
	res <- sccomp_estimate(sub, formula_composition=~H2Ocontent, .sample=orig.ident, .cell_group=state, cores=5) %>% sccomp_test
	pdf(paste0("/home/moliva/root_sc_paper/analyses/drought/composition_plot_ERS_",i,".pdf"), width=3, height=6)
	print(res %>% plot)
	dev.off()
	saveRDS(res,paste0("/home/moliva/root_sc_paper/analyses/drought/sccomp_res_",i,".rds"))
	
}
