library(Seurat)
library(sccomp)
library(dplyr)

col <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")
Idents(col) <- "orig.ident"
Idents(col) <- "orig.ident"
col <- subset(col, idents=c("RL650_col_cvi_r3","RL970_col_ler_r4","RL971_col_c24_r5","RL2757_col_fullMS_r1","RL2758_col_fullMS_r2","RL4207_col_fullMS_r3"))
col$genotype <- "col"

c24 <- readRDS("/home/moliva/root_sc_paper/data/integrated/c24/integrated_c24_predicted_id.rds")
colnames(c24@meta.data)[6] <- "lineage"
colnames(c24@meta.data)[7] <- "identity"
c24$genotype <- "c24"

cvi <- readRDS("/home/moliva/root_sc_paper/data/integrated/cvi/integrated_cvi_predicted_id.rds")
colnames(cvi@meta.data)[6] <- "lineage"
colnames(cvi@meta.data)[7] <- "identity"
cvi$genotype <- "cvi"

ler <- readRDS("/home/moliva/root_sc_paper/data/integrated/ler/integrated_ler_predicted_id.rds")
colnames(ler@meta.data)[5] <- "lineage"
colnames(ler@meta.data)[6] <- "identity"
ler$genotype <- "ler"

ws2 <- readRDS("/home/moliva/root_sc_paper/data/integrated/ws2/integrated_ws2_predicted_id.rds")
colnames(ws2@meta.data)[6] <- "lineage"
colnames(ws2@meta.data)[7] <- "identity"
ws2$genotype <- "ws2"

obj <- merge(col, list(c24,cvi,ler,ws2))
rm(col,c24,cvi,ler,ws2)

ids <- grep("ERS",obj$identity,value=TRUE)
Idents(obj) <- "identity"
obj_sub <- subset(obj, idents = ids)
obj_sub$celltype <- gsub(" ERS[+-]", "", obj_sub$identity)
obj_sub$state <- gsub(".*(ERS[+-])", "\\1", obj_sub$identity)

obj_sub$genotype <- factor(obj_sub$genotype, levels=c("col","c24","cvi","ler","ws2"))
sccomp_res <- sccomp_estimate(obj_sub, formula_composition=~genotype, .sample=orig.ident, .cell_group=state, cores=5) %>% sccomp_test
pdf("/home/moliva/root_sc_paper/analyses/ecotype/ecotype_boxplot_ERS_total.pdf", width=3, height=6)
sccomp_res %>% sccomp_boxplot(factor="genotype")
dev.off()
saveRDS(sccomp_res,"/home/moliva/root_sc_paper/analyses/ecotype/sccomp_res_total.rds")
