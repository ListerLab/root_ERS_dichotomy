library(Seurat)
library(hdWGCNA)
library(matrixStats)
library(ggplot2)


obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

modules <- GetModules(obj)
genes <- row.names(subset(modules, module %in% c("Module 13","Module 14","Module 20","Module 21")))

Idents(obj) <- "lineage"
obj <- subset(obj, idents = c("xpp"))
obj$state <- gsub(".*(ERS[+-])", "\\1", obj$identity)
Idents(obj) <- "state"

data <- GetAssayData(obj, slot="data")
data <- data[genes,]

avg_ERS_pos <- rowMeans(as.matrix(data[,WhichCells(obj, idents="ERS+")]))
avg_ERS_neg <- rowMeans(as.matrix(data[,WhichCells(obj, idents="ERS-")]))
delta <- (avg_ERS_pos - avg_ERS_neg)

gene_data <- data.frame(delta=delta, avg_ERS_pos=avg_ERS_pos, module=modules[genes,"module"])
gene_data[row.names(subset(modules, module %in% c("Module 13","Module 14"))),"type"] <- "symmetric"
gene_data[row.names(subset(modules, module %in% c("Module 20","Module 21"))),"type"] <- "asymmetric"

pdf("/home/moliva/root_sc_paper/analyses/QC/expression_levels_ERS/expression_ERS_neg_pos.pdf", width=8, height=6)
ggplot(gene_data,aes(x=avg_ERS_pos,y=delta, color=type, shape=module))+geom_point()+theme_classic()+geom_smooth(method='lm', aes(group=type))+scale_x_continuous(transform='sqrt')
dev.off()