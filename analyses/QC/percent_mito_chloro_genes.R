library(Seurat)
library(ggplot2)


obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")
obj$percent.mito <- PercentageFeatureSet(obj,pattern="^ATM")
obj$percent.chloro <- PercentageFeatureSet(obj,pattern="^ATC")

obj$state <- "no ERS"
obj@meta.data[row.names(subset(obj@meta.data, identity %in% grep("ERS+",unique(obj$identity),value=TRUE))),"state"] <- "ERS+"
obj@meta.data[row.names(subset(obj@meta.data, identity %in% grep("ERS-",unique(obj$identity),value=TRUE))),"state"] <- "ERS-"
Idents(obj) <- "state"

pdf("/home/moliva/root_sc_paper/analyses/QC/mito_chloro_genes/percent_mito_per_state.pdf")
VlnPlot(obj,"percent.mito")+geom_hline(yintercept=5,linetype="dashed")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/QC/mito_chloro_genes/percent_chloro_per_state.pdf")
VlnPlot(obj,"percent.chloro")+geom_hline(yintercept=5,linetype="dashed")
dev.off()

