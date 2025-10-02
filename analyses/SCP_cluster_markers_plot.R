library(Seurat)

obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col.rds")

obj_scp <- obj
obj_scp[["RNA"]] <- as(obj[["RNA"]], Class="Assay")

detach('package:Seurat', unload=TRUE)
detach('package:SeuratObject', unload=TRUE)
library(Seurat, lib.loc="/group/ll004/moliva/bin/R_lib_Seuratv4/")
library(SCP)


obj_scp$clusters <- factor(obj_scp$clusters, levels=c(24,19,1,0,2,9,21,5,26,4,11,23,13,25,8,20,3,6,12,15,7,10,14,17,22,16,18))


pdf("/home/moliva/root_sc_paper/analyses/markers/cluster_markers.pdf",width=10, height=6)
FeatureStatPlot(obj_scp, stat.by=c("AT2G04025","AT1G44760","AT1G05010","AT1G26820","AT2G39530","AT4G40090","AT1G12090","AT2G14900","AT2G18800","AT4G23410","AT4G14650","AT3G15353","AT1G14730","AT1G68810","AT1G43790"),fill.by="group",plot_type="box", group.by="clusters",bg.by="lineage",stack=TRUE)
dev.off()