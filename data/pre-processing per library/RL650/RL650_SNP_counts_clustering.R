library(scDblFinder)
library(dbscan)
library(ggplot2)
library(Seurat)

snps <- read.csv("RL650_col_cvi_r3_snps.tsv", sep='\t', row.names=1)
sce <- readRDS("../doublet_detection//RL650_col_cvi_r3_scDblFinder_out.rds")

snps <- snps[colnames(sce)[ sce$scDblFinder.class == "singlet" ],]
snps$log_ref <- log10 (snps$reference_count +1)
snps$log_alt <- log10(snps$alternate_count +1)

db <- dbscan(snps[,c("log_ref","log_alt")], eps=0.25, minPts=200)
snps$cluster <- db$cluster
g <- ggplot(snps, aes(x=log_ref, y=log_alt, colour=factor(cluster)))+geom_point()
ggsave("snps_RL650_clusters.pdf", plot=g, device="pdf", width=9, height=9)

obj <- readRDS("../seurat_object/RL650_col_cvi_r3_singlets.rds")
obj_ref <- subset(obj, cells = row.names(subset(snps, cluster == 1)))
obj_alt <- subset(obj, cells = row.names(subset(snps, cluster == 2)))

saveRDS(obj_ref, "../seurat_object/RL650_col_singlets.rds")
saveRDS(obj_alt, "../seurat_object/RL650_cvi_singlets.rds")