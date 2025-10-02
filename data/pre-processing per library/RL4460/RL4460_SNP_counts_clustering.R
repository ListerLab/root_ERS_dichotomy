library(scDblFinder)
library(dbscan)
library(ggplot2)
library(Seurat)

snps <- read.csv("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4460_cvi_CDFs/genotyping/RL4460_cvi_CDFs_snps.tsv", sep='\t', row.names=1)
sce <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4460_cvi_CDFs/doublet_detection/RL4460_cvi_CDFs_scDblFinder_out.rds")

snps <- snps[colnames(sce)[ sce$scDblFinder.class == "singlet" ],]
snps$log_ref <- log10 (snps$reference_count +1)
snps$log_alt <- log10(snps$alternate_count +1)

db <- dbscan(snps[,c("log_ref","log_alt")], eps=0.2, minPts=50)
snps$cluster <- db$cluster
g <- ggplot(snps, aes(x=log_ref, y=log_alt, colour=factor(cluster)))+geom_point()
ggsave("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4460_cvi_CDFs/genotyping/snps_RL4460_clusters.pdf", plot=g, device="pdf", width=9, height=9)

obj <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4460_cvi_CDFs/seurat_object/RL4460_cvi_CDFs_singlets.rds")
obj <- subset(obj, cells = row.names(subset(snps, cluster == 0)), invert = TRUE)
obj@meta.data[row.names(subset(snps, cluster == 2)), "genotype"] <- "cvi"
obj@meta.data[row.names(subset(snps, cluster == 1)), "genotype"] <- "CDFs"

saveRDS(obj, "/home/moliva/root_sc_paper/data/raw_data_per_library/RL4460_cvi_CDFs/seurat_object/RL4460_cvi_CDFs_singlets_genotype.rds")
