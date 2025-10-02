library(scDblFinder)
library(dbscan)
library(ggplot2)
library(Seurat)

snps <- read.csv("/home/moliva/root_sc_paper/data/raw_data_per_library/RL1465_c24_N681600/genotyping/RL1465_c24_N681600_snps.tsv", sep='\t', row.names=1)
sce <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL1465_c24_N681600/doublet_detection/RL1465_c24_N681600_scDblFinder_out.rds")

snps <- snps[colnames(sce)[ sce$scDblFinder.class == "singlet" ],]
snps$log_ref <- log10 (snps$reference_count +1)
snps$log_alt <- log10(snps$alternate_count +1)

db <- dbscan(snps[,c("log_ref","log_alt")], eps=0.2, minPts=80)
snps$cluster <- db$cluster
g <- ggplot(snps, aes(x=log_ref, y=log_alt, colour=factor(cluster)))+geom_point()
ggsave("/home/moliva/root_sc_paper/data/raw_data_per_library/RL1465_c24_N681600/genotyping/snps_RL1465_clusters.pdf", plot=g, device="pdf", width=9, height=9)

obj <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL1465_c24_N681600/seurat_object/RL1465_c24_N681600_singlets.rds")
obj <- subset(obj, cells = row.names(subset(snps, cluster == 0)), invert = TRUE)
obj@meta.data[row.names(subset(snps, cluster == 1)), "genotype"] <- "c24"
obj@meta.data[row.names(subset(snps, cluster == 2)), "genotype"] <- "N681600"

saveRDS(obj, "/home/moliva/root_sc_paper/data/raw_data_per_library/RL1465_c24_N681600/seurat_object/RL1465_c24_N681600_singlets_genotype.rds")
