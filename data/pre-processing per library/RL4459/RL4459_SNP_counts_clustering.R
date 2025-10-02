library(scDblFinder)
library(dbscan)
library(ggplot2)
library(Seurat)

snps <- read.csv("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4459_c24_pyl11458/genotyping/RL4459_c24_pyl11458_snps.tsv", sep='\t', row.names=1)
sce <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4459_c24_pyl11458/doublet_detection/RL4459_c24_pyl11458_scDblFinder_out.rds")

snps <- snps[colnames(sce)[ sce$scDblFinder.class == "singlet" ],]
snps$log_ref <- log10 (snps$reference_count +1)
snps$log_alt <- log10(snps$alternate_count +1)

db <- dbscan(snps[,c("log_ref","log_alt")], eps=0.2, minPts=50)
snps$cluster <- db$cluster
g <- ggplot(snps, aes(x=log_ref, y=log_alt, colour=factor(cluster)))+geom_point()
ggsave("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4459_c24_pyl11458/genotyping/snps_RL4459_clusters.pdf", plot=g, device="pdf", width=9, height=9)

obj <- readRDS("/home/moliva/root_sc_paper/data/raw_data_per_library/RL4459_c24_pyl11458/seurat_object/RL4459_c24_pyl11458_singlets.rds")
obj <- subset(obj, cells = row.names(subset(snps, cluster == 0)), invert = TRUE)
obj@meta.data[row.names(subset(snps, cluster == 2)), "genotype"] <- "c24"
obj@meta.data[row.names(subset(snps, cluster == 1)), "genotype"] <- "pyl11458"

saveRDS(obj, "/home/moliva/root_sc_paper/data/raw_data_per_library/RL4459_c24_pyl11458/seurat_object/RL4459_c24_pyl11458_singlets_genotype.rds")
