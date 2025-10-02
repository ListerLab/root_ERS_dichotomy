library(Seurat)

obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

pdf("/home/moliva/root_sc_paper/analyses/ERS_features/umap_col_ngenes.pdf", width=6, height=6)
FeaturePlot(obj,"nFeature_RNA", min.cutoff="q5", max.cutoff="q95")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/ERS_features/umap_col_nUMI.pdf", width=6, height=6)
FeaturePlot(obj,"nCount_RNA", min.cutoff="q10", max.cutoff="q90")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/ERS_features/AT1G10140.pdf", width=6, height=6)
FeaturePlot(obj,"AT1G10140", slot="counts", min.cutoff="q10", max.cutoff="q90")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/ERS_features/AT4G26080.pdf", width=6, height=6)
FeaturePlot(obj,"AT4G26080", min.cutoff="q10", max.cutoff="q90")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/ERS_features/AT5G56550.pdf", width=6, height=6)
FeaturePlot(obj,"AT5G56550", min.cutoff="q10", max.cutoff="q90")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/ERS_features/AT2G41660.pdf", width=6, height=6)
FeaturePlot(obj,"AT2G41660", min.cutoff="q10", max.cutoff="q90")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/ERS_features/AT2G16720.pdf", width=6, height=6)
FeaturePlot(obj,"AT2G16720", min.cutoff="q10", max.cutoff="q90")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/ERS_features/AT1G43160.pdf", width=6, height=6)
FeaturePlot(obj,"AT1G43160", min.cutoff="q10", max.cutoff="q90")
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/ERS_features/AT3G04010.pdf", width=6, height=6)
FeaturePlot(obj,"AT3G04010", min.cutoff="q10", max.cutoff="q90")
dev.off()


