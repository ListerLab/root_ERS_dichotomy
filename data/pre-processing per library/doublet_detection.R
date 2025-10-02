suppressMessages(library(scDblFinder))
suppressMessages(library(Seurat))

lib <- commandArgs(trailingOnly = TRUE)[1]
obj <- readRDS(paste0("/home/moliva/root_sc_paper/data/raw_data_per_library/",lib,"/seurat_object/",lib,".rds"))

set.seed(1234)
sce <- scDblFinder(LayerData(obj, assay="RNA", layer="counts"), artificialDoublets=25000, dims=30, dbr.sd=1, propRandom=1)
singlets <- colnames(sce)[sce$scDblFinder.class == "singlet"]
obj_sub <- subset(obj, cells=singlets)

dir.create(paste0("/home/moliva/root_sc_paper/data/raw_data_per_library/",lib,"/doublet_detection/"))
saveRDS(sce, paste0("/home/moliva/root_sc_paper/data/raw_data_per_library/",lib,"/doublet_detection/",lib,"_scDblFinder_out.rds"))
saveRDS(obj_sub, paste0("/home/moliva/root_sc_paper/data/raw_data_per_library/",lib,"/seurat_object/",lib,"_singlets.rds"))