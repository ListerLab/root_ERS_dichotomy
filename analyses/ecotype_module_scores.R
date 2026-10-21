library(Seurat)
library(hdWGCNA)
library(ggplot2)
library(patchwork)



c24 <- readRDS("/home/moliva/root_sc_paper/data/integrated/c24/integrated_c24_predicted_id.rds")
cvi <- readRDS("/home/moliva/root_sc_paper/data/integrated/cvi/integrated_cvi_predicted_id.rds")
ler <- readRDS("/home/moliva/root_sc_paper/data/integrated/ler/integrated_ler_predicted_id.rds")
ws2 <- readRDS("/home/moliva/root_sc_paper/data/integrated/ws2/integrated_ws2_predicted_id.rds")
ref <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

objs <- list(c24=c24,cvi=cvi,ler=ler,ws2=ws2)

objs <- lapply(names(objs), function(nm){

	x <- objs[[nm]]
	x <- ProjectModules(x, seurat_ref= ref, wgcna_name="root_col_modules", wgcna_name_proj="projected")
	x <- ModuleExprScore(x, n_genes="all", method="Seurat")

	plot_list <- ModuleFeaturePlot(x, features="scores", order=TRUE, reduction="umap.rpca")
	pdf(paste0("/home/moliva/root_sc_paper/analyses/ecotype/",nm,"_projected_module_scores.pdf"), width=36, height=30)
	print(wrap_plots(plot_list, ncol=6))
	dev.off()

	scores <- GetModuleScores(x)
	x$Module_20 <- scores$'Module 20'
	x$Module_21 <- scores$'Module 21'

	identity_order <- c("root cap early", "columella","root cap tip","border cells","root cap lateral","epidermis early","hair cells early","hair cells ERS+","hair cells ERS-","non-hair cells","cortex ERS+","cortex ERS-","endodermis ERS+", "endodermis ERS-","ppp ERS+","ppp ERS-","xpp ERS+","xpp ERS-","procambium ERS+","procambium ERS-","metaphloem ERS+","metaphloem ERS-","protophloem ERS+","protophloem ERS-","metaxylem ERS+","metaxylem ERS-","protoxylem early","protoxylem ERS+","protoxylem ERS-")
	Idents(x) <- "predicted.identity"
	Idents(x) <- factor(Idents(x), levels=identity_order)

	pdf(paste0("/home/moliva/root_sc_paper/analyses/ecotype/",nm,"_module20_scores.pdf"), width=10, height=6)
	print(VlnPlot(x,"Module_20"))
	dev.off()
	pdf(paste0("/home/moliva/root_sc_paper/analyses/ecotype/",nm,"_module21_scores.pdf"), width=10, height=6)
	print(VlnPlot(x,"Module_21"))
	dev.off()

})
