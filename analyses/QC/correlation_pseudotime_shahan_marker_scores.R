library(Seurat)
library(tidyr)
library(ggplot2)


obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")

markers <- read.csv("/home/moliva/root_sc_paper/published_datasets/shahan/stage_markers_shahan.csv")

obj$identity <- gsub(" ERS\\+$","", obj$identity)
obj$identity <- gsub(" ERS-$","", obj$identity)
obj$identity <- gsub(" early","", obj$identity)
Idents(obj) <- "identity"


for(i in c("non-hair cells","hair cells","cortex","endodermis","xpp","ppp","procambium","metaphloem","protophloem")){
	sub <- subset(obj, idents=i)
	sub <- AddModuleScore(sub, list(subset(markers, group == paste0("meristem_",i))$gene,subset(markers, group == paste0("elongation_",i))$gene,subset(markers, group == paste0("maturation_",i))$gene), name="score")
	meta <- sub@meta.data[,c("pseudotime","score1","score2","score3")]
	colnames(meta) <- c("pseudotime","meristem","elongation","maturation")
	meta <- meta %>% pivot_longer(cols = c(meristem,elongation,maturation), names_to = "score_type", values_to = "score_value")
	pdf(paste0("/home/moliva/root_sc_paper/analyses/QC/trajectories/correlation_pseudotime_shahan_markers/",i,".pdf"))
	print(ggplot(meta, aes(x = pseudotime, y = score_value, color = score_type)) +
  		geom_point(alpha=0.1) +
  		theme_minimal() +
  		labs(x = "Pseudotime", y = "Score", color = "Score Type") +
  		geom_smooth( aes(group=score_type)) +
  		ggtitle(i))
  	dev.off()
}


for(i in c("metaxylem","protoxylem")){
	sub <- subset(obj, idents=i)
	sub <- AddModuleScore(sub, list(subset(markers, group == paste0("meristem_",i))$gene,subset(markers, group == paste0("elongation_",i))$gene), name="score")
	meta <- sub@meta.data[,c("pseudotime","score1","score2")]
	colnames(meta) <- c("pseudotime","meristem","elongation")
	meta <- meta %>% pivot_longer(cols = c(meristem,elongation), names_to = "score_type", values_to = "score_value")
	pdf(paste0("/home/moliva/root_sc_paper/analyses/QC/trajectories/correlation_pseudotime_shahan_markers/",i,".pdf"))
	print(ggplot(meta, aes(x = pseudotime, y = score_value, color = score_type)) +
  		geom_point(alpha=0.1) +
  		theme_minimal() +
  		labs(x = "Pseudotime", y = "Score", color = "Score Type") +
  		geom_smooth( aes(group=score_type)) +
  		ggtitle(i))
  	dev.off()
}


