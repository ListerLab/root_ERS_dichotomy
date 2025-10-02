library(Seurat)
library(hdWGCNA)
library(tidyverse)
library(cowplot)
library(patchwork)
library(pheatmap)
library(xgboost)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 30)

obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

targets <- readRDS("/home/moliva/root_sc_paper/analyses/TFs/dapseq_targets.rds")

all_genes <- rownames(obj)
row_indices <- integer(0)
col_indices <- integer(0)
gene_to_index <- setNames(seq_along(all_genes), all_genes)

for (tf in names(targets)) {
  targets[[tf]] <- targets[[tf]][targets[[tf]] %in% all_genes]
  row_indices <- c(row_indices, gene_to_index[targets[[tf]]])
  col_indices <- c(col_indices, rep(which(names(targets) == tf), length(targets[[tf]])))
}

sparse_matrix <- Matrix::sparseMatrix(i = row_indices, j = col_indices, dims = c(length(all_genes), length(targets)), dimnames = list(all_genes, names(targets)))

tf_df <- data.frame(motif_name=names(targets), motif_ID=names(targets),gene_name=names(targets))

n_targets <- vector()
for (tf in names(targets)){
	n_targets <- c(n_targets, length(targets[[tf]]))
}

tf_df$n_targets <- n_targets

obj <- SetMotifMatrix(obj, sparse_matrix)
obj <- SetMotifTargets(obj, targets)
obj <- SetMotifs(obj, tf_df)

tf_genes <- names(targets)
modules <- GetModules(obj)
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)
genes_use <- intersect(genes_use,row.names(obj))

obj <-SetDatExpr(obj, features=genes_use)

model_params <- list(objective = 'reg:squarederror', max_depth = 1, eta = 0.1, nthread=30, alpha=0.5)
obj <- ConstructTFNetwork(obj, model_params)
obj <- AssignTFRegulons(obj, strategy="A", reg_thresh=0.01, n_tfs=10)
obj <- RegulonScores(obj, target_type="positive",ncores=30)
obj <- RegulonScores(obj, target_type="negative",ncores=30)

pos_regulon_scores <- GetRegulonScores(obj, target_type="positive")
pos_regulon_scores$identity <- obj@meta.data[row.names(pos_regulon_scores),"identity"]
avg_regulon_scores <- pos_regulon_scores %>% group_by(identity) %>% summarize(across(everything(),mean))
avg_regulon_scores <- as.data.frame(avg_regulon_scores)
row.names(avg_regulon_scores) <- avg_regulon_scores$identity
avg_regulon_scores$identity <- NULL

genenames <- read.csv("/home/moliva/root_sc_paper/genome/gene_descriptions_tair.txt", sep = "\t")
genenames$Primary.Gene.Symbol <- gsub(".*\\((.*)\\).*", "\\1", genenames$Primary.Gene.Symbol)
genenames[genenames$Primary.Gene.Symbol=="","Primary.Gene.Symbol"] <- as.character(genenames[genenames$Primary.Gene.Symbol=="","Locus.Identifier"])
genenames$Primary.Gene.Symbol <- as.character(genenames$Primary.Gene.Symbol)
genenames$Locus.Identifier <- as.character(genenames$Locus.Identifier)
colnames(avg_regulon_scores) <- genenames[match(colnames(avg_regulon_scores), genenames$Locus.Identifier),"Primary.Gene.Symbol"]

identity_order <- c("root cap early", "columella","root cap tip","border cells","root cap lateral","epidermis early","hair cells early","hair cells ERS+","hair cells ERS-","non-hair cells","cortex ERS+","cortex ERS-","endodermis ERS+", "endodermis ERS-","ppp ERS+","ppp ERS-","xpp ERS+","xpp ERS-","procambium ERS+","procambium ERS-","metaphloem ERS+","metaphloem ERS-","protophloem ERS+","protophloem ERS-","metaxylem ERS+","metaxylem ERS-","protoxylem early","protoxylem ERS+","protoxylem ERS-")
pdf("/home/moliva/root_sc_paper/analyses/TFs/heatmap_avg_regulon_scores_2.pdf",width=24, height=6)
pheatmap(avg_regulon_scores[identity_order,], cluster_rows=FALSE,cluster_cols=TRUE, scale="column",cutree_cols=8,clustering_method="ward.D2") 
dev.off()

saveRDS(obj,"/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules_regulons_2.rds")
