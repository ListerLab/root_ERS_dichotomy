library(Seurat)
library(tibble)
library(purrr)
library(dplyr)
library(tidyr)
library(pheatmap)


col <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")
Idents(col) <- "orig.ident"
Idents(col) <- "orig.ident"
col <- subset(col, idents=c("RL650_col_cvi_r3","RL970_col_ler_r4","RL971_col_c24_r5","RL2757_col_fullMS_r1","RL2758_col_fullMS_r2","RL4207_col_fullMS_r3"))
col$genotype <- "col"

c24 <- readRDS("/home/moliva/root_sc_paper/data/integrated/c24/integrated_c24_predicted_id.rds")
colnames(c24@meta.data)[6] <- "lineage"
colnames(c24@meta.data)[7] <- "identity"
c24$genotype <- "c24"

cvi <- readRDS("/home/moliva/root_sc_paper/data/integrated/cvi/integrated_cvi_predicted_id.rds")
colnames(cvi@meta.data)[6] <- "lineage"
colnames(cvi@meta.data)[7] <- "identity"
cvi$genotype <- "cvi"

ler <- readRDS("/home/moliva/root_sc_paper/data/integrated/ler/integrated_ler_predicted_id.rds")
colnames(ler@meta.data)[5] <- "lineage"
colnames(ler@meta.data)[6] <- "identity"
ler$genotype <- "ler"

ws2 <- readRDS("/home/moliva/root_sc_paper/data/integrated/ws2/integrated_ws2_predicted_id.rds")
colnames(ws2@meta.data)[6] <- "lineage"
colnames(ws2@meta.data)[7] <- "identity"
ws2$genotype <- "ws2"

obj <- merge(col, list(c24,cvi,ler,ws2))
obj <- JoinLayers(obj)
rm(col,c24,cvi,ler,ws2)

pseudo_obj <- AggregateExpression(obj, assays="RNA", return.seurat=T, group.by=c("genotype","identity","orig.ident"))
colnames(pseudo_obj@meta.data) <- "gentoype"
ids <- unique(obj$identity)
ids <- ids[order(ids)]
pseudo_obj$identity <- c(rep(ids,each=3),rep(ids,each=6),rep(rep(ids,each=3),3))
pseudo_obj$geno_id <- paste(pseudo_obj$gentoype, pseudo_obj$identity, sep='_')
Idents(pseudo_obj) <- "geno_id"

comb <- combn(c("c24","col","cvi","ler","ws2"),2)

for(i in 1:ncol(comb)){

  l <- lapply(ids, function(x){
        
	de <- FindMarkers(pseudo_obj, ident.1=paste(comb[1,i], x, sep='_'), ident.2=paste(comb[2,i], x, sep='_'), test.use = "DESeq2")
	de <- subset(de, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
    return(de)
  })
  
  names(l) <- ids
  
  assign(paste0(comb[1,i],"_vs_",comb[2,i]),l)
  
  saveRDS(l,paste0("/home/moliva/root_sc_paper/analyses/ecotype/deg_list_pairwise_per_identity/",paste0(comb[1,i],"_vs_",comb[2,i],".rds")))
  
}



c24_up <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_col[[x]], avg_log2FC >0)), row.names(subset(c24_vs_cvi[[x]], avg_log2FC >0)), row.names(subset(c24_vs_ler[[x]], avg_log2FC >0)), row.names(subset(c24_vs_ws2[[x]], avg_log2FC >0))))
	rm <- unique(c(row.names(col_vs_cvi[[x]]),row.names(col_vs_ler[[x]]),row.names(col_vs_ws2[[x]]),row.names(cvi_vs_ler[[x]]),row.names(cvi_vs_ws2[[x]]),row.names(ler_vs_ws2[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

c24_down <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_col[[x]], avg_log2FC <0)), row.names(subset(c24_vs_cvi[[x]], avg_log2FC <0)), row.names(subset(c24_vs_ler[[x]], avg_log2FC <0)), row.names(subset(c24_vs_ws2[[x]], avg_log2FC <0))))
	rm <- unique(c(row.names(col_vs_cvi[[x]]),row.names(col_vs_ler[[x]]),row.names(col_vs_ws2[[x]]),row.names(cvi_vs_ler[[x]]),row.names(cvi_vs_ws2[[x]]),row.names(ler_vs_ws2[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

col_up <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_col[[x]], avg_log2FC <0)), row.names(subset(col_vs_cvi[[x]], avg_log2FC >0)), row.names(subset(col_vs_ler[[x]], avg_log2FC >0)), row.names(subset(col_vs_ws2[[x]], avg_log2FC >0))))
	rm <- unique(c(row.names(c24_vs_cvi[[x]]),row.names(c24_vs_ler[[x]]),row.names(c24_vs_ws2[[x]]),row.names(cvi_vs_ler[[x]]),row.names(cvi_vs_ws2[[x]]),row.names(ler_vs_ws2[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

col_down <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_col[[x]], avg_log2FC >0)), row.names(subset(col_vs_cvi[[x]], avg_log2FC <0)), row.names(subset(col_vs_ler[[x]], avg_log2FC <0)), row.names(subset(col_vs_ws2[[x]], avg_log2FC <0))))
	rm <- unique(c(row.names(c24_vs_cvi[[x]]),row.names(c24_vs_ler[[x]]),row.names(c24_vs_ws2[[x]]),row.names(cvi_vs_ler[[x]]),row.names(cvi_vs_ws2[[x]]),row.names(ler_vs_ws2[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

cvi_up <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_cvi[[x]], avg_log2FC <0)), row.names(subset(col_vs_cvi[[x]], avg_log2FC <0)), row.names(subset(cvi_vs_ler[[x]], avg_log2FC >0)), row.names(subset(cvi_vs_ws2[[x]], avg_log2FC >0))))
	rm <- unique(c(row.names(c24_vs_col[[x]]),row.names(c24_vs_ler[[x]]),row.names(c24_vs_ws2[[x]]),row.names(col_vs_ler[[x]]),row.names(col_vs_ws2[[x]]),row.names(ler_vs_ws2[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

cvi_down <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_cvi[[x]], avg_log2FC >0)), row.names(subset(col_vs_cvi[[x]], avg_log2FC >0)), row.names(subset(cvi_vs_ler[[x]], avg_log2FC <0)), row.names(subset(cvi_vs_ws2[[x]], avg_log2FC <0))))
	rm <- unique(c(row.names(c24_vs_col[[x]]),row.names(c24_vs_ler[[x]]),row.names(c24_vs_ws2[[x]]),row.names(col_vs_ler[[x]]),row.names(col_vs_ws2[[x]]),row.names(ler_vs_ws2[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

ler_up <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_ler[[x]], avg_log2FC <0)), row.names(subset(col_vs_ler[[x]], avg_log2FC <0)), row.names(subset(cvi_vs_ler[[x]], avg_log2FC <0)), row.names(subset(ler_vs_ws2[[x]], avg_log2FC >0))))
	rm <- unique(c(row.names(c24_vs_col[[x]]),row.names(c24_vs_cvi[[x]]),row.names(c24_vs_ws2[[x]]),row.names(col_vs_cvi[[x]]),row.names(col_vs_ws2[[x]]),row.names(cvi_vs_ws2[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

ler_down <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_ler[[x]], avg_log2FC >0)), row.names(subset(col_vs_ler[[x]], avg_log2FC >0)), row.names(subset(cvi_vs_ler[[x]], avg_log2FC >0)), row.names(subset(ler_vs_ws2[[x]], avg_log2FC <0))))
	rm <- unique(c(row.names(c24_vs_col[[x]]),row.names(c24_vs_cvi[[x]]),row.names(c24_vs_ws2[[x]]),row.names(col_vs_cvi[[x]]),row.names(col_vs_ws2[[x]]),row.names(cvi_vs_ws2[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

ws2_up <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_ws2[[x]], avg_log2FC <0)), row.names(subset(col_vs_ws2[[x]], avg_log2FC <0)), row.names(subset(cvi_vs_ws2[[x]], avg_log2FC <0)), row.names(subset(ler_vs_ws2[[x]], avg_log2FC <0))))
	rm <- unique(c(row.names(c24_vs_col[[x]]),row.names(c24_vs_cvi[[x]]),row.names(c24_vs_ler[[x]]),row.names(col_vs_cvi[[x]]),row.names(col_vs_ler[[x]]),row.names(cvi_vs_ler[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

ws2_down <- lapply(ids, function(x){
	keep <- Reduce(intersect,list(row.names(subset(c24_vs_ws2[[x]], avg_log2FC >0)), row.names(subset(col_vs_ws2[[x]], avg_log2FC >0)), row.names(subset(cvi_vs_ws2[[x]], avg_log2FC >0)), row.names(subset(ler_vs_ws2[[x]], avg_log2FC >0))))
	rm <- unique(c(row.names(c24_vs_col[[x]]),row.names(c24_vs_cvi[[x]]),row.names(c24_vs_ler[[x]]),row.names(col_vs_cvi[[x]]),row.names(col_vs_ler[[x]]),row.names(cvi_vs_ler[[x]])))
	genes <- setdiff(keep,rm)
	return(genes)
})

de_lists <- list( c24_up = c24_up, c24_down = c24_down, col_up = col_up, col_down = col_down, cvi_up = cvi_up, cvi_down = cvi_down, ler_up = ler_up, ler_down = ler_down, ws2_up = ws2_up, ws2_down = ws2_down)

df <- imap_dfr(de_lists, function(de_list, list_name) {
  tibble(
    CellType = ids,
    GeneCount = map_int(de_list, ~ if (is.vector(.x)) length(.x) else nrow(.x)),
    ListName = list_name
  )
})

df <- df %>%
  pivot_wider(names_from = CellType, values_from = GeneCount) %>%
  column_to_rownames("ListName")
  

pdf("/home/moliva/root_sc_paper/analyses/ecotype/ecotype_specific_deg_per_identity.pdf", width=10, height=4)
pheatmap(df, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()


modules <- read.csv("/home/moliva/root_sc_paper/analyses/hdWGCNA/modules.tsv",sep='\t')
modules <- subset(modules, module != 'grey')
gene_modules <- setNames(modules$module, modules$gene_name)

de_by_list <- imap(de_lists, function(de_list, list_name) {
  all_genes <- unique(Reduce(c, de_list))
  genes_with_module <- intersect(all_genes, names(gene_modules))
  module_counts <- table(gene_modules[genes_with_module])
  module_totals <- table(gene_modules)
  prop <- (module_counts / module_totals[names(module_counts)])*100
  as.numeric(prop) %>% setNames(names(prop))
})

prop_df <- bind_rows(de_by_list, .id = "ListName") %>%
  replace(is.na(.), 0)
  
rownames(prop_df) <- prop_df$ListName
prop_df$ListName <- NULL
prop_df <- as.data.frame(prop_df)
prop_df <- prop_df[,paste0("Module ",1:25)]
row.names(prop_df) <- names(de_lists)

pdf("/home/moliva/root_sc_paper/analyses/ecotype/prop_module_ecotype_specific_deg.pdf", width=10, height=4)
pheatmap(prop_df, cluster_rows=FALSE, cluster_cols=FALSE)
dev.off()
