library(GenomicRanges)
library(ChIPpeakAnno)
library(TxDb.Athaliana.BioMart.plantsmart51)
library(rtracklayer)
library(dplyr)


cisbp <- read.csv("/home/moliva/root_sc_paper/data/TFs/cisBP2/TF_Information_all_motifs_plus.txt", sep = '\t')

motifs <- import("/home/moliva/root_sc_paper/data/TFs/dap_data_v4/family_cluster_01_master_chr1-5_GEM_events.narrowPeak.fimo0_1e-4.bed7", format = "BED")
motifs$name <- gsub("\\.","", motifs$name)
motifs$name <- gsub(".*tnt","", motifs$name)
motifs$name <- gsub("_m1","", motifs$name)
motifs$name <- gsub("MYB74_colamp_a","MYB74_col", motifs$name)
motifs$name <- gsub("At5g58900_colamp_a","At5g58900_col_a", motifs$name)

strand(motifs) <- "*"

promoters <- promoters(genes(TxDb.Athaliana.BioMart.plantsmart51), upstream = 500, downstream = 50)
promoters <- trim(promoters)
seqlevelsStyle(promoters) <- "UCSC"

targets <- sapply(unique(motifs$name), function(x){
  m <- subset(motifs, name == x)
  tg <- subsetByOverlaps(promoters,m)

  return(length(unique(tg$gene_id)))
})

targets <- as.data.frame(targets)


cisbp <- merge(cisbp, targets, by.x="DBID.1", by.y="row.names")

cisbp_d <- cisbp %>%
  group_by(DBID) %>%
  filter(TF_Status == "D") %>%
  filter(targets == max(targets))

cisbp_i <- cisbp %>%
  filter(!(DBID %in% cisbp_d$DBID)) %>%
  group_by(DBID) %>%
  filter(targets == max(targets))

tfs <- rbind(cisbp_d,cisbp_i)
tfs <- as.data.frame(tfs)

targets <- lapply(tfs$DBID, function(x){
  m <- subset(tfs, DBID ==x)$DBID.1
  m <- subset(motifs, name == m)
  tg <- subsetByOverlaps(promoters,m)

  return(unique(tg$gene_id))
})

names(targets) <- tfs$DBID


saveRDS(targets,"/home/moliva/root_sc_paper/analyses/TFs/dapseq_targets.rds")