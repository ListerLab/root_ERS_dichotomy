library(dplyr)
library(ggplot2)

proto_1 <- read.csv('/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/all_genes_rep1.tsv',sep='\t')
proto_2 <- read.csv('/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/all_genes_rep2.tsv',sep='\t')

modules <- read.csv("/home/moliva/root_sc_paper/analyses/hdWGCNA/modules.tsv",sep='\t')

proto_1 <- subset(proto_1, padj < 0.05)
proto_2 <- subset(proto_2, padj < 0.05)
proto_1$AGI <- row.names(proto_1)
proto_2$AGI <- row.names(proto_2)

proto_de <- merge(proto_1,proto_2,all=FALSE, suffixes=c("_1","_2"), by="AGI")
proto_de$log2FC <- (proto_de$log2FoldChange_1 + proto_de$log2FoldChange_2)/2

colnames(modules)[1] <- "AGI"

proto_de <- merge(proto_de,modules,all=FALSE, by="AGI")

pdf("/home/moliva/root_sc_paper/analyses/QC/protoplasting/module_genes_de_protoplasting.pdf", width=8, height=5)
proto_de %>% filter(module != "grey") %>%
  ggplot(aes(x=factor(module, levels = paste0("Module ",1:25)), y=log2FC, color=module)) +
  theme(panel.grid.major = element_line(colour = ),panel.background = element_rect(fill = NA, colour = "black"), axis.text.x = element_text(angle = 90,vjust=0.4, hjust = 1),axis.text.y = element_text(vjust = 0.4, hjust = 1), legend.key = element_rect(fill=NA), legend.position = "none")+
  xlab(label = "gene module")+
  geom_boxplot(color="black",outlier.alpha = 0) +
  geom_point(size=2, alpha=0.5)+
  ylim(-2.5,2.5)
dev.off()