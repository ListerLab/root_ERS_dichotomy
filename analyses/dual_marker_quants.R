library(mclust)
library(ggplot2)
library(Seurat)


quant <- read.csv("/home/moliva/root_sc_paper/analyses/segmentation/Per_cell_quant_AT3G01720_AT5G56550_after_substract.csv",header=1)
quant$ratio <- quant$meanC2/(quant$meanC2 + quant$meanC3)

fit <- Mclust(quant$ratio)

quant$cluster <- fit$classification

pdf("/home/moliva/root_sc_paper/analyses/segmentation/scatter_quant_AT3G01720_AT5G56550.pdf", width=7, height=6)
ggplot(quant, aes(x = meanC2, y=meanC3, color=factor(cluster))) +
  geom_point() +
  labs(x = "quantification AT3G01720",
       y = "quantification AT5G56550") +
  theme_classic()
dev.off() 

means <- fit$parameters$mean
sds   <- sqrt(fit$parameters$variance$sigmasq)
props <- fit$parameters$pro


pdf("/home/moliva/root_sc_paper/analyses/segmentation/AT3G01720_AT5G56550_quant_ratio_histogram.pdf")
ggplot(quant, aes(x = ratio)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 20,
                 fill = "grey85",
                 color = "white") +
  stat_function(fun = function(x)
    props[1] * dnorm(x, means[1], sds[1]),
    color = "blue", linewidth = 1) +
  stat_function(fun = function(x)
    props[2] * dnorm(x, means[2], sds[2]),
    color = "red", linewidth = 1) +
  labs(x = "Marker ratio",
       y = "Density") +
  theme_classic()
dev.off() 
  
  
obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories.rds")

ids <- grep("hair cells",unique(obj$identity),value=TRUE)

Idents(obj) <- "identity"

cells <- WhichCells(obj, idents=ids)

df <-data.frame(C2=GetAssayData(obj, slot = "data")["AT3G01720",cells],C3=GetAssayData(obj, slot = "data")["AT5G56550",cells])
df <- df[df$C2 + df$C3 >0,]
df$score <- df$C2-df$C3

fit_sc <- Mclust(df$score)

means <- fit_sc$parameters$mean
sds   <- sqrt(fit_sc$parameters$variance$sigmasq)
props <- fit_sc$parameters$pro


df$identity <- obj@meta.data[row.names(df),"identity"]
df$group <- "H"
df[grep("non-hair",df$identity),"group"] <- "N"
df$group <- factor(df$group, levels = c("H", "N"))

pdf("/home/moliva/root_sc_paper/analyses/segmentation/AT3G01720_AT5G56550_sc_score_histogram_per_epi_identity.pdf")
ggplot(df, aes(x = score, fill = group)) +
  geom_histogram(aes(y = after_stat(count / sum(count) / width)),
                 bins = 20,
                 colour = "white",
  				 linewidth = 0.2,
                 position = "stack") +
  stat_function(fun = function(x)
    props[1] * dnorm(x, means[1], sds[1]),
    color = "blue", linewidth = 1) +
  stat_function(fun = function(x)
    props[2] * dnorm(x, means[2], sds[2]),
    color = "red", linewidth = 1)+
  scale_fill_manual(values = c(H = "grey80",N = "grey50")) +
  theme_classic()
dev.off()