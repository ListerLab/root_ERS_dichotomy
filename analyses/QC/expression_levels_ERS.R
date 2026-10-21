library(Seurat)
library(hdWGCNA)
library(matrixStats)
library(ggplot2)


obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

modules <- GetModules(obj)
genes <- row.names(subset(modules, module %in% c("Module 13","Module 14","Module 20","Module 21")))

Idents(obj) <- "lineage"
obj <- subset(obj, idents = c("xpp"))
obj$state <- gsub(".*(ERS[+-])", "\\1", obj$identity)
Idents(obj) <- "state"

data <- GetAssayData(obj, slot="data")
data <- data[genes,]

avg_ERS_pos <- rowMeans(as.matrix(data[,WhichCells(obj, idents="ERS+")]))
avg_ERS_neg <- rowMeans(as.matrix(data[,WhichCells(obj, idents="ERS-")]))
delta <- (avg_ERS_pos - avg_ERS_neg)

gene_data <- data.frame(delta=delta, avg_ERS_pos=avg_ERS_pos, avg_ERS_neg=avg_ERS_neg, module=modules[genes,"module"])
gene_data[row.names(subset(modules, module %in% c("Module 13","Module 14"))),"type"] <- "symmetric"
gene_data[row.names(subset(modules, module %in% c("Module 20","Module 21"))),"type"] <- "asymmetric"

pdf("/home/moliva/root_sc_paper/analyses/QC/expression_levels_ERS/expression_ERS_pos_delta.pdf", width=8, height=6)
ggplot(gene_data,aes(x=avg_ERS_pos,y=delta, color=type, shape=module))+geom_point()+theme_classic()+geom_smooth(method='lm', aes(group=type))+scale_x_continuous(transform='sqrt')
dev.off()


range_min <- max(
  quantile(subset(gene_data, type=="asymmetric")$avg_ERS_pos, 0.10),
  quantile(subset(gene_data, type=="symmetric")$avg_ERS_pos, 0.10)
)

range_max <- min(
  quantile(subset(gene_data, type=="asymmetric")$avg_ERS_pos, 0.90),
  quantile(subset(gene_data, type=="symmetric")$avg_ERS_pos, 0.90)
)

gene_data_overlap <- subset(
  gene_data,
  avg_ERS_pos >= range_min & avg_ERS_pos <= range_max
)

fit <- lm(avg_ERS_neg ~ avg_ERS_pos,
          data=subset(gene_data, type=="symmetric"))

gene_data$resid_expected <- gene_data$avg_ERS_neg -
  predict(fit, newdata=gene_data)

gene_data_overlap$resid_expected <- gene_data_overlap$avg_ERS_neg -
  predict(fit, newdata=gene_data_overlap)

wtest <- wilcox.test(resid_expected ~ type, data=gene_data_overlap)

tapply(gene_data_overlap$resid_expected, gene_data_overlap$type, median)


pdf("/home/moliva/root_sc_paper/analyses/QC/expression_levels_ERS/expression_ERS_neg_pos.pdf")
ggplot(gene_data, aes(x=avg_ERS_pos, y=avg_ERS_neg, colour=type)) +
  annotate("rect",
           xmin=range_min, xmax=range_max,
           ymin=-Inf, ymax=Inf,
           alpha=0.05, fill="grey70") +
  geom_vline(xintercept=c(range_min, range_max), linetype="dashed", colour="grey40") +
  geom_point(alpha=0.65, size=1.3) +
  geom_abline(slope=1, intercept=0, linetype="dotted", colour="grey50") +
  geom_smooth(
    data=subset(gene_data, type=="symmetric"),
    method="lm", se=FALSE, colour="black", linewidth=0.8
  ) +
  labs(
    x="Average expression in ERS+",
    y="Average expression in ERS-",
    colour=NULL
  ) +
  theme_classic() +
  theme(legend.position="top")
dev.off()


pdf("/home/moliva/root_sc_paper/analyses/QC/expression_levels_ERS/expression_ERS_pos_distribution_symmetry.pdf")
ggplot(gene_data, aes(x=type, y=avg_ERS_pos, fill=type)) +
  annotate("rect",
           xmin=0.5, xmax=2.5,
           ymin=range_min, ymax=range_max,
           alpha=0.05, fill="grey70") +
  geom_hline(yintercept=c(range_min, range_max), linetype="dashed", colour="grey40") +
  geom_violin(trim=FALSE, alpha=0.5, colour=NA) +
  geom_boxplot(width=0.16, outlier.shape=NA, alpha=0.85) +
  geom_jitter(width=0.12, alpha=0.25, size=0.8) +
  labs(
    x=NULL,
    y="Average expression in ERS+"
  ) +
  theme_classic() +
  theme(legend.position="none")
dev.off()



p_label <- if (wtest$p.value < 2.2e-16) {
  "P < 0.0001"
} else {
  paste0("P = ", signif(wtest$p.value, 3))
}

med_asym <- median(subset(gene_data_overlap, type=="asymmetric")$resid_expected)
med_sym  <- median(subset(gene_data_overlap, type=="symmetric")$resid_expected)

med_label <- paste0(
  "median:\nERS+-specific=", round(med_asym, 3),
  "\nShared=", round(med_sym, 3)
)

y_top <- max(gene_data_overlap$resid_expected, na.rm=TRUE)

pdf("/home/moliva/root_sc_paper/analyses/QC/expression_levels_ERS/residuals_distribution_symmetry.pdf")
ggplot(gene_data_overlap, aes(x=type, y=resid_expected, fill=type)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey40") +
  geom_violin(trim=FALSE, alpha=0.5, colour=NA) +
  geom_boxplot(width=0.16, outlier.shape=NA, alpha=0.85) +
  geom_jitter(width=0.12, alpha=0.25, size=0.8) +
  annotate("text", x=1.5, y=y_top, label=p_label, vjust=-0.2, size=3.5) +
  annotate("text", x=1.5, y=y_top - 0.08 * diff(range(gene_data_overlap$resid_expected)),
           label=med_label, vjust=1, size=3.2) +
  labs(
    x=NULL,
    y="Residual ERS- expression\n(observed - expected)"
  ) +
  theme_classic() +
  theme(legend.position="none")
dev.off()




wilcox_results <- data.frame(
  statistic = wtest$statistic,
  p_value = wtest$p.value,
  median_asymmetric = median(subset(gene_data_overlap, type=="asymmetric")$resid_expected),
  median_symmetric  = median(subset(gene_data_overlap, type=="symmetric")$resid_expected)
)
write.csv(wilcox_results, "wilcox_residual_test.csv", row.names = FALSE)