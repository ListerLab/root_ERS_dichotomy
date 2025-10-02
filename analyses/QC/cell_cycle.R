library(Seurat)
library(dplyr)
library(ggplot2)
options(future.globals.maxSize = 2500 * 1024^2)

obj <- readRDS("/home/moliva/root_sc_paper/data/integrated/col/integrated_col_trajectories_modules.rds")

cell_cycle <- read.csv("/home/moliva/root_sc_paper/data/cell_cycle/cell_cycle_genes.csv")

obj <- CellCycleScoring(obj, s.features = subset(cell_cycle, phase == "S_genes")$AGI, g2m.features = subset(cell_cycle, phase == "G2_M_genes")$AGI, set.ident = TRUE)


meta <- obj@meta.data %>%
  dplyr::select(identity, Phase) %>%
  group_by(identity, Phase) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))
  
identity_order <- c("root cap early", "columella","root cap tip","border cells","root cap lateral","epidermis early","hair cells early","hair cells ERS+","hair cells ERS-","non-hair cells","cortex ERS+","cortex ERS-","endodermis ERS+", "endodermis ERS-","ppp ERS+","ppp ERS-","xpp ERS+","xpp ERS-","procambium ERS+","procambium ERS-","metaphloem ERS+","metaphloem ERS-","protophloem ERS+","protophloem ERS-","metaxylem ERS+","metaxylem ERS-","protoxylem early","protoxylem ERS+","protoxylem ERS-")


pdf("/home/moliva/root_sc_paper/analyses/QC/cell_cycle/proportion_cell_cycle_phase_identity.pdf", width=8,height=6)
ggplot(meta, aes(x=factor(identity,levels=identity_order), y = perc*100, fill=Phase)) + geom_bar(stat="identity")+ labs(x="identity",y= "proportion (%)", fill="cell cycle phase") + theme(panel.grid.major = element_line(colour = NA),panel.background = element_rect(fill = "white", colour = "black"),axis.text.x = element_text(angle = 90))
dev.off()

pdf("/home/moliva/root_sc_paper/analyses/QC/cell_cycle/umap_cell_cycle_phase.pdf", width=7, height=6)
DimPlot(obj, group.by="Phase")
dev.off()

obj_reg <- obj
obj_reg[["RNA"]] <- split(obj_reg[["RNA"]], f = obj_reg$orig.ident)
obj_reg <- ScaleData(obj, vars.to.regress=c("S.Score", "G2M.Score"))

obj_reg <- FindVariableFeatures(obj_reg, nfeatures=25000,verbose=F)
proto <- read.csv("/home/moliva/root_sc_paper/data/protoplasting_sensitive_genes/de_genes.txt", header=FALSE)
mito <- grep("ATM", row.names(obj_reg), value=TRUE)
chloro <- grep("ATC", row.names(obj_reg), value=TRUE)
VariableFeatures(obj_reg) <- setdiff(VariableFeatures(obj_reg),c(proto$V1,mito,chloro))

obj_reg <- RunPCA(obj_reg, npcs=50,verbose=F)
obj_reg <- IntegrateLayers(obj_reg, method= RPCAIntegration, orig.reduction="pca", new.reduction="integrated.rpca", dims=1:50, reference=1:3, verbose=FALSE)
obj_reg <- RunUMAP(obj_reg, reduction="integrated.rpca", dims=1:45, reduction.name="umap.reg",verbose=F)

pdf("/home/moliva/root_sc_paper/analyses/QC/cell_cycle/umap_integrated_after_cellcycle_regression.pdf", width=9, height=6)
DimPlot(obj_reg, group.by="identity",reduction="umap.reg")
dev.off()
