# Maresin scGEX Analysis.R
# Jacqueline Larouche
# Fall 2022

# Contains the code used for scGEX data analysis of the maresin 1 and vehicle-treated TA VML defects (2mm) at 7dpi. Mice received bilateral 2mm defects to the TA muscle on day 0. Tissues were treated via IM injection of 100ng Maresin 1 in 0.1% EtOH in saline, or 0.1% EtOH in saline at days 1, 3, and 5 post injury. TAs from three mice were pooled, labeled with one CMO per condition (Maresin vs Vehicle), the two conditions were pooled, viable cells were enriched via FACS, and then submitted to the AGC for sequencing.

# Workflow: 
# 1. Initialize variables
# 2. Load ST Data and Create List of Seurat Objects
# 3. Plots for Figure 5 (UMAPs, DotPlot, StackedBar)
# 4. Plots for Supp Figs

#### 0. Initialize Environment ####
library(Seurat) #v3.2.2
library(dplyr) #v1.0.2
library(ggplot2) # v3.2.1
library(tidyverse) # v1.3.0
library(tibble) # v3.0.1
library(Matrix)
library(data.table)
library(RColorBrewer)
library(tibble)
library(cowplot)
library(reshape)
library(dittoSeq)
library(scales) # for hue_pal
library(NCmisc) #v1.1.5 for p.toZ
library(patchwork) # plot_layout

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/TA_VML/Maresin1_Metabolomics/scRNAseq")

#### 1. Initialize Variables ####
sample_dict <- list()
sample_dict[["6629-JL-1_CMO302"]] = "Vehicle"
sample_dict[["6629-JL-1_CMO303"]] = "Maresin"
samples <- names(sample_dict)

Endothelial <- '#5F9C04'
Macrophage <- '#F92B39'
Macrophage_M2 <- '#9C1B64'
Macrophage_Fibrotic <- '#EB2BF9'
Basophil <- '#F9842B'
B_Cell <- '#F9EB2B'
T_conv <- '#5993FF'
Treg <- '#7259FF'
CD8_Tcell <- '#FF59E6'
NK_cell <- '#FF7259'
Granulocyte <- '#BF436E'
Myoblast <- '#00FFAF'
Fibroblast <- '#0050FF'
col_vec <- c(Macrophage_M2, T_conv, NK_cell, B_Cell, Treg, Fibroblast, Endothelial, Granulocyte, Myoblast, Macrophage_Fibrotic, CD8_Tcell, Macrophage, Basophil)
col_vec_reduced <- c(Macrophage_M2, Treg, B_Cell, Fibroblast, Endothelial, Granulocyte, Myoblast)

cells.combined <- readRDS(file = "ProcessedData/maresin_merged_sct_09282022.RDS")

#### 2. Load ST Data and Create List of Seurat Objects ####
veh.data <- Read10X_h5(paste0(getwd(), "/RawData/Sample_6629-JL-1/per_sample_outs/6629-JL-1_CMO302/count/sample_filtered_feature_bc_matrix.h5"))
mar.data <- Read10X_h5(paste0(getwd(), "/RawData/Sample_6629-JL-1/per_sample_outs/6629-JL-1_CMO303/count/sample_filtered_feature_bc_matrix.h5"))
veh <- CreateSeuratObject(counts = veh.data$`Gene Expression`, project = "Vehicle", min.cells = 3, min.features = 200)
mar <- CreateSeuratObject(counts = mar.data$`Gene Expression`, project = "Maresin", min.cells = 3, min.features = 200)

cells.combined <- merge(veh, mar)

cells.combined <- PercentageFeatureSet(cells.combined, pattern = "^mt-", col.name = "percent.mt") 

cells.combined <- cells.combined %>%
  SCTransform(vars.to.regress = "percent.mt") %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30)
cells.combined <- FindClusters(cells.combined, resolution = 0.2)

dittoDimPlot(cells.combined, "seurat_clusters", split.by = c("orig.ident"), color.panel = c(hue_pal()(13)), do.label = TRUE)

Idents(cells.combined) <- 'seurat_clusters'
cluster.degs <- FindAllMarkers(cells.combined, logfc.threshold = 1, only.pos = TRUE)

# Map Clusters to Celltypes
Idents(cells.combined) <- 'seurat_clusters'
current.cluster.ids <- 0:12
celltype.cluster.ids <- c("Macrophage_M2", "T_conv", "NK_cell", "B_cell", "T_reg", "Fibroblast", "Endothelial", "Granulocyte", "Myoblast", "Macrophage_Fibrotic", "CD8_Tcell", "Dendritic", "Basophil")
reduced.cluster.ids <- c("MonocyteDerived", "TNKcell", "TNKcell", "Bcell", "TNKcell", "Fibroblast", "Endothelial", "Granulocyte", "Myoblast", "MonocyteDerived", "TNKcell", "MonocyteDerived", "Granulocyte")
cells.combined$celltype <- plyr::mapvalues(x = cells.combined$seurat_clusters, from = current.cluster.ids, to = celltype.cluster.ids)
cells.combined$celltype_reduced <- plyr::mapvalues(x = cells.combined$seurat_clusters, from = current.cluster.ids, to = reduced.cluster.ids)

saveRDS(cells.combined, file = "ProcessedData/maresin_merged_sct_09282022.RDS")

#### 3. Plots for Figure 5 (UMAPs, DotPlot, StackedBar) ####

pdf(file="Plots/umap_celltype_reduced.pdf", height = 4, width = 6)
DimPlot(cells.combined, group.by = 'celltype_reduced', reduction = "umap", label = FALSE, pt.size = 2, cols = col_vec_reduced)
dev.off()
pdf(file="Plots/umap_celltype.pdf", height = 4, width = 6)
DimPlot(cells.combined, group.by = 'celltype', reduction = "umap", label = FALSE, pt.size = 2, cols = col_vec)
dev.off()

pdf(file="Plots/umap_celltype_reduced_splitby_tx.pdf", height = 2.5, width = 5)
dittoDimPlot(cells.combined, "celltype_reduced", split.by = c("orig.ident"), color.panel = col_vec_reduced)
dev.off()
pdf(file="Plots/umap_celltype_splitby_tx.pdf", height = 2.5, width = 6)
dittoDimPlot(cells.combined, "celltype", split.by = c("orig.ident"), color.panel = col_vec)
dev.off()


Idents(cells.combined) <- 'celltype'
my_levels <- c('Macrophage_M2', 'Dendritic', 'Macrophage_Fibrotic', 'T_conv', 'CD8_Tcell', 'T_reg', 'NK_cell', 'B_cell', 'Fibroblast', 'Endothelial', 'Granulocyte', 'Myoblast', 'Basophil')
Idents(cells.combined) <- factor(Idents(cells.combined), levels= my_levels)

pdf(file="Plots/dotplot_celltype_markers.pdf", height = 5, width = 6)
DotPlot(cells.combined, features = c("Cd68", "Siglech",  "Spp1", "Trem2","Cd3d", "Cd8a", "Ctla4", "Areg", "Gzma", "Ncr1", "Cd79a", "Ms4a1", "Acta2", "Col1a1", "Cdh5", "Fabp4", "S100a8", "Tnnt3", "Myog", "Fcer1a", "Mcpt8")) + coord_flip() + RotatedAxis()
dev.off()
pdf(file="Plots/dotplot_celltype_reduced_markers.pdf", height = 4, width = 5)
DotPlot(cells.combined, group.by = 'celltype_reduced', features = c("Cd68", "Cd209a", "Cd3d", "Cd8a", "Cd4", "Ncr1", "Cd79a", "Ms4a1", "Acta2", "Col1a1", "Cdh5", "Fabp4", "S100a8", "Mcpt8", "Tnnt3", "Myog")) + coord_flip() + RotatedAxis()
dev.off()

pdf("Plots/stackedbar_tx_celltype.pdf", width = 4, height = 4)
dittoBarPlot(cells.combined, group.by = "orig.ident", "celltype",
                   color.panel = c(B_Cell, Basophil, CD8_Tcell, Macrophage,
                                   Endothelial, Fibroblast,
                                   Granulocyte, Macrophage_Fibrotic, 
                                   Macrophage_M2, Myoblast, NK_cell, T_conv,  Treg)) +
  coord_cartesian(ylim=c(0,1)) +
  labs(x="", y="Fraction of Cells", main = "") +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black"),
        axis.title.x = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, colour = "black"))
dev.off()

# Fold change differences in abundance
t1 <- table(cells.combined$orig.ident, cells.combined$celltype_reduced)
t2 <- table(cells.combined$orig.ident)
t3 <- sweep(t1, 1, t2, FUN = '/')
fc <- (t3[1,]-t3[2,])/t3[2,]
fc
df <- data.frame(fc)
df$celltype <- rownames(df)

t4 <- table(cells.combined$orig.ident, cells.combined$celltype)
t5 <- table(cells.combined$orig.ident)
t6 <- sweep(t4, 1, t5, FUN = '/')
fc2 <- (t6[1,]-t6[2,])/t6[2,]
df2 <- data.frame(fc2)
df2$celltype <- rownames(df2)


p1 <- ggplot(df, aes(x = celltype, y = fc, fill = celltype)) + 
  geom_bar(stat = "identity") + 
  labs(x="", y="Fold Change") +
  scale_fill_manual(values = c(B_Cell, Endothelial, Fibroblast, Granulocyte, 
                               Macrophage_M2, Myoblast, Treg)) +
  theme_classic() +
  #scale_x_discrete(limits = c("Vehicle", "Maresin")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,1)) + 
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black"),
        axis.title.x = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
p1
p2 <- ggplot(df2, aes(x = celltype, y = fc2, fill = celltype)) + 
  geom_bar(stat = "identity") + 
  labs(x="", y="Fold Change") +
  scale_fill_manual(values = c(Macrophage, Macrophage_Fibrotic, Macrophage_M2)) +
  theme_classic() +
  scale_x_discrete(limits = c("Dendritic", "Macrophage_Fibrotic", "Macrophage_M2")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.5)) + 
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black"),
        axis.title.x = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
pdf("Plots/barplot_abundance_fold_change.pdf", width = 4, height = 4)
p1
dev.off()
pdf("Plots/barplot_abundance_fold_change_MoDer.pdf", width = 2.5, height = 4)
p2
dev.off()

for (i in 1:ncol(t1)){
  print(colnames(t1)[i])
  res <- prop.test(x = t1[,i], n = t2, correct = FALSE)
  print(res)
}

for (i in c(1, 10, 12)){
  print(colnames(t4)[i])
  res <- prop.test(x = t4[,i], n = t2, correct = FALSE)
  print(res)
}

#### Monocyte Derived Cells Only ####
Idents(cells.combined) <- 'celltype_reduced'
macs <- subset(cells.combined, idents = 'MonocyteDerived')

pdf(file="Plots/umap_macrophages_celltype.pdf", height = 4, width = 6)
DimPlot(macs, group.by = 'celltype', cols = c(Macrophage_M2, Macrophage_Fibrotic, Macrophage), pt.size = 2)
dev.off()

#### 4. Plots for Supp Figs ####
pdf(file="Plots/umap_qc.pdf", height = 4, width = 10)
p1 <- FeaturePlot(cells.combined, features = c("nFeature_RNA"), min.cutoff = "q10", max.cutoff = "q90")
p2 <- FeaturePlot(cells.combined, features = c("nCount_RNA"), min.cutoff = "q10", max.cutoff = "q90")
p1+p2 + plot_layout(ncol = 2)
dev.off()
