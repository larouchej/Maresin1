# Maresin scGEX Analysis.R
# Jacqueline Larouche
# Spring 2023

# Contains the code used for scGEX data analysis of the maresin 1 and vehicle-treated TA VML defects (2mm) at 7dpi. Mice received bilateral 2mm defects to the TA muscle on day 0. Tissues were treated via IM injection of 100ng Maresin 1 in 0.1% EtOH in saline, or 0.1% EtOH in saline at days 1, 3, and 5 post injury. TAs from three mice were pooled for each sample, and each condition was repeated twice (for a total of 6 mice and 12 muscles across two datasets per treatment)


#### 1. Initialize Environment & Variables ####
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
cells.combined <- readRDS(file = "ProcessedData/maresin_merged_integrated_03312023.RDS")

col_vec <- dittoColors(13)
Macrophage <- col_vec[1]
Neutrophil <- col_vec[2]
Tcell <- col_vec[3]
NKcell <- col_vec[4]
MuSC <- col_vec[5]
Myonuclei <- col_vec[6]
Bcell <- col_vec[7]
Endothelial <- col_vec[8]
Dendritic <- col_vec[9]
Fibroblast <- col_vec[10]
Erythrocyte <- col_vec[11]
SchwannCell <- col_vec[12]
Basophil <- col_vec[13]

#### 2. Load GEX Data and Create Seurat Object ####
veh.data <- Read10X_h5(paste0(getwd(), "/RawData/Sample_6629-JL-1/per_sample_outs/6629-JL-1_CMO302/count/sample_filtered_feature_bc_matrix.h5"))
veh2.data <- Read10X_h5(paste0(getwd(), "/RawData/Sample_7698-JL-1/filtered_feature_bc_matrix.h5"))
mar.data <- Read10X_h5(paste0(getwd(), "/RawData/Sample_6629-JL-1/per_sample_outs/6629-JL-1_CMO303/count/sample_filtered_feature_bc_matrix.h5"))
mar2.data <- Read10X_h5(paste0(getwd(), "/RawData/Sample_7698-JL-2/filtered_feature_bc_matrix.h5"))

veh1 <- CreateSeuratObject(counts = veh.data$`Gene Expression`, project = "Vehicle1", min.cells = 3, min.features = 200)
veh2 <- CreateSeuratObject(counts = veh2.data, project = "Vehicle2", min.cells = 3, min.features = 200)
mar1 <- CreateSeuratObject(counts = mar.data$`Gene Expression`, project = "Maresin1", min.cells = 3, min.features = 200)
mar2 <- CreateSeuratObject(counts = mar2.data, project = "Maresin2", min.cells = 3, min.features = 200)
cells.combined <- merge(veh1, c(veh2, mar1, mar2))
table(cells.combined$orig.ident)
# Maresin1 Maresin2 Vehicle1 Vehicle2 
# 2421     8590     1421    10228 
Idents(cells.combined) <- 'orig.ident'
current.cluster.ids <- c('Maresin1', 'Maresin2', 'Vehicle1', 'Vehicle2')
tx.cluster.ids <- c('Maresin', 'Maresin', 'Vehicle', 'Vehicle')
cells.combined$treatment <- plyr::mapvalues(x = cells.combined$orig.ident, from = current.cluster.ids, to = tx.cluster.ids)
table(cells.combined$treatment)

cells.combined <- PercentageFeatureSet(cells.combined, pattern = "^mt-", col.name = "percent.mt") 
VlnPlot(cells.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(cells.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cells.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf(file="Plots/QC_violin.pdf", height = 4, width = 8)
plot1 + plot2
dev.off()

cells.combined <- subset(cells.combined, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & percent.mt < 10)
table(cells.combined$treatment)
# Maresin Vehicle 
# 10462   11217 

#### 3. Integrate using Seurat ####
cells.list <- SplitObject(cells.combined, split.by = "orig.ident")
cells.list <- lapply(X = cells.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = cells.list)
anchors <- FindIntegrationAnchors(object.list = cells.list, anchor.features = features)
cells.combined <- IntegrateData(anchorset = anchors)

DefaultAssay(cells.combined) <- "integrated"
cells.combined <- cells.combined %>% 
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
cells.combined <- FindClusters(cells.combined, resolution = 0.1)

DimPlot(cells.combined, group.by = c('orig.ident', 'seurat_clusters')) 
FeaturePlot(cells.combined, features = c('nCount_RNA', 'nFeature_RNA')) 
table(cells.combined$orig.ident, cells.combined$seurat_clusters)

pdf(file="Plots/UMAP_sample_cluster.pdf", height = 4, width = 10)
DimPlot(cells.combined, group.by = c('seurat_clusters'), label = T) 
dev.off()

saveRDS(cells.combined, "ProcessedData/maresin_merged_integrated_03312023.RDS")

#### 4. Celltype Annotation ####

Idents(cells.combined) <- 'seurat_clusters'
cluster.degs <- FindAllMarkers(cells.combined, logfc.threshold = 1, only.pos = TRUE)
genes <- c("Adgre1", "Siglech", "Clec4d", "Cd3d", "Cd8a", "Cd4", "Foxp3", "Gzma", "Ncr1", "Cd79a", "Ms4a1","Pax7", "Acta2", "Col1a1", "Cdh5", "Fabp4", "Tnnt3", "Myh4", "Pdgfra", "Fcer1a", "Mcpt8", 'Flt3', 'Sox10')
VlnPlot(cells.combined, features = genes, pt.size = 0)

# Map Clusters to Celltypes
Idents(cells.combined) <- 'seurat_clusters'
current.cluster.ids <- 0:16
celltype.cluster.ids <- c("Macrophage", "Neutrophil", "CD8_Tcell", "NKcell", "MuSC", "Myonuclei", "Bcell", "Endothelial", "Treg", "ImmuneProgenitor", "Tcell", "Myeloid", "Dendritic", "Fibroblast", "Erythrocyte", "SchwannCell", "Basophil")
reduced.cluster.ids <- c("Macrophage", "Neutrophil", "Tcell", "NKcell", "MuSC", "Myonuclei", "Bcell", "Endothelial", "Tcell", "Macrophage", "Tcell", "Macrophage", "Dendritic", "Fibroblast", "Erythrocyte", "SchwannCell", "Basophil")
cells.combined$celltype <- plyr::mapvalues(x = cells.combined$seurat_clusters, from = current.cluster.ids, to = celltype.cluster.ids)
cells.combined$celltype_reduced <- plyr::mapvalues(x = cells.combined$seurat_clusters, from = current.cluster.ids, to = reduced.cluster.ids)

saveRDS(cells.combined, file = "ProcessedData/maresin_merged_integrated_03312023.RDS")

#### 5. Plots for main figure (UMAPs, StackedBar, MAST DEGs) ####
pdf(file="Plots/umap_celltype.pdf", height = 4, width = 6)
DimPlot(cells.combined, group.by = 'celltype', reduction = "umap", label = FALSE, pt.size = 2, cols = dittoColors(13))
dev.off()

pdf(file="Plots/umap_celltype_reduced.pdf", height = 4, width = 6)
DimPlot(cells.combined, group.by = 'celltype_reduced', reduction = "umap", label = FALSE, pt.size = 2, cols = dittoColors(13))
dev.off()

pdf(file="Plots/umap_celltype_splitby_tx.pdf", height = 5, width = 2.5)
dittoDimPlot(cells.combined, "celltype", split.by = c("treatment"), split.ncol = 1) + NoLegend()
dev.off()

pdf(file="Plots/umap_celltype_reduced_splitby_tx.pdf", height = 5, width = 2.75)
dittoDimPlot(cells.combined, "celltype_reduced", split.by = c("treatment"), split.ncol = 1)# + NoLegend()
dev.off()

Idents(cells.combined) <- 'celltype'
my_levels <- c("Macrophage", "Dendritic", "Neutrophil", "CD8_Tcell", "Treg", "Tcell", "NKcell", "Bcell", "ImmuneProgenitor", "Myeloid", "Basophil", "Myonuclei", "MuSC", "Fibroblast", "Endothelial",  "SchwannCell", 'Erythrocyte')
Idents(cells.combined) <- factor(Idents(cells.combined), levels= my_levels)
genes <- rev(c("Adgre1", "Cd68", "Siglech", "Clec4d", "Csf3r", "Cd3d", "Cd8a", "Cd4", "Foxp3", "Gzma", "Ncr1", "Cd79a", "Ms4a1", 'Flt3', "Fcer1a", "Mcpt8", "Tnnt3", "Myh4", "Pax7", "Myf5", "Acta2", "Col1a1", "Cdh5", "Fabp4",  'Sox10', 'S100b', 'Hba-a2', 'Hbb-bs'))

# MAST DE by celltype across treatments
cellTx <- paste0(cells.combined$celltype_reduced, "_", cells.combined$treatment)
m <- data.frame("CelltypeReduced_Treatment" = cellTx)
rownames(m) <- rownames(cells.combined@meta.data)
cells.combined <- AddMetaData(object = cells.combined, metadata = m)
Idents(cells.combined) <- "CelltypeReduced_Treatment"

CELLTYPES <- unique(cells.combined$celltype_reduced)

mast_list <- vector(mode="list", length = 17)
names(mast_list) = CELLTYPES

# Compare treatments for each celltype, save each matrix in a list. 
# This step takes a few hours on 1 core
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  # Find  cluster marker genes
  cells.mast.de <- FindMarkers(object = cells.combined,
                               ident.1 = paste0(CELLTYPE, "_Maresin"),
                               ident.2 = paste0(CELLTYPE, "_Vehicle"),
                               only.pos = FALSE, 
                               min.pct = 0,
                               logfc.threshold = 0,
                               test.use = "MAST")
  cells.mast.de <- as.data.frame(cells.mast.de)
  cells.mast.de <- rownames_to_column(cells.mast.de, var = "gene")
  cells.mast.de$celltype <- rep(CELLTYPE, dim(cells.mast.de)[1])
  mast_list[[CELLTYPE]] <- cells.mast.de
}

# Combine all matrices into one dataframe
mast_df <- data.frame()
for (CELLTYPE in CELLTYPES) {
  print(CELLTYPE)
  mast_df  <- rbind(mast_df, mast_list[[CELLTYPE]])
}
dim(mast_df); head(mast_df) # 34000 by 7

# Save differential expression results
write.csv(mast_df, "Plots/mast_df_reduced.csv")

df <- read_csv("Plots/mast_df_reduced.csv")
df$z <- p.to.Z(df$p_val) * sign(df$avg_log2FC)
df$z.adj <- p.to.Z(df$p_val_adj) * sign(df$avg_log2FC)
df <- df[sample(nrow(df)), ]
df$celltype_factor <- factor(df$celltype,  levels=my_levels, ordered=T)

# Make custom color column to facilitate grey coloring by threshold.
mast_col <- c(Macrophage, Dendritic, Neutrophil, Tcell, NKcell, Bcell, Basophil, Myonuclei, MuSC, Fibroblast, Endothelial, SchwannCell, Erythrocyte)
col <- mast_col[df$celltype_factor]
col[df$p_val_adj > 0.05] <- "#D3D3D3" # grey
df$col <- as.factor(col)

q <- ggplot(df, aes(x = celltype_factor, y = z.adj, color = col)) +
  geom_jitter(width = 0.40, alpha = .55, size = 1) +
  coord_cartesian(ylim=c(-50,50)) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size = 15, color = "black"), 
        axis.title.x = element_blank(), panel.border = element_blank(), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), panel.background = element_blank()) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(legend.position="none") +
  labs(y = "Z-score") +
  scale_color_manual(values = levels(df$col)) +
  geom_hline(aes(yintercept=0), color="darkgrey", linetype="dashed")
q
pdf("Plots/mast_degs_reduced.pdf", width = 7, height = 5)
q
dev.off()

# Fold change differences in abundance
t1 <- table(cells.combined$treatment, cells.combined$celltype_reduced)
t2 <- table(cells.combined$treatment)
t3 <- sweep(t1, 1, t2, FUN = '/')
fc <- (t3[1,]-t3[2,])/t3[2,]
fc
df <- data.frame(fc)
df$celltype <- rownames(df)

p1 <- ggplot(df, aes(x = celltype, y = fc, fill = celltype)) + 
  geom_bar(stat = "identity") + 
  labs(x="", y="Fold Change") +
  scale_fill_manual(values = c(Basophil, Bcell, Dendritic, Endothelial, Erythrocyte, Fibroblast, Macrophage, MuSC, Myonuclei, Neutrophil, NKcell, SchwannCell, Tcell)) +
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

pdf("Plots/barplot_abundance_fold_change.pdf", width = 6, height = 5)
p1
dev.off()

for (i in 1:ncol(t1)){
  print(colnames(t1)[i])
  res <- prop.test(x = t1[,i], n = t2, correct = FALSE)
  print(res)
}


#### 6. Plots for Macrophages Only ####
Idents(cells.combined) <- 'celltype_reduced'
macs <- subset(cells.combined, idents = 'Macrophage')

# Re-cluster
macs <- macs %>% ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30)
macs <- FindClusters(macs, resolution = 0.1)
VlnPlot(macs, features = c('Adgre1', 'Cd68'), pt.size = 0)
macs <- subset(macs, idents = '3', invert = T)

Idents(macs) <- 'treatment'
deg.macs.tx <- FindMarkers(macs, ident.1 = 'Maresin', logfc.threshold = 0, test.use = 'MAST')
deg.macs.tx$gene <- rownames(deg.macs.tx)


goi <- c('Cxcl2', 'Cxcl9', 'Ccl2', 'Chil3', 'Ly6c2', 'Tgfbi', 'Hmox1', 'Thbs1', 'Csf3r', 
         'Mmp9', 'Cd83', 'Egr1', 'Igf1', 'Fscn1', 'Marcks', 'Uap1', 'Naaa')

library(EnhancedVolcano)
vol <- EnhancedVolcano(deg.macs.tx,
                         lab = rownames(deg.macs.tx),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         title = '',
                         subtitle = '',
                         xlab = bquote(~Log[2]~ ("fold change")),
                         ylab = bquote(~-Log[10]("p-adj")),
                         axisLabSize = 25,
                         caption = NULL,
                         selectLab = goi,
                         xlim = c(-.5, .5),
                         pCutoff = 0.05,
                         FCcutoff = 0.0585,
                         pointSize = 1.0,
                         labSize = 6.0,
                         boxedLabels = FALSE,
                         colAlpha = 0.5,
                         legendPosition = 'none',
                         drawConnectors = TRUE,
                         widthConnectors = 0.2,
                         colConnectors = 'black',
                         col = c('#999999', '#009E73', '#56B4E9', '#E69F00'),
                         max.overlaps = 25)
vol

pdf("Plots/volcano_macs_pbs_v_mar1_v2.pdf", width = 7, height = 7)
vol
dev.off()

#### 7. Supp Figs ####
pdf(file="Plots/umap_qc.pdf", height = 4, width = 10)
p1 <- FeaturePlot(cells.combined, features = c("nFeature_RNA"), min.cutoff = "q10", max.cutoff = "q90")
p2 <- FeaturePlot(cells.combined, features = c("nCount_RNA"), min.cutoff = "q10", max.cutoff = "q90")
p1+p2 + plot_layout(ncol = 2)
dev.off()

Idents(cells.combined) <- 'celltype_reduced'
my_levels <- c("Macrophage", "Dendritic", "Neutrophil", "Tcell",  "NKcell", "Bcell", "Basophil", "Myonuclei", "MuSC", "Fibroblast", "Endothelial",  "SchwannCell", 'Erythrocyte')
Idents(cells.combined) <- factor(Idents(cells.combined), levels= my_levels)
genes <- rev(c("Adgre1", "Cd68", "Siglech", "Clec4d", "Csf3r", "Cd3d", "Cd8a", "Cd4", "Foxp3", "Gzma", "Ncr1", "Cd79a", "Ms4a1", "Fcer1a", "Mcpt8", "Tnnt3", "Myh4", "Pax7", "Myf5", "Acta2", "Col1a1", "Cdh5", "Fabp4",  'Sox10', 'S100b', 'Hba-a2', 'Hbb-bs'))
pdf(file="Plots/dotplot_celltype_reduced_markers.pdf", height = 7, width = 9)
DotPlot(cells.combined, features = genes, assay = 'RNA') + coord_flip() + RotatedAxis()
dev.off()
