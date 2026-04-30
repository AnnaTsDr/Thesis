library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(monocle3)
library(clustree)
library(SeuratWrappers)
library(base)
library(scater)
library(cowplot)
library(patchwork)
library(dbplyr)
library(data.table)
library(EnhancedVolcano)
library(SeuratDisk)
library(ggrepel)
library(clustifyr)
library(clustifyrdata)

setwd("E:/Immunaging")
setwd("F:/Immunaging")
memory.limit(1024000)#memory in MB

library(MazamaCoreUtils)
library(pbapply)
library(prob)
source(file = "D:/Anna/HezisData/AssessNodesAdaptedfromSeurat.R")
source(file = "AssessNodesAdaptedfromSeurat.R")

#GSE152522 - Imbalance of Regulatory and Cytotoxic SARS-CoV-2-Reactive CD4+ T Cells in COVID-19
#cd4t0n6
cd4t0n6.raw <- read.csv(file = "GSE152522/cd4t0n6_umi.csv", sep = "\t", row.names = 1)
#cd4t0n6.raw <- read.delim(file = "GSE152522/GSE152522_cd4t0n6_umi.txt.gz", row.names = 1)
object.size(cd4t0n6.raw)
dim(cd4t0n6.raw)

cd4t0n6 <- CreateSeuratObject(counts = cd4t0n6.raw, project = "cd4t0n6", min.cells = 3, min.features = 200)
object.size(cd4t0n6)

saveRDS(cd4t0n6, file = "cd4t0n6_seurat_object.rds")
cd4t0n6 <- readRDS(file = "cd4t0n6_seurat_object.rds")

Metadata.cd4t0n6 <- read.delim(file = "GSE152522/GSE152522_cd4t0n6_annotation.txt.gz")
#Metadata.cd4t0n6 <- read.table(file = "GSE152522/GSE152522_cd4t0n6_annotation.txt.gz")

cd4t0n6 <- AddMetaData(cd4t0n6, metadata = Metadata.cd4t0n6)

saveRDS(cd4t0n6, file = "cd4t0n6.rds")

#cd4t6
cd4t6.raw <- read.csv(file = "GSE152522/cd4t6_umi.csv", sep = "\t", row.names = 1)
cd4t6.raw <- read.table(file = "GSE152522/GSE152522_cd4t6_umi.txt.gz", row.names = 1)
object.size(cd4t6.raw)
dim(cd4t6.raw)

cd4t6 <- CreateSeuratObject(counts = cd4t6.raw, project = "cd4t6", min.cells = 3, min.features = 200)
object.size(cd4t0n6)

saveRDS(cd4t6, file = "cd4t6_seurat_object.rds")
cd4t6 <- readRDS(file = "cd4t6_seurat_object.rds")

Metadata.cd4t6 <- read.delim(file = "GSE152522/GSE152522_cd4t6_annotation.txt.gz")
#Metadata.cd4t6 <- read.table(file = "GSE152522_cd4t6_annotation.txt.gz", row.names = 1)
cd4t6 <- AddMetaData(cd4t6, metadata = Metadata.cd4t6)

#cd4t24
cd4t24.raw <- read.csv(file = "GSE152522/cd4t24_umi.csv", sep = "\t", row.names = 1)
cd4t24.raw <- read.table(file = "GSE152522/GSE152522_cd4t24_umi.txt.gz", row.names = 1)
object.size(cd4t24.raw)
dim(cd4t24.raw)

cd4t24 <- CreateSeuratObject(counts = cd4t24.raw, project = "cd4t24", min.cells = 3, min.features = 200)
object.size(cd4t24)

saveRDS(cd4t24, file = "cd4t24_seurat_object.rds")
cd4t24 <- readRDS(file = "cd4t24_seurat_object.rds")

Metadata.cd4t24 <- read.delim(file = "GSE152522/GSE152522_cd4t24_annotation.txt.gz")
cd4t24 <- AddMetaData(cd4t24, metadata = Metadata.cd4t24)

Covid_data <- merge(cd4t0n6, c(cd4t6, cd4t24))

saveRDS(Covid_data, file = "Covid_data.rds")
Covid_data <- readRDS(file = "Covid_data.rds")

#Hashimoto
genes <- read.csv(file = "Centenarians/genes.csv", row.names = 1)
Supercentenerians.raw <- read.table(file = "Centenarians/01.UMI.txt.gz")
row.names(Supercentenerians.raw) <- make.names(genes[,2], unique = TRUE)
Supercentenerians <- CreateSeuratObject(counts = Supercentenerians.raw, project = "Supercentenerians", min.cells = 3, min.features = 200)
object.size(Supercentenerians.raw)
dim(Supercentenerians.raw)

write.csv(rownames(Supercentenerians), file = "genes_supercentenerians.csv")
genes <- read.csv(file = "genes_supercentenerians.csv", row.names = 2)
row.names(Supercentenerians@assays$RNA) <- row.names(genes)

saveRDS(Supercentenerians, file = "Supercentenerians_seurat_object.rds")
Supercentenerians <- readRDS(file = "Supercentenerians_seurat_object.rds")

Metadata.Supercentenerians <- read.table(file = "Centenarians/03.Cell.Barcodes.txt", row.names = 1)

#GSE149652 - Intratumoral CD4+ T Cells Mediate Anti-tumor Cytotoxicity in Human Bladder Cancer
#Anti-PD-L1_A_tumor_CD4_droplet_count_matrice
Anti_PD_L1_A_tumor.raw <- read.csv(file = "GSE149652/Anti-PD-L1_A_tumor_CD4_droplet_count_matrice.csv", sep = ",", row.names = 1)
object.size(Anti_PD_L1_A_tumor.raw)
dim(Anti_PD_L1_A_tumor.raw)

Anti_PD_L1_A_tumor <- CreateSeuratObject(counts = Anti_PD_L1_A_tumor.raw, project = "Anti_PD_L1_A_tumor", min.cells = 3, min.features = 200)
object.size(Anti_PD_L1_A_tumor)

saveRDS(Anti_PD_L1_A_tumor, file = "Anti_PD_L1_A_tumor_seurat_object.rds")
Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_seurat_object.rds")

#Anti-PD-L1_B_tumor_CD4_droplet_count_matrice
Anti_PD_L1_B_tumor.raw <- read.csv(file = "GSE149652/Anti-PD-L1_B_tumor_CD4_droplet_count_matrice.csv", sep = ",", row.names = 1)
object.size(Anti_PD_L1_B_tumor.raw)
dim(Anti_PD_L1_B_tumor.raw)

Anti_PD_L1_B_tumor <- CreateSeuratObject(counts = Anti_PD_L1_B_tumor.raw, project = "Anti_PD_L1_B_tumor", min.cells = 3, min.features = 200)
object.size(Anti_PD_L1_B_tumor)

saveRDS(Anti_PD_L1_B_tumor, file = "Anti_PD_L1_B_tumor_seurat_object.rds")
Anti_PD_L1_B_tumor <- readRDS(file = "Anti_PD_L1_B_tumor_seurat_object.rds")

#Anti-PD-L1_C_tumor_CD4_droplet_count_matrice
Anti_PD_L1_C_tumor.raw <- read.csv(file = "GSE149652/Anti-PD-L1_C_tumor_CD4_droplet_count_matrice.csv", sep = ",", row.names = 1)
object.size(Anti_PD_L1_C_tumor.raw)
dim(Anti_PD_L1_C_tumor.raw)

Anti_PD_L1_C_tumor <- CreateSeuratObject(counts = Anti_PD_L1_C_tumor.raw, project = "Anti_PD_L1_C_tumor", min.cells = 3, min.features = 200)
object.size(Anti_PD_L1_C_tumor)

saveRDS(Anti_PD_L1_C_tumor, file = "Anti_PD_L1_C_tumor_seurat_object.rds")
Anti_PD_L1_C_tumor <- readRDS(file = "Anti_PD_L1_C_tumor_seurat_object.rds")

#Anti-PD-L1_D_tumor_CD4_droplet_count_matrice
Anti_PD_L1_D_tumor.raw <- read.csv(file = "GSE149652/Anti-PD-L1_D_tumor_CD4_droplet_count_matrice.csv", sep = ",", row.names = 1)
object.size(Anti_PD_L1_D_tumor.raw)
dim(Anti_PD_L1_D_tumor.raw)

Anti_PD_L1_D_tumor <- CreateSeuratObject(counts = Anti_PD_L1_D_tumor.raw, project = "Anti_PD_L1_D_tumor", min.cells = 3, min.features = 200)
object.size(Anti_PD_L1_D_tumor)

saveRDS(Anti_PD_L1_D_tumor, file = "Anti_PD_L1_D_tumor_seurat_object.rds")
Anti_PD_L1_D_tumor <- readRDS(file = "Anti_PD_L1_D_tumor_seurat_object.rds")

#Chemo_tumor_CD4_droplet_count_matrice
Chemo_tumor.raw <- read.csv(file = "GSE149652/Chemo_tumor_CD4_droplet_count_matrice.csv", sep = ",", row.names = 1)
object.size(Chemo_tumor.raw)
dim(Chemo_tumor.raw)

Chemo_tumor <- CreateSeuratObject(counts = Chemo_tumor.raw, project = "Chemo_tumor", min.cells = 3, min.features = 200)
object.size(Chemo_tumor)

saveRDS(Chemo_tumor, file = "Chemo_tumor_seurat_object.rds")
Chemo_tumor <- readRDS(file = "Chemo_tumor_seurat_object.rds")

#Untreated_A_tumor_CD4_droplet_count_matrice
Untreated_A_tumor.raw <- read.csv(file = "GSE149652/Untreated_A_tumor_CD4_droplet_count_matrice.csv", sep = ",", row.names = 1)
object.size(Untreated_A_tumor.raw)
dim(Untreated_A_tumor.raw)

Untreated_A_tumor <- CreateSeuratObject(counts = Untreated_A_tumor.raw, project = "Untreated_A_tumor", min.cells = 3, min.features = 200)
object.size(Untreated_A_tumor)

saveRDS(Untreated_A_tumor, file = "Untreated_A_tumor_seurat_object.rds")
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_seurat_object.rds")

#Chemo_tumor_CD4_droplet_count_matrice
Untreated_B_tumor.raw <- read.csv(file = "GSE149652/Untreated_B_tumor_CD4_droplet_count_matrice.csv", sep = ",", row.names = 1)
object.size(Untreated_B_tumor.raw)
dim(Untreated_B_tumor.raw)

Untreated_B_tumor <- CreateSeuratObject(counts = Untreated_B_tumor.raw, project = "Untreated_B_tumor", min.cells = 3, min.features = 200)
object.size(Untreated_B_tumor)

saveRDS(Untreated_B_tumor, file = "Untreated_B_tumor_seurat_object.rds")
Untreated_B_tumor <- readRDS(file = "Untreated_B_tumor_seurat_object.rds")

tumor_data <- merge(Anti_PD_L1_A_tumor, c(Anti_PD_L1_B_tumor, Anti_PD_L1_C_tumor, Anti_PD_L1_D_tumor, Untreated_A_tumor, Untreated_B_tumor, Chemo_tumor))

saveRDS(tumor_data, file = "tumor_data.rds")
tumor_data <- readRDS(file = "tumor_data.rds")

#merge data
cd4_big_data <- merge(Supercentenerians, c(tumor_data, Covid_data))

#cd4_big_data <- CreateSeuratObject(counts = cd4_big_data@assays$RNA@counts, min.cells = 3, min.features = 200)

saveRDS(cd4_big_data, file = "cd4_big_data_seurat_object.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_seurat_object.rds")

cd4_big_data$condition <- plyr::mapvalues(x = cd4_big_data$orig.ident, from = c("cd4t0n6", "cd4t6", "cd4t24", "Anti_PD_L1_A_tumor",
                                                                                "Anti_PD_L1_B_tumor", "Anti_PD_L1_C_tumor",
                                                                                "Anti_PD_L1_D_tumor", "Chemo_tumor", 
                                                                                "Untreated_A_tumor", "Untreated_B_tumor"), 
                                          to = c("COVID-19", "COVID-19", "COVID-19",
                                                 "Tumor", "Tumor", "Tumor", "Tumor", "Tumor", "Tumor", "Tumor"))

cd4_big_data[["percent.mt"]] <- PercentageFeatureSet(cd4_big_data, pattern = "^MT-")
cd4_big_data[['percent.ribo']] <- PercentageFeatureSet(cd4_big_data, pattern = "^RP[SL]")
VlnPlot(cd4_big_data, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(cd4_big_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cd4_big_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cd4_big_data <- subset(cd4_big_data, subset = nCount_RNA < 20000)

cd4_big_data <- NormalizeData(cd4_big_data, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(cd4_big_data, file = "cd4_big_data_Normalized.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_Normalized.rds")

#cd4_big_data <- subset(cd4_big_data, subset = CD3D > 0 | CD3G > 0 | CD3E > 0 | TRAC > 0)
#cd4_big_data <- subset(cd4_big_data, subset = CD8A == 0 | CD8B == 0)

cd4_big_data <- FindVariableFeatures(cd4_big_data, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 2))

top10 <- head(VariableFeatures(cd4_big_data), 10)
plot1 <- VariableFeaturePlot(cd4_big_data)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

saveRDS(cd4_big_data, file = "cd4_big_data_HVG.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_HVG.rds")

cd4_big_data <- ScaleData(cd4_big_data)

saveRDS(cd4_big_data, file = "cd4_big_data_ScaleData.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_ScaleData.rds")

cd4_big_data <- RunPCA(cd4_big_data, npcs = 20)

DimPlot(cd4_big_data, reduction = "pca")

ElbowPlot(cd4_big_data, ndims = 20)

saveRDS(cd4_big_data, file = "cd4_big_data_PCA.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_PCA.rds")

cd4_big_data <- RunTSNE(object = cd4_big_data, dims = 1:20, perplexity = 50, check_duplicates = FALSE)
TSNEPlot(object = cd4_big_data, label = TRUE, pt.size = 0.5)
FeaturePlot(cd4_big_data, features = "nCount_RNA")

cd4_big_data <- RunUMAP(object = cd4_big_data, dims = 1:20, check_duplicates = FALSE)
DimPlot(cd4_big_data, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(cd4_big_data, file = "cd4_big_data_TSNE.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_TSNE.rds")

cd4_big_data <- FindNeighbors(cd4_big_data, dims = 1:20)
res <- c(0, 0.4, 0.6, 0.8, 1.2, 1.4, 1.8, 2, 3)
cd4_big_data <- FindClusters(cd4_big_data, resolution = res, n.start = 100)

clustree(cd4_big_data, prefix = "RNA_snn_res.")#to choose the proper resolution I will use clustree and some of cell markers
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD14", node_colour_aggr = "median")#phages
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "HBB", node_colour_aggr = "median")#rbc
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "MS4A1", node_colour_aggr = "median")#b cells, cd20
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD19", node_colour_aggr = "median")#b cells
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells


clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "sc3_stability", node_colour_aggr = "median")
clustree_overlay(cd4_big_data, prefix = "RNA_snn_res.",red_dim = "tsne", x_value = "tsne1", y_value = "tsne2")
clustree_overlay(cd4_big_data, prefix = "RNA_snn_res.",red_dim = "umap", x_value = "umap1", y_value = "umap2")
clustree_overlay(cd4_big_data, prefix = "RNA_snn_res.",red_dim = "pca", x_value = "pca1", y_value = "pca2")

saveRDS(cd4_big_data, file = "cd4_big_data_Clusters.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_Clusters.rds")

cd4_big_data_heatmap <- clustify(cd4_big_data, cbmc_ref, cluster_col = "RNA_snn_res.1.2", obj_out = FALSE)
plot_cor_heatmap(cd4_big_data_heatmap)

cd4_big_data <- clustify(cd4_big_data, cbmc_ref, cluster_col = "RNA_snn_res.1.2", obj_out = TRUE, threshold = 0.15)
DimPlot(cd4_big_data, reduction = "tsne", group.by = "type")

cd4_big_data <- subset(cd4_big_data, subset = type == "CD4 T")

saveRDS(cd4_big_data, file = "cd4_big_data_res_1_2.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_res_1_2.rds")

cd4_big_data <- FindVariableFeatures(cd4_big_data, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 2))

top10 <- head(VariableFeatures(cd4_big_data), 10)
plot1 <- VariableFeaturePlot(cd4_big_data)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

saveRDS(cd4_big_data, file = "cd4_big_data_HVG_reclustered.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_HVG_reclustered.rds")

cd4_big_data <- ScaleData(cd4_big_data)

saveRDS(cd4_big_data, file = "cd4_big_data_ScaleData_reclustered.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_ScaleData_reclustered.rds")

cd4_big_data <- RunPCA(cd4_big_data, npcs = 20)

DimPlot(cd4_big_data, reduction = "pca")

ElbowPlot(cd4_big_data, ndims = 20)

saveRDS(cd4_big_data, file = "cd4_big_data_PCA_reclustered.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_PCA_reclustered.rds")

cd4_big_data <- FindNeighbors(cd4_big_data, dims = 1:20)
cd4_big_data <- FindClusters(cd4_big_data, resolution = res, n.start = 100)

clustree(cd4_big_data, prefix = "RNA_snn_res.")
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "sc3_stability", node_colour_aggr = "median")
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells

Idents(cd4_big_data) <- cd4_big_data$RNA_snn_res.0.8

cd4_big_data.markers <- list()

for (n in 1:13) {
  cd4_big_data.markers[[n]] <- FindMarkers(cd4_big_data, ident.1 = n, only.pos = TRUE)
}

cd4_big_data.cons.markers <- list()

for (n in 1:13) {
  cd4_big_data.cons.markers[[n]] <- FindConservedMarkers(cd4_big_data, ident.1 = n, grouping.var = "")
}

cd4_big_data_heatmap <- clustify(cd4_big_data, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(cd4_big_data_heatmap)

cd4_big_data <- clustify(cd4_big_data, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = TRUE)
DimPlot(cd4_big_data, reduction = "tsne", group.by = "type.clustify")

cd4_big_data <- subset(cd4_big_data, subset = type.clustify == "CD4 T" | type.clustify == "unassigned")

saveRDS(cd4_big_data, file = "cd4_big_data_recluster.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_recluster.rds")

cd4_big_data <- FindVariableFeatures(cd4_big_data, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 2))

top10 <- head(VariableFeatures(cd4_big_data), 10)
plot1 <- VariableFeaturePlot(cd4_big_data)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

saveRDS(cd4_big_data, file = "cd4_big_data_HVG_rereclustered.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_HVG_rereclustered.rds")

cd4_big_data <- ScaleData(cd4_big_data)

saveRDS(cd4_big_data, file = "cd4_big_data_ScaleData_rereclustered.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_ScaleData_rereclustered.rds")

cd4_big_data <- RunPCA(cd4_big_data, npcs = 20)

DimPlot(cd4_big_data, reduction = "pca")

ElbowPlot(cd4_big_data, ndims = 20)

saveRDS(cd4_big_data, file = "cd4_big_data_PCA_rereclustered.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_PCA_rereclustered.rds")

cd4_big_data <- FindNeighbors(cd4_big_data, dims = 1:20)
cd4_big_data <- FindClusters(cd4_big_data, resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 
                                                          1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9), n.start = 100)

saveRDS(cd4_big_data, file = "cd4_big_data_rerecluster.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_rerecluster.rds")

clustree(cd4_big_data, prefix = "RNA_snn_res.")
clustree(cd4_big_data, prefix = "RNA_snn_res.", node_colour = "sc3_stability", node_colour_aggr = "median")

cd4_big_data_heatmap <- clustify(cd4_big_data, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(cd4_big_data_heatmap)

cd4_big_data <- clustify(cd4_big_data, cbmc_ref, cluster_col = "RNA_snn_res.1.2", obj_out = TRUE, threshold = 0.15)
DimPlot(cd4_big_data, reduction = "tsne", group.by = "type")



cd4_big_data.markers <- FindAllMarkers(cd4_big_data)
EnhancedVolcano(cd4_big_data.markers, lab = rownames(cd4_big_data.markers), x = "", y = "")

Cd4_big_data <- clustify(cd4_big_data, ref_hema_microarray, cluster_col = "RNA_snn_res.0.8", obj_out = TRUE)
DimPlot(Cd4_big_data, reduction = "tsne", group.by = "type")

markers.ctl <- FindMarkers(Cd4_big_data, ident.1 = "CD8+ Effector Memory")
Idents(Cd4_big_data) <- Cd4_big_data$type

current.cluster.ids <- c("CD4+ Effector Memory", "CD8+ Effector Memory", "Monocyte", "unassined")
new.cluster.ids <- c("CD4+ Effector Memory", "CD4+ CTL", "Monocyte", "unassined")
Idents(object = Cd4_big_data) <- plyr::mapvalues(x = Idents(object = Cd4_big_data), from = current.cluster.ids, to = new.cluster.ids)

saveRDS(Cd4_big_data, file = "Cd4_big_data_recluster.rds")
Cd4_big_data <- readRDS(file = "Cd4_big_data_recluster.rds")

#cd4_big_data <- RenameIdents(cd4_big_data)