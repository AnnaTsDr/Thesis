library(MazamaCoreUtils)
library(pbapply)
library(prob)
source(file = "AssessNodesAdaptedfromSeurat.R")

#Analysis of cd4t0n6 dataset, covid-19 paper data, cd4 t cells only
cd4t0n6 <- readRDS(file = "cd4t0n6_seurat_object.rds")

cd4t0n6$species <- "Human"
cd4t0n6$condition <- "Bladder_Cancer"

cd4t0n6[["percent.mt"]] <- PercentageFeatureSet(cd4t0n6, pattern = "^MT-")
cd4t0n6[['percent.ribo']] <- PercentageFeatureSet(cd4t0n6, pattern = "^RP[SL]")
VlnPlot(cd4t0n6, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))
FeatureScatter(cd4t0n6, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cd4t0n6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cd4t0n6 <- NormalizeData(cd4t0n6, normalization.method = "LogNormalize", scale.factor = 10000)
cd4t0n6 <- FindVariableFeatures(cd4t0n6, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                          mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(cd4t0n6), 10)
plot1 <- VariableFeaturePlot(cd4t0n6)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

cd4t0n6 <- ScaleData(cd4t0n6)
cd4t0n6 <- RunPCA(cd4t0n6, npcs = 20)
DimPlot(cd4t0n6, reduction = "pca")
ElbowPlot(cd4t0n6, ndims = 20)

cd4t0n6 <- RunTSNE(object = cd4t0n6, dims = 1:20)
TSNEPlot(object = cd4t0n6, label = TRUE, pt.size = 0.5)
FeaturePlot(cd4t0n6, features = "nCount_RNA")

cd4t0n6 <- RunUMAP(object = cd4t0n6, dims = 1:20)
DimPlot(cd4t0n6, label = TRUE, pt.size = 0.5, reduction = "umap")

cd4t0n6 <- FindNeighbors(cd4t0n6, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
cd4t0n6 <- FindClusters(cd4t0n6, resolution = res)

clustree(cd4t0n6, prefix = "RNA_snn_res.")
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(cd4t0n6, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(cd4t0n6, "cd4t0n6_clustered.rds")
cd4t0n6 <- readRDS(file = "cd4t0n6_clustered.rds")

Idents(cd4t0n6) <- cd4t0n6$RNA_snn_res.0.8

cd4t0n6 <- BuildClusterTree(object = cd4t0n6, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = cd4t0n6)

cd4t0n6.markers <- list()

for (n in 1:17) {
  cd4t0n6.markers[[n]] <- FindMarkers(cd4t0n6, ident.1 = n, only.pos = TRUE)
}

table(Idents(cd4t0n6))
TSNEPlot(object = cd4t0n6, label = TRUE, pt.size = 0.5)

cd4t0n6 <- RenameIdents(cd4t0n6, "21" = "11", "20" = "11", "19" = "11", "18" = "11", "17" = "11", "16" = "11",
                        "15" = "11", "14" = "11", "13" = "11", "12" = "11", "10" = "9", "7" = "4", "6" = "4",
                        "5" = "4")
cd4t0n6 <- BuildClusterTree(object = cd4t0n6, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = cd4t0n6)

cd4t0n6$my_clusters <- cd4t0n6@active.ident

cd4t0n6.markers <- FindAllMarkers(cd4t0n6)

cd4t0n6_heatmap <- clustify(cd4t0n6, cbmc_ref, cluster_col = "my_clusters", obj_out = FALSE)
plot_cor_heatmap(cd4t0n6_heatmap)

cd4t6 <- readRDS(file = "cd4t6_seurat_object.rds")

cd4t6$species <- "Human"
cd4t6$condition <- "Bladder_Cancer"

cd4t6[["percent.mt"]] <- PercentageFeatureSet(cd4t6, pattern = "^MT-")
cd4t6[['percent.ribo']] <- PercentageFeatureSet(cd4t6, pattern = "^RP[SL]")
VlnPlot(cd4t6, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(cd4t6, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cd4t6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cd4t6 <- NormalizeData(cd4t6, normalization.method = "LogNormalize", scale.factor = 10000)
cd4t6 <- FindVariableFeatures(cd4t6, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(cd4t6), 10)
plot1 <- VariableFeaturePlot(cd4t6)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

cd4t6 <- ScaleData(cd4t6)
cd4t6 <- RunPCA(cd4t6, npcs = 20)
DimPlot(cd4t6, reduction = "pca")
ElbowPlot(cd4t6, ndims = 20)

cd4t6 <- RunTSNE(object = cd4t6, dims = 1:20)
TSNEPlot(object = cd4t6, label = TRUE, pt.size = 0.5)
FeaturePlot(cd4t6, features = "nCount_RNA")

cd4t6 <- RunUMAP(object = cd4t6, dims = 1:20)
DimPlot(cd4t6, label = TRUE, pt.size = 0.5, reduction = "umap")

cd4t6 <- FindNeighbors(cd4t6, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
cd4t6 <- FindClusters(cd4t6, resolution = res)

clustree(cd4t6, prefix = "RNA_snn_res.")
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(cd4t6, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(cd4t6, "cd4t6_clustered.rds")
cd4t6 <- readRDS(file = "cd4t6_clustered.rds")

Idents(cd4t6) <- cd4t6$RNA_snn_res.0.8

cd4t6 <- BuildClusterTree(object = cd4t6, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = cd4t6)

cd4t6.markers <- list()

for (n in 1:13) {
  cd4t6.markers[[n]] <- FindMarkers(cd4t6, ident.1 = n, only.pos = TRUE)
}

table(Idents(cd4t6))
TSNEPlot(object = cd4t6, label = TRUE, pt.size = 0.5)

cd4t6.markers <- FindAllMarkers(cd4t6)

cd4t6 <- RenameIdents(cd4t6, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

cd4t6_heatmap <- clustify(cd4t6, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(cd4t6_heatmap)

cd4t24 <- readRDS(file = "cd4t24_seurat_object.rds")

cd4t24$species <- "Human"
cd4t24$condition <- "Bladder_Cancer"

cd4t24[["percent.mt"]] <- PercentageFeatureSet(cd4t24, pattern = "^MT-")
cd4t24[['percent.ribo']] <- PercentageFeatureSet(cd4t24, pattern = "^RP[SL]")
VlnPlot(cd4t24, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(cd4t24, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cd4t24, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cd4t24 <- NormalizeData(cd4t24, normalization.method = "LogNormalize", scale.factor = 10000)
cd4t24 <- FindVariableFeatures(cd4t24, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(cd4t24), 10)
plot1 <- VariableFeaturePlot(cd4t24)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

cd4t24 <- ScaleData(cd4t24)
cd4t24 <- RunPCA(cd4t24, npcs = 20)
DimPlot(cd4t24, reduction = "pca")
ElbowPlot(cd4t24, ndims = 20)

cd4t24 <- RunTSNE(object = cd4t24, dims = 1:20)
TSNEPlot(object = cd4t24, label = TRUE, pt.size = 0.5)
FeaturePlot(cd4t24, features = "nCount_RNA")

cd4t24 <- RunUMAP(object = cd4t24, dims = 1:20)
DimPlot(cd4t24, label = TRUE, pt.size = 0.5, reduction = "umap")

cd4t24 <- FindNeighbors(cd4t24, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
cd4t24 <- FindClusters(cd4t24, resolution = res)

clustree(cd4t24, prefix = "RNA_snn_res.")
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(cd4t24, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(cd4t24, "cd4t24_clustered.rds")
cd4t24 <- readRDS(file = "cd4t24_clustered.rds")

Idents(cd4t24) <- cd4t24$RNA_snn_res.0.8

cd4t24 <- BuildClusterTree(object = cd4t24, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = cd4t24)

cd4t24.markers <- list()

for (n in 1:13) {
  cd4t24.markers[[n]] <- FindMarkers(cd4t24, ident.1 = n, only.pos = TRUE)
}

table(Idents(cd4t24))
TSNEPlot(object = cd4t24, label = TRUE, pt.size = 0.5)

cd4t24.markers <- FindAllMarkers(cd4t24)

cd4t24 <- RenameIdents(cd4t24, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

cd4t24_heatmap <- clustify(cd4t24, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(cd4t24_heatmap)

Supercentenerians <- readRDS(file = "Supercentenerians_seurat_object.rds")

Supercentenerians$species <- "Human"

Supercentenerians <- AddMetaData(Supercentenerians, metadata = Metadata.Supercentenerians)

Supercentenerians$condition <- Supercentenerians$V3

Supercentenerians$condition <- plyr::mapvalues(x = Supercentenerians$condition, from = c("SC", "CT"), 
                                          to = c("Supercentenarians", "Control"))
Supercentenerians$V3 <- NULL
Supercentenerians$V2 <- NULL

Supercentenerians[["percent.mt"]] <- PercentageFeatureSet(Supercentenerians, pattern = "^MT-")
Supercentenerians[['percent.ribo']] <- PercentageFeatureSet(Supercentenerians, pattern = "^RP[SL]")
VlnPlot(Supercentenerians, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Supercentenerians, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Supercentenerians, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Supercentenerians <- NormalizeData(Supercentenerians, normalization.method = "LogNormalize", scale.factor = 10000)
Supercentenerians <- FindVariableFeatures(Supercentenerians, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Supercentenerians), 10)
plot1 <- VariableFeaturePlot(Supercentenerians)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

Supercentenerians <- ScaleData(Supercentenerians)
Supercentenerians <- RunPCA(Supercentenerians, npcs = 20)
DimPlot(Supercentenerians, reduction = "pca")
ElbowPlot(Supercentenerians, ndims = 20)

Supercentenerians <- RunTSNE(object = Supercentenerians, dims = 1:20)
TSNEPlot(object = Supercentenerians, label = TRUE, pt.size = 0.5)
FeaturePlot(Supercentenerians, features = "nCount_RNA")

Supercentenerians <- RunUMAP(object = Supercentenerians, dims = 1:20)
DimPlot(Supercentenerians, label = TRUE, pt.size = 0.5, reduction = "umap")

Supercentenerians <- FindNeighbors(Supercentenerians, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
Supercentenerians <- FindClusters(Supercentenerians, resolution = res)

clustree(Supercentenerians, prefix = "RNA_snn_res.")
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "CD14", node_colour_aggr = "median")#phages, cd14
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "CD19", node_colour_aggr = "median")#b cells
clustree(Supercentenerians, prefix = "RNA_snn_res.", node_colour = "MS4A1", node_colour_aggr = "median")#b cells, cd20

saveRDS(Supercentenerians, "Supercentenerians_clustered.rds")
Supercentenerians <- readRDS(file = "Supercentenerians_clustered.rds")

Idents(Supercentenerians) <- Supercentenerians$RNA_snn_res.0.8

Supercentenerians <- BuildClusterTree(object = Supercentenerians, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = Supercentenerians)

Supercentenerians.markers <- list()

for (n in 1:21) {
  Supercentenerians.markers[[n]] <- FindMarkers(Supercentenerians, ident.1 = n, only.pos = TRUE)
}

table(Idents(Supercentenerians))
TSNEPlot(object = Supercentenerians, label = TRUE, pt.size = 0.5)

Supercentenerians.markers <- FindAllMarkers(Supercentenerians)

Supercentenerians <- RenameIdents(Supercentenerians, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

Supercentenerians_heatmap <- clustify(Supercentenerians, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(Supercentenerians_heatmap)

Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_seurat_object.rds")

Anti_PD_L1_A_tumor$species <- "Human"
Anti_PD_L1_A_tumor$condition <- "Covid-19"

Anti_PD_L1_A_tumor[["percent.mt"]] <- PercentageFeatureSet(Anti_PD_L1_A_tumor, pattern = "^MT-")
Anti_PD_L1_A_tumor[['percent.ribo']] <- PercentageFeatureSet(Anti_PD_L1_A_tumor, pattern = "^RP[SL]")
VlnPlot(Anti_PD_L1_A_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Anti_PD_L1_A_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Anti_PD_L1_A_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Anti_PD_L1_A_tumor <- NormalizeData(Anti_PD_L1_A_tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Anti_PD_L1_A_tumor <- FindVariableFeatures(Anti_PD_L1_A_tumor, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Anti_PD_L1_A_tumor), 10)
plot1 <- VariableFeaturePlot(Anti_PD_L1_A_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

Anti_PD_L1_A_tumor <- ScaleData(Anti_PD_L1_A_tumor)
Anti_PD_L1_A_tumor <- RunPCA(Anti_PD_L1_A_tumor, npcs = 20)
DimPlot(Anti_PD_L1_A_tumor, reduction = "pca")
ElbowPlot(Anti_PD_L1_A_tumor, ndims = 20)

Anti_PD_L1_A_tumor <- RunTSNE(object = Anti_PD_L1_A_tumor, dims = 1:20)
TSNEPlot(object = Anti_PD_L1_A_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Anti_PD_L1_A_tumor, features = "nCount_RNA")

Anti_PD_L1_A_tumor <- RunUMAP(object = Anti_PD_L1_A_tumor, dims = 1:20)
DimPlot(Anti_PD_L1_A_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

Anti_PD_L1_A_tumor <- FindNeighbors(Anti_PD_L1_A_tumor, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
Anti_PD_L1_A_tumor <- FindClusters(Anti_PD_L1_A_tumor, resolution = res)

clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.")
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(Anti_PD_L1_A_tumor, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(Anti_PD_L1_A_tumor, "Anti_PD_L1_A_tumor_clustered.rds")
Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_clustered.rds")

Idents(Anti_PD_L1_A_tumor) <- Anti_PD_L1_A_tumor$RNA_snn_res.0.8

Anti_PD_L1_A_tumor <- BuildClusterTree(object = Anti_PD_L1_A_tumor, dims = 1:20, reorder = TRUE, 
                                       reorder.numeric = TRUE)
PlotClusterTree(object = Anti_PD_L1_A_tumor)

Anti_PD_L1_A_tumor.markers <- list()

for (n in 1:4) {
  Anti_PD_L1_A_tumor.markers[[n]] <- FindMarkers(Anti_PD_L1_A_tumor, ident.1 = n, only.pos = TRUE)
}

table(Idents(Anti_PD_L1_A_tumor))
TSNEPlot(object = Anti_PD_L1_A_tumor, label = TRUE, pt.size = 0.5)

Anti_PD_L1_A_tumor.markers <- FindAllMarkers(Anti_PD_L1_A_tumor)

Anti_PD_L1_A_tumor <- RenameIdents(Anti_PD_L1_A_tumor, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

Anti_PD_L1_A_tumor_heatmap <- clustify(Anti_PD_L1_A_tumor, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(Anti_PD_L1_A_tumor_heatmap)

Anti_PD_L1_B_tumor <- readRDS(file = "Anti_PD_L1_B_tumor_seurat_object.rds")

Anti_PD_L1_B_tumor$species <- "Human"
Anti_PD_L1_B_tumor$condition <- "Covid-19"

Anti_PD_L1_B_tumor[["percent.mt"]] <- PercentageFeatureSet(Anti_PD_L1_B_tumor, pattern = "^MT-")
Anti_PD_L1_B_tumor[['percent.ribo']] <- PercentageFeatureSet(Anti_PD_L1_B_tumor, pattern = "^RP[SL]")
VlnPlot(Anti_PD_L1_B_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Anti_PD_L1_B_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Anti_PD_L1_B_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Anti_PD_L1_B_tumor <- NormalizeData(Anti_PD_L1_B_tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Anti_PD_L1_B_tumor <- FindVariableFeatures(Anti_PD_L1_B_tumor, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Anti_PD_L1_B_tumor), 10)
plot1 <- VariableFeaturePlot(Anti_PD_L1_B_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

Anti_PD_L1_B_tumor <- ScaleData(Anti_PD_L1_B_tumor)
Anti_PD_L1_B_tumor <- RunPCA(Anti_PD_L1_B_tumor, npcs = 20)
DimPlot(Anti_PD_L1_B_tumor, reduction = "pca")
ElbowPlot(Anti_PD_L1_B_tumor, ndims = 20)

Anti_PD_L1_B_tumor <- RunTSNE(object = Anti_PD_L1_B_tumor, dims = 1:20)
TSNEPlot(object = Anti_PD_L1_B_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Anti_PD_L1_B_tumor, features = "nCount_RNA")

Anti_PD_L1_B_tumor <- RunUMAP(object = Anti_PD_L1_B_tumor, dims = 1:20)
DimPlot(Anti_PD_L1_B_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

Anti_PD_L1_B_tumor <- FindNeighbors(Anti_PD_L1_B_tumor, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
Anti_PD_L1_B_tumor <- FindClusters(Anti_PD_L1_B_tumor, resolution = res)

clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.")
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(Anti_PD_L1_B_tumor, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(Anti_PD_L1_B_tumor, "Anti_PD_L1_B_tumor_clustered.rds")
Anti_PD_L1_B_tumor <- readRDS(file = "Anti_PD_L1_B_tumor_clustered.rds")

Idents(Anti_PD_L1_B_tumor) <- Anti_PD_L1_B_tumor$RNA_snn_res.0.8

Anti_PD_L1_B_tumor <- BuildClusterTree(object = Anti_PD_L1_B_tumor, dims = 1:20, reorder = TRUE, 
                                       reorder.numeric = TRUE)
PlotClusterTree(object = Anti_PD_L1_B_tumor)

Anti_PD_L1_B_tumor.markers <- list()

for (n in 1:4) {
  Anti_PD_L1_B_tumor.markers[[n]] <- FindMarkers(Anti_PD_L1_B_tumor, ident.1 = n, only.pos = TRUE)
}

table(Idents(Anti_PD_L1_B_tumor))
TSNEPlot(object = Anti_PD_L1_B_tumor, label = TRUE, pt.size = 0.5)

Anti_PD_L1_B_tumor.markers <- FindAllMarkers(Anti_PD_L1_B_tumor)

Anti_PD_L1_B_tumor <- RenameIdents(Anti_PD_L1_B_tumor, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

Anti_PD_L1_B_tumor_heatmap <- clustify(Anti_PD_L1_B_tumor, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(Anti_PD_L1_B_tumor_heatmap)

Anti_PD_L1_C_tumor <- readRDS(file = "Anti_PD_L1_C_tumor_seurat_object.rds")

Anti_PD_L1_C_tumor$species <- "Human"
Anti_PD_L1_C_tumor$condition <- "Covid-19"

Anti_PD_L1_C_tumor[["percent.mt"]] <- PercentageFeatureSet(Anti_PD_L1_C_tumor, pattern = "^MT-")
Anti_PD_L1_C_tumor[['percent.ribo']] <- PercentageFeatureSet(Anti_PD_L1_C_tumor, pattern = "^RP[SL]")
VlnPlot(Anti_PD_L1_C_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Anti_PD_L1_C_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Anti_PD_L1_C_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Anti_PD_L1_C_tumor <- NormalizeData(Anti_PD_L1_C_tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Anti_PD_L1_C_tumor <- FindVariableFeatures(Anti_PD_L1_C_tumor, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Anti_PD_L1_C_tumor), 10)
plot1 <- VariableFeaturePlot(Anti_PD_L1_C_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

Anti_PD_L1_C_tumor <- ScaleData(Anti_PD_L1_C_tumor)
Anti_PD_L1_C_tumor <- RunPCA(Anti_PD_L1_C_tumor, npcs = 20)
DimPlot(Anti_PD_L1_C_tumor, reduction = "pca")
ElbowPlot(Anti_PD_L1_C_tumor, ndims = 20)

Anti_PD_L1_C_tumor <- RunTSNE(object = Anti_PD_L1_C_tumor, dims = 1:20)
TSNEPlot(object = Anti_PD_L1_C_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Anti_PD_L1_C_tumor, features = "nCount_RNA")

Anti_PD_L1_C_tumor <- RunUMAP(object = Anti_PD_L1_C_tumor, dims = 1:20)
DimPlot(Anti_PD_L1_C_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

Anti_PD_L1_C_tumor <- FindNeighbors(Anti_PD_L1_C_tumor, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
Anti_PD_L1_C_tumor <- FindClusters(Anti_PD_L1_C_tumor, resolution = res)

clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.")
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(Anti_PD_L1_C_tumor, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(Anti_PD_L1_C_tumor, "Anti_PD_L1_C_tumor_clustered.rds")
Anti_PD_L1_C_tumor <- readRDS(file = "Anti_PD_L1_C_tumor_clustered.rds")

Idents(Anti_PD_L1_C_tumor) <- Anti_PD_L1_C_tumor$RNA_snn_res.0.8

Anti_PD_L1_C_tumor <- BuildClusterTree(object = Anti_PD_L1_C_tumor, dims = 1:20, reorder = TRUE, 
                                       reorder.numeric = TRUE)
PlotClusterTree(object = Anti_PD_L1_C_tumor)

Anti_PD_L1_C_tumor.markers <- list()

for (n in 1:4) {
  Anti_PD_L1_C_tumor.markers[[n]] <- FindMarkers(Anti_PD_L1_C_tumor, ident.1 = n, only.pos = TRUE)
}

table(Idents(Anti_PD_L1_C_tumor))
TSNEPlot(object = Anti_PD_L1_C_tumor, label = TRUE, pt.size = 0.5)

Anti_PD_L1_C_tumor.markers <- FindAllMarkers(Anti_PD_L1_C_tumor)

Anti_PD_L1_C_tumor <- RenameIdents(Anti_PD_L1_C_tumor, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

Anti_PD_L1_C_tumor_heatmap <- clustify(Anti_PD_L1_C_tumor, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(Anti_PD_L1_C_tumor_heatmap)

Anti_PD_L1_D_tumor <- readRDS(file = "Anti_PD_L1_D_tumor_seurat_object.rds")

Anti_PD_L1_D_tumor$species <- "Human"
Anti_PD_L1_D_tumor$condition <- "Covid-19"

Anti_PD_L1_D_tumor[["percent.mt"]] <- PercentageFeatureSet(Anti_PD_L1_D_tumor, pattern = "^MT-")
Anti_PD_L1_D_tumor[['percent.ribo']] <- PercentageFeatureSet(Anti_PD_L1_D_tumor, pattern = "^RP[SL]")
VlnPlot(Anti_PD_L1_D_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Anti_PD_L1_D_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Anti_PD_L1_D_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Anti_PD_L1_D_tumor <- NormalizeData(Anti_PD_L1_D_tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Anti_PD_L1_D_tumor <- FindVariableFeatures(Anti_PD_L1_D_tumor, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Anti_PD_L1_D_tumor), 10)
plot1 <- VariableFeaturePlot(Anti_PD_L1_D_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

Anti_PD_L1_D_tumor <- ScaleData(Anti_PD_L1_D_tumor)
Anti_PD_L1_D_tumor <- RunPCA(Anti_PD_L1_D_tumor, npcs = 20)
DimPlot(Anti_PD_L1_D_tumor, reduction = "pca")
ElbowPlot(Anti_PD_L1_D_tumor, ndims = 20)

Anti_PD_L1_D_tumor <- RunTSNE(object = Anti_PD_L1_D_tumor, dims = 1:20)
TSNEPlot(object = Anti_PD_L1_D_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Anti_PD_L1_D_tumor, features = "nCount_RNA")

Anti_PD_L1_D_tumor <- RunUMAP(object = Anti_PD_L1_D_tumor, dims = 1:20)
DimPlot(Anti_PD_L1_D_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

Anti_PD_L1_D_tumor <- FindNeighbors(Anti_PD_L1_D_tumor, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
Anti_PD_L1_D_tumor <- FindClusters(Anti_PD_L1_D_tumor, resolution = res)

clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.")
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(Anti_PD_L1_D_tumor, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(Anti_PD_L1_D_tumor, "Anti_PD_L1_D_tumor_clustered.rds")
Anti_PD_L1_D_tumor <- readRDS(file = "Anti_PD_L1_D_tumor_clustered.rds")

Idents(Anti_PD_L1_D_tumor) <- Anti_PD_L1_D_tumor$RNA_snn_res.0.8

Anti_PD_L1_D_tumor <- BuildClusterTree(object = Anti_PD_L1_D_tumor, dims = 1:20, reorder = TRUE, 
                                       reorder.numeric = TRUE)
PlotClusterTree(object = Anti_PD_L1_D_tumor)

Anti_PD_L1_D_tumor.markers <- list()

for (n in 1:6) {
  Anti_PD_L1_D_tumor.markers[[n]] <- FindMarkers(Anti_PD_L1_D_tumor, ident.1 = n, only.pos = TRUE)
}

table(Idents(Anti_PD_L1_D_tumor))
TSNEPlot(object = Anti_PD_L1_D_tumor, label = TRUE, pt.size = 0.5)

Anti_PD_L1_D_tumor.markers <- FindAllMarkers(Anti_PD_L1_D_tumor)

Anti_PD_L1_D_tumor <- RenameIdents(Anti_PD_L1_D_tumor, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

Anti_PD_L1_D_tumor_heatmap <- clustify(Anti_PD_L1_D_tumor, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(Anti_PD_L1_D_tumor_heatmap)

Chemo_tumor <- readRDS(file = "Chemo_tumor_seurat_object.rds")

Chemo_tumor$species <- "Human"
Chemo_tumor$condition <- "Covid-19"

Chemo_tumor[["percent.mt"]] <- PercentageFeatureSet(Chemo_tumor, pattern = "^MT-")
Chemo_tumor[['percent.ribo']] <- PercentageFeatureSet(Chemo_tumor, pattern = "^RP[SL]")
VlnPlot(Chemo_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Chemo_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Chemo_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Chemo_tumor <- NormalizeData(Chemo_tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Chemo_tumor <- FindVariableFeatures(Chemo_tumor, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Chemo_tumor), 10)
plot1 <- VariableFeaturePlot(Chemo_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

Chemo_tumor <- ScaleData(Chemo_tumor)
Chemo_tumor <- RunPCA(Chemo_tumor, npcs = 20)
DimPlot(Chemo_tumor, reduction = "pca")
ElbowPlot(Chemo_tumor, ndims = 20)

Chemo_tumor <- RunTSNE(object = Chemo_tumor, dims = 1:20)
TSNEPlot(object = Chemo_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Chemo_tumor, features = "nCount_RNA")

Chemo_tumor <- RunUMAP(object = Chemo_tumor, dims = 1:20)
DimPlot(Chemo_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

Chemo_tumor <- FindNeighbors(Chemo_tumor, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
Chemo_tumor <- FindClusters(Chemo_tumor, resolution = res)

clustree(Chemo_tumor, prefix = "RNA_snn_res.")
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(Chemo_tumor, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(Chemo_tumor, "Chemo_tumor_clustered.rds")
Chemo_tumor <- readRDS(file = "Chemo_tumor_clustered.rds")

Idents(Chemo_tumor) <- Chemo_tumor$RNA_snn_res.0.8

Chemo_tumor <- BuildClusterTree(object = Chemo_tumor, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = Chemo_tumor)

Chemo_tumor.markers <- list()

for (n in 1:4) {
  Chemo_tumor.markers[[n]] <- FindMarkers(Chemo_tumor, ident.1 = n, only.pos = TRUE)
}

table(Idents(Chemo_tumor))
TSNEPlot(object = Chemo_tumor, label = TRUE, pt.size = 0.5)

Chemo_tumor.markers <- FindAllMarkers(Chemo_tumor)

Chemo_tumor <- RenameIdents(Chemo_tumor, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

Chemo_tumor_heatmap <- clustify(Chemo_tumor, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(Chemo_tumor_heatmap)

Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_seurat_object.rds")

Untreated_A_tumor$species <- "Human"
Untreated_A_tumor$condition <- "Covid-19"

Untreated_A_tumor[["percent.mt"]] <- PercentageFeatureSet(Untreated_A_tumor, pattern = "^MT-")
Untreated_A_tumor[['percent.ribo']] <- PercentageFeatureSet(Untreated_A_tumor, pattern = "^RP[SL]")
VlnPlot(Untreated_A_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Untreated_A_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Untreated_A_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Untreated_A_tumor <- NormalizeData(Untreated_A_tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Untreated_A_tumor <- FindVariableFeatures(Untreated_A_tumor, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Untreated_A_tumor), 10)
plot1 <- VariableFeaturePlot(Untreated_A_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

Untreated_A_tumor <- ScaleData(Untreated_A_tumor)
Untreated_A_tumor <- RunPCA(Untreated_A_tumor, npcs = 20)
DimPlot(Untreated_A_tumor, reduction = "pca")
ElbowPlot(Untreated_A_tumor, ndims = 20)

Untreated_A_tumor <- RunTSNE(object = Untreated_A_tumor, dims = 1:20)
TSNEPlot(object = Untreated_A_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Untreated_A_tumor, features = "nCount_RNA")

Untreated_A_tumor <- RunUMAP(object = Untreated_A_tumor, dims = 1:20)
DimPlot(Untreated_A_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

Untreated_A_tumor <- FindNeighbors(Untreated_A_tumor, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
Untreated_A_tumor <- FindClusters(Untreated_A_tumor, resolution = res)

clustree(Untreated_A_tumor, prefix = "RNA_snn_res.")
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(Untreated_A_tumor, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(Untreated_A_tumor, "Untreated_A_tumor.rds")
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor.rds")

Idents(Untreated_A_tumor) <- Untreated_A_tumor$RNA_snn_res.0.8

Untreated_A_tumor <- BuildClusterTree(object = Untreated_A_tumor, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = Untreated_A_tumor)

Untreated_A_tumor.markers <- list()

for (n in 1:6) {
  Untreated_A_tumor.markers[[n]] <- FindMarkers(Untreated_A_tumor, ident.1 = n, only.pos = TRUE)
}

table(Idents(Untreated_A_tumor))
TSNEPlot(object = Untreated_A_tumor, label = TRUE, pt.size = 0.5)

Untreated_A_tumor.markers <- FindAllMarkers(Untreated_A_tumor)

Untreated_A_tumor <- RenameIdents(Untreated_A_tumor, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

Untreated_A_tumor_heatmap <- clustify(Untreated_A_tumor, cbmc_ref, cluster_col = "RNA_snn_res.0.8", obj_out = FALSE)
plot_cor_heatmap(Untreated_A_tumor_heatmap)

Untreated_B_tumor <- readRDS(file = "Untreated_B_tumor_seurat_object.rds")

Untreated_B_tumor$species <- "Human"
Untreated_B_tumor$condition <- "Covid-19"

Untreated_B_tumor[["percent.mt"]] <- PercentageFeatureSet(Untreated_B_tumor, pattern = "^MT-")
Untreated_B_tumor[['percent.ribo']] <- PercentageFeatureSet(Untreated_B_tumor, pattern = "^RP[SL]")
VlnPlot(Untreated_B_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Untreated_B_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Untreated_B_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Untreated_B_tumor <- NormalizeData(Untreated_B_tumor, normalization.method = "LogNormalize", scale.factor = 10000)
Untreated_B_tumor <- FindVariableFeatures(Untreated_B_tumor, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Untreated_B_tumor), 10)
plot1 <- VariableFeaturePlot(Untreated_B_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

Untreated_B_tumor <- ScaleData(Untreated_B_tumor)
Untreated_B_tumor <- RunPCA(Untreated_B_tumor, npcs = 20)
DimPlot(Untreated_B_tumor, reduction = "pca")
ElbowPlot(Untreated_B_tumor, ndims = 20)

Untreated_B_tumor <- RunTSNE(object = Untreated_B_tumor, dims = 1:20)
TSNEPlot(object = Untreated_B_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Untreated_B_tumor, features = "nCount_RNA")

Untreated_B_tumor <- RunUMAP(object = Untreated_B_tumor, dims = 1:20)
DimPlot(Untreated_B_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

Untreated_B_tumor <- FindNeighbors(Untreated_B_tumor, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
Untreated_B_tumor <- FindClusters(Untreated_B_tumor, resolution = res)

clustree(Untreated_B_tumor, prefix = "RNA_snn_res.")
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "PTPRC", node_colour_aggr = "median")#cd45 all wbc, all cells in the dataset are wbc
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD3D", node_colour_aggr = "median")# t cells, all cels in this dataset are t cells, t cells refers to cd4, cd8 and nk
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD3G", node_colour_aggr = "median")# t cells
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD3E", node_colour_aggr = "median")#t cells
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "KLRF1", node_colour_aggr = "median")#nk
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD8A", node_colour_aggr = "median")#cd8 t cells, no cd8 t cells in this dataset
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD8B", node_colour_aggr = "median")#cd8 t cells
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "CD4", node_colour_aggr = "median")#cd4 t cells
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "FCGR3A", node_colour_aggr = "median")#phages, cd16
clustree(Untreated_B_tumor, prefix = "RNA_snn_res.", node_colour = "HBA1", node_colour_aggr = "median")#rbc

saveRDS(Untreated_B_tumor, "Untreated_B_tumor.rds")
Untreated_B_tumor <- readRDS(file = "Untreated_B_tumor.rds")

Idents(Untreated_B_tumor) <- Untreated_B_tumor$RNA_snn_res.0.8

Untreated_B_tumor <- BuildClusterTree(object = Untreated_B_tumor, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = Untreated_B_tumor)

Untreated_B_tumor.markers <- list()

for (n in 1:2) {
  Untreated_B_tumor.markers[[n]] <- FindMarkers(Untreated_B_tumor, ident.1 = n, only.pos = TRUE)
}

table(Idents(Untreated_B_tumor))
TSNEPlot(object = Untreated_B_tumor, label = TRUE, pt.size = 0.5)

Untreated_B_tumor.markers <- FindAllMarkers(Untreated_B_tumor)

Untreated_B_tumor <- RenameIdents(Untreated_B_tumor, "12" = "4", "18" = "4", "7" = "17", "1" = "0", "10" = "2", "13" = "3", "14" = "5",
                        "9" = "5")

Untreated_B_tumor_heatmap <- clustify(Untreated_B_tumor, ref_hema_microarray, cluster_col = "seurat_clusters", obj_out = FALSE)
plot_cor_heatmap(Untreated_B_tumor_heatmap)

#Hezis data
raw_counts <- read.csv(file = "RawData1.csv", sep = ",", row.names = 1, header = TRUE)
object.size(raw_counts)
dim(raw_counts)

agingdata <- CreateSeuratObject(counts = raw_counts, project = "Aging_mouse")
object.size(agingdata)

saveRDS(agingdata, file = "AgingDataSeuratObject.rds")
agingdata <- readRDS("AgingDataSeuratObject.rds")

genes <- rownames(agingdata)
write.csv(genes, "aging_genes.csv")
genes <- read.csv("aging_genes.csv")
row.names(agingdata@assays$RNA) <- genes[,3]

agingdata$species <- "Mouse"

agingdata[["percent.mt"]] <- PercentageFeatureSet(agingdata, pattern = "^MT-")
agingdata[['percent.ribo']] <- PercentageFeatureSet(agingdata, pattern = "^RP[SL]")
p1 <- VlnPlot(agingdata, features = "nFeature_RNA", pt.size = 0) + labs(title = "nGenes", tag = "A") + NoLegend()
p2 <- VlnPlot(agingdata, features = "nCount_RNA", pt.size = 0) + labs(title = "nUMI", tag = "B") + NoLegend()
p3 <- VlnPlot(agingdata, features = "percent.mt", pt.size = 0) + labs(title = "Mito percent", tag = "C") + NoLegend()
p4 <- VlnPlot(agingdata, features = "percent.ribo", pt.size = 0) + labs(title = "Ribo percent", tag = "C") + NoLegend()
p1 <- AugmentPlot(plot = p1) + labs(y = "nGenes")
p2 <- AugmentPlot(plot = p2) + labs(y = "nUMI")
p3 <- AugmentPlot(plot = p3) + labs(y = "Mito percrnt")
p4 <- AugmentPlot(plot = p4) + labs(y = "Ribo percrnt")
p1 + p2 + p3 + p4
VlnPlot(agingdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Metadata.aging <- read.csv(file = "Metadata.csv")

Mouse <- as.array(Metadata.aging$Mouse)
Batch <- as.array(Metadata.aging$Batch)
Age_group <- as.array(Metadata.aging$Age_group)
Subset <- as.array(Metadata.aging$Subset)

agingdata <- AddMetaData(agingdata, metadata = Mouse, col.name = "Mouse")
agingdata <- AddMetaData(agingdata, metadata = Batch, col.name = "Batch")
agingdata <- AddMetaData(agingdata, metadata = Age_group, col.name = "Age_group")
agingdata <- AddMetaData(agingdata, metadata = Subset, col.name = "Subset_Idan")

agingdata$condition <- agingdata$Age_group

saveRDS(agingdata, file = "AgingDataMetadataAdded.rds")
agingdata <- readRDS("AgingDataMetadataAdded.rds")

plot1 <- FeatureScatter(agingdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(agingdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

agingdata <- NormalizeData(agingdata, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(agingdata, file = "AgingDataNormalized.rds")
agingdata <- readRDS(file = "AgingDataNormalized.rds")

agingdata <- FindVariableFeatures(agingdata, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                  mean.cutoff = c(0.0125, 2))

top10 <- head(VariableFeatures(agingdata), 10)
plot1 <- VariableFeaturePlot(agingdata)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

saveRDS(agingdata, file = "AgingDataHVG.rds")
write.csv(VariableFeatures(agingdata), file = "AgingDataHVG.csv")
agingdata <- readRDS(file = "AgingDataHVG.rds")

all.genes <- rownames(agingdata)
agingdata <- ScaleData(agingdata, features = all.genes)

saveRDS(agingdata, file = "AgingDataScaled.rds")
agingdata <- readRDS(file = "AgingDataScaled.rds")

agingdata <- RunPCA(agingdata, npcs = 20)
DimPlot(agingdata, reduction = "pca")
ElbowPlot(agingdata, ndims = 20)

agingdata <- RunTSNE(object = agingdata, dims = 1:20)
TSNEPlot(object = agingdata, label = TRUE, pt.size = 0.5)
FeaturePlot(agingdata, features = "nCount_RNA")

agingdata <- RunUMAP(agingdata, dims = 1:20)
DimPlot(agingdata, reduction = "umap")

saveRDS(agingdata, file = "AgingData&tSNE&UMAP.rds")
agingdata <- readRDS("AgingData&tSNE&UMAP.rds")

agingdata <- FindNeighbors(agingdata, dims = 1:20)
agingdata <- FindClusters(agingdata, resolution = 0.8)

saveRDS(agingdata, file = "AgingDataClustered.rds")

agingdata <- BuildClusterTree(agingdata, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = agingdata)
