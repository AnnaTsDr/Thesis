#Analysis of cd4t24 dataset, covid-19 paper data, pbmc
cd4t24 <- readRDS(file = "cd4t24_seurat_object.rds")
cd4t24.metadata <- read.delim("GSE152522/cd4t24_annotation.txt")

cd4t24[["percent.mt"]] <- PercentageFeatureSet(cd4t24, pattern = "^MT-")
cd4t24[['percent.ribo']] <- PercentageFeatureSet(cd4t24, pattern = "^RP[SL]")
VlnPlot(cd4t24, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))
FeatureScatter(cd4t24, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cd4t24, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cell_name <- colnames(cd4t24)
rownames(cd4t24.metadata) <- cell_name

cd4t24 <- AddMetaData(cd4t24, metadata = cd4t24.metadata)

#The data is after QC
cd4t24 <- NormalizeData(cd4t24, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(cd4t24, file = "cd4t24_Normalized.rds")
cd4t24 <- readRDS(file = "cd4t24_Normalized.rds")

cd4t24 <- subset(cd4t24, subset = CD8A == 0)
cd4t24 <- subset(cd4t24, subset = CD79A == 0)
cd4t24 <- subset(cd4t24, subset = LYZ == 0)
cd4t24 <- subset(cd4t24, subset = CST3 == 0)

cd4t24 <- FindVariableFeatures(cd4t24, selection.method = "mean.var.plot", dispersion.cutoff = c(0.5, Inf), 
                              mean.cutoff = c(0.0125, 5))
top10 <- head(VariableFeatures(cd4t24), 10)
plot1 <- VariableFeaturePlot(cd4t24)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

write.csv(VariableFeatures(cd4t24), "cd4t24_HVG.csv")

saveRDS(cd4t24, file = "cd4t24_HVG.rds")
cd4t24 <- readRDS(file = "cd4t24_HVG.rds")

all.genes <- rownames(cd4t24)
cd4t24 <- ScaleData(cd4t24, features = all.genes)
cd4t24 <- RunPCA(cd4t24, npcs = 20)
DimPlot(cd4t24, reduction = "pca")
ElbowPlot(cd4t24, ndims = 20)

saveRDS(cd4t24, file = "cd4t24_PCA.rds")
cd4t24 <- readRDS(file = "cd4t24_PCA.rds")

cd4t24 <- RunTSNE(object = cd4t24, dims = 1:20)
TSNEPlot(object = cd4t24, label = TRUE, pt.size = 0.5)
FeaturePlot(cd4t24, features = "nCount_RNA", cols = c("green", "black"))

cd4t24 <- RunUMAP(object = cd4t24, dims = 1:20)
DimPlot(cd4t24, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(cd4t24, file = "cd4t24_TSNE_UMAP.rds")
cd4t24 <- readRDS(file = "cd4t24_TSNE_UMAP.rds")

cd4t24 <- FindNeighbors(cd4t24, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
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

Idents(cd4t24) <- cd4t24$RNA_snn_res.0.3

cd4t24 <- BuildClusterTree(object = cd4t24, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = cd4t24)

cd4t24.markers <- FindAllMarkers(cd4t24)
m <- cd4t24.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(cd4t24, features = m$gene)
write.csv(cd4t24.markers, "cd4t24_all_markers_before_names.csv")

cd4t24.table_idents <- as.data.frame(table(Idents(cd4t24)))
write.csv(cd4t24.table_idents, "cd4t24_table_idents_before_names.csv")

cd4t24 <- RenameIdents(cd4t24, "3" = "Cytotoxic CD4", "2" = "Cytotoxic CD4", "1" = "Cytotoxic CD4", "5" = "Tfh",
                       "9" = "Tregs", "8" = "Tregs", "4" = "TCM", "6" = "Th17", "7" = "Th1")
cd4t24 <- BuildClusterTree(object = cd4t24, dims = 1:20, reorder = TRUE)
PlotClusterTree(object = cd4t24)

cd4t24$my_clusters <- cd4t24@active.ident

cd4t24.markers <- FindAllMarkers(cd4t24)
write.csv(cd4t24.markers, "cd4t24_all_markers_names.csv")

saveRDS(cd4t24, "cd4t24_clustered_names.rds")
cd4t24 <- readRDS(file = "cd4t24_clustered_names.rds")

cd4t24_cytotoxic <- subset(cd4t24, subset == my_clusters = "Cytotoxic CD4")

saveRDS(cd4t24_cytotoxic, "cd4t24_cytotoxic.rds")
cd4t24_cytotoxic <- readRDS(file = "cd4t24_cytotoxic.rds")

cd4t24_heatmap <- clustify(cd4t24, ref_hema_microarray, cluster_col = "my_clusters", obj_out = FALSE)
plot_cor_heatmap(cd4t24_heatmap)
