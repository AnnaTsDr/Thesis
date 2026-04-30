#Analysis of cd4t6 dataset, covid-19 paper data, pbmc
cd4t6 <- readRDS(file = "cd4t6_seurat_object.rds")
cd4t6.metadata <- read.delim("GSE152522/cd4t6_annotation.txt")

cd4t6[["percent.mt"]] <- PercentageFeatureSet(cd4t6, pattern = "^MT-")
cd4t6[['percent.ribo']] <- PercentageFeatureSet(cd4t6, pattern = "^RP[SL]")
VlnPlot(cd4t6, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))
FeatureScatter(cd4t6, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cd4t6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cell_name <- colnames(cd4t6)
rownames(cd4t6.metadata) <- cell_name

cd4t6 <- AddMetaData(cd4t6, metadata = cd4t6.metadata)

#The data is after QC
cd4t6 <- NormalizeData(cd4t6, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(cd4t6, file = "cd4t6_Normalized.rds")
cd4t6 <- readRDS(file = "cd4t6_Normalized.rds")

cd4t6 <- subset(cd4t6, subset = CD8A == 0)
cd4t6 <- subset(cd4t6, subset = CD79A == 0)
cd4t6 <- subset(cd4t6, subset = LYZ == 0)
cd4t6 <- subset(cd4t6, subset = CST3 == 0)

cd4t6 <- FindVariableFeatures(cd4t6, selection.method = "mean.var.plot", dispersion.cutoff = c(0.7, Inf), 
                                mean.cutoff = c(0.0125, 5))
top10 <- head(VariableFeatures(cd4t6), 10)
plot1 <- VariableFeaturePlot(cd4t6)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

write.csv(VariableFeatures(cd4t6), "cd4t6_HVG.csv")

saveRDS(cd4t6, file = "cd4t6_HVG.rds")
cd4t6 <- readRDS(file = "cd4t6_HVG.rds")

all.genes <- rownames(cd4t6)
cd4t6 <- ScaleData(cd4t6, features = all.genes)
cd4t6 <- RunPCA(cd4t6, npcs = 20)
DimPlot(cd4t6, reduction = "pca")
ElbowPlot(cd4t6, ndims = 20)

saveRDS(cd4t6, file = "cd4t6_PCA.rds")
cd4t6 <- readRDS(file = "cd4t6_PCA.rds")

cd4t6 <- RunTSNE(object = cd4t6, dims = 1:20)
TSNEPlot(object = cd4t6, label = TRUE, pt.size = 0.5)
FeaturePlot(cd4t6, features = "nCount_RNA", cols = c("green", "black"))

cd4t6 <- RunUMAP(object = cd4t6, dims = 1:20)
DimPlot(cd4t6, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(cd4t6, file = "cd4t6_TSNE_UMAP.rds")
cd4t6 <- readRDS(file = "cd4t6_TSNE_UMAP.rds")

cd4t6 <- FindNeighbors(cd4t6, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
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

Idents(cd4t6) <- cd4t6$RNA_snn_res.0.3

cd4t6 <- BuildClusterTree(object = cd4t6, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = cd4t6)

cd4t6.markers <- FindAllMarkers(cd4t6)
m <- cd4t6.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(cd4t6, features = m$gene)
write.csv(cd4t6.markers, "cd4t6_all_markers_before_names.csv")

cd4t6.table_idents <- as.data.frame(table(Idents(cd4t6)))
write.csv(cd4t6.table_idents, "cd4t6_table_idents_before_names.csv")

cd4t6 <- RenameIdents(cd4t6, "3" = "Cytotoxic CD4", "2" = "Th1", "5" = "Th1", "4" = "Tfh", "1" = "Tregs", "6" = "TCM",
                      "7" = "Tfh")
cd4t6 <- BuildClusterTree(object = cd4t6, dims = 1:20, reorder = TRUE)
PlotClusterTree(object = cd4t6)

cd4t6$my_clusters <- cd4t6@active.ident

cd4t6.markers <- FindAllMarkers(cd4t6)
write.csv(cd4t6.markers, "cd4t6_all_markers_names.csv")

saveRDS(cd4t6, "cd4t6_clustered_names.rds")
cd4t6 <- readRDS(file = "cd4t6_clustered_names.rds")

cd4t6_cytotoxic <- subset(cd4t6, subset = my_clusters == "Cytotoxic CD4")

saveRDS(cd4t6_cytotoxic, "cd4t6_cytotoxic.rds")
cd4t6_cytotoxic <- readRDS(file = "cd4t6_cytotoxic.rds")

cd4t6_heatmap <- clustify(cd4t6, cbmc_ref, cluster_col = "my_clusters", obj_out = FALSE)
plot_cor_heatmap(cd4t6_heatmap)
