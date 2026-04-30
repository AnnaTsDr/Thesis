#Analysis of cd4t0n6 dataset, covid-19 paper data, pbmc
cd4t0n6 <- readRDS(file = "cd4t0n6_seurat_object.rds")
cd4t0n6.metadata <- read.delim("GSE152522/cd4t0n6_annotation.csv")

cd4t0n6[["percent.mt"]] <- PercentageFeatureSet(cd4t0n6, pattern = "^MT-")
cd4t0n6[['percent.ribo']] <- PercentageFeatureSet(cd4t0n6, pattern = "^RP[SL]")
VlnPlot(cd4t0n6, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))
FeatureScatter(cd4t0n6, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cd4t0n6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cell_name <- colnames(cd4t0n6)
rownames(cd4t0n6.metadata) <- cell_name

cd4t0n6 <- AddMetaData(cd4t0n6, metadata = cd4t0n6.metadata)

#The data is after QC
cd4t0n6 <- NormalizeData(cd4t0n6, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(cd4t0n6, file = "cd4t0n6_Normalized.rds")
cd4t0n6 <- readRDS(file = "cd4t0n6_Normalized.rds")

cd4t0n6 <- subset(cd4t0n6, subset = CD8A == 0)
cd4t0n6 <- subset(cd4t0n6, subset = CD79A == 0)
cd4t0n6 <- subset(cd4t0n6, subset = LYZ == 0)
cd4t0n6 <- subset(cd4t0n6, subset = CST3 == 0)

cd4t0n6 <- FindVariableFeatures(cd4t0n6, selection.method = "mean.var.plot", dispersion.cutoff = c(0.7, Inf), 
                                mean.cutoff = c(0.0125, 5))
top10 <- head(VariableFeatures(cd4t0n6), 10)
plot1 <- VariableFeaturePlot(cd4t0n6)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

write.csv(VariableFeatures(cd4t0n6), "cd4t0n6_HVG.csv")

saveRDS(cd4t0n6, file = "cd4t0n6_HVG.rds")
cd4t0n6 <- readRDS(file = "cd4t0n6_HVG.rds")

all.genes <- rownames(cd4t0n6)
cd4t0n6 <- ScaleData(cd4t0n6, features = all.genes)
cd4t0n6 <- RunPCA(cd4t0n6, npcs = 20)
DimPlot(cd4t0n6, reduction = "pca")
ElbowPlot(cd4t0n6, ndims = 20)

saveRDS(cd4t0n6, file = "cd4t0n6_PCA.rds")
cd4t0n6 <- readRDS(file = "cd4t0n6_PCA.rds")

cd4t0n6 <- RunTSNE(object = cd4t0n6, dims = 1:20)
TSNEPlot(object = cd4t0n6, label = TRUE, pt.size = 0.5)
FeaturePlot(cd4t0n6, features = "nCount_RNA", cols = c("green", "black"))

cd4t0n6 <- RunUMAP(object = cd4t0n6, dims = 1:20)
DimPlot(cd4t0n6, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(cd4t0n6, file = "cd4t0n6_TSNE_UMAP.rds")
cd4t0n6 <- readRDS(file = "cd4t0n6_TSNE_UMAP.rds")

cd4t0n6 <- FindNeighbors(cd4t0n6, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
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

Idents(cd4t0n6) <- cd4t0n6$RNA_snn_res.0.5

cd4t0n6 <- BuildClusterTree(object = cd4t0n6, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = cd4t0n6)

cd4t0n6.markers <- FindAllMarkers(cd4t0n6)
m <- cd4t0n6.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(cd4t0n6, features = m$gene)
write.csv(cd4t0n6.markers, "cd4t0n6_all_markers_before_names.csv")

cd4t0n6.table_idents <- as.data.frame(table(Idents(cd4t0n6)))
write.csv(cd4t0n6.table_idents, "cd4t0n6_table_idents_before_names.csv")

cd4t0n6 <- RenameIdents(cd4t0n6, "10" = "7", "9" = "7", "8" = "5", "2" = "1", "11" = "6", "15" = "12", "14" = "12", 
                        "13" = "12")
cd4t0n6 <- BuildClusterTree(object = cd4t0n6, dims = 1:20, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = cd4t0n6)

cd4t0n6 <- RenameIdents(cd4t0n6, "1" = "MKI67", "2" = "Cytotoxic CD4", "3" = "IFN responce", "4" = "Th1", "5" = "Tfh" ,
                        "6" = "TCM", "7" = "Th17")

cd4t0n6 <- BuildClusterTree(object = cd4t0n6, dims = 1:20, reorder = TRUE)
PlotClusterTree(object = cd4t0n6)

cd4t0n6$my_clusters <- cd4t0n6@active.ident

cd4t0n6.markers <- FindAllMarkers(cd4t0n6)
write.csv(cd4t0n6.markers, "cd4t0n6_all_markers_names.csv")

saveRDS(cd4t0n6, "cd4t0n6_clustered_names.rds")
cd4t0n6 <- readRDS(file = "cd4t0n6_clustered_names.rds")

cd4t0n6_cytotoxic <- subset(cd4t0n6, subset = my_clusters == "Cytotoxic CD4")

saveRDS(cd4t0n6_cytotoxic, "cd4t0n6_cytotoxic.rds")
cd4t0n6_cytotoxic <- readRDS(file = "cd4t0n6_cytotoxic.rds")

cd4t0n6_heatmap <- clustify(cd4t0n6, cbmc_ref, cluster_col = "my_clusters", obj_out = FALSE)
plot_cor_heatmap(cd4t0n6_heatmap)
