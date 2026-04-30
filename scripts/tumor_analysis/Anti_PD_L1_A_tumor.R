Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_seurat_object.rds")
Anti_PD_L1_A_tumor.metadata <- read.csv("GSE149652/Anti-PD-L1_A_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Anti_PD_L1_A_tumor.metadata) <- Anti_PD_L1_A_tumor.metadata[,1]
Anti_PD_L1_A_tumor.index <- rownames(Anti_PD_L1_A_tumor.metadata)
rownames(Anti_PD_L1_A_tumor.metadata) <- gsub("-", ".", Anti_PD_L1_A_tumor.index)

Anti_PD_L1_A_tumor <- AddMetaData(Anti_PD_L1_A_tumor, metadata = Anti_PD_L1_A_tumor.metadata)

Anti_PD_L1_A_tumor$nCount_RNA <- Anti_PD_L1_A_tumor$n_counts
Anti_PD_L1_A_tumor$nFeature_RNA <- Anti_PD_L1_A_tumor$n_genes

#Anti_PD_L1_A_tumor[["percent.mt"]] <- PercentageFeatureSet(Anti_PD_L1_A_tumor, pattern = "^MT-")
#Anti_PD_L1_A_tumor[['percent.ribo']] <- PercentageFeatureSet(Anti_PD_L1_A_tumor, pattern = "^RP[SL]")
VlnPlot(Anti_PD_L1_A_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"))#, "percent.mt", "percent.ribo"
FeatureScatter(Anti_PD_L1_A_tumor, feature1 = "nCount_RNA", feature2 = "percent_mito")
FeatureScatter(Anti_PD_L1_A_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Anti_PD_L1_A_tumor <- NormalizeData(Anti_PD_L1_A_tumor, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(Anti_PD_L1_A_tumor, file = "Anti_PD_L1_A_tumor_Normalized.rds")
Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_Normalized.rds")

Anti_PD_L1_A_tumor <- subset(Anti_PD_L1_A_tumor, subset = CD8A == 0)
Anti_PD_L1_A_tumor <- subset(Anti_PD_L1_A_tumor, subset = CD79A == 0)
Anti_PD_L1_A_tumor <- subset(Anti_PD_L1_A_tumor, subset = LYZ == 0)
Anti_PD_L1_A_tumor <- subset(Anti_PD_L1_A_tumor, subset = CST3 == 0)

Anti_PD_L1_A_tumor <- FindVariableFeatures(Anti_PD_L1_A_tumor, selection.method = "mean.var.plot", 
                                           dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Anti_PD_L1_A_tumor), 10)
plot1 <- VariableFeaturePlot(Anti_PD_L1_A_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

write.csv(VariableFeatures(Anti_PD_L1_A_tumor), "Anti_PD_L1_A_tumor_HVG.csv")

saveRDS(Anti_PD_L1_A_tumor, file = "Anti_PD_L1_A_tumor_HVG.rds")
Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_HVG.rds")

all.genes <- rownames(Anti_PD_L1_A_tumor)
Anti_PD_L1_A_tumor <- ScaleData(Anti_PD_L1_A_tumor, features = all.genes)
Anti_PD_L1_A_tumor <- RunPCA(Anti_PD_L1_A_tumor, npcs = 300)
DimPlot(Anti_PD_L1_A_tumor, reduction = "pca")
ElbowPlot(Anti_PD_L1_A_tumor, ndims = 300)

saveRDS(Anti_PD_L1_A_tumor, file = "Anti_PD_L1_A_tumor_PCA.rds")
Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_PCA.rds")

Anti_PD_L1_A_tumor <- RunTSNE(object = Anti_PD_L1_A_tumor, dims = 1:5)
TSNEPlot(object = Anti_PD_L1_A_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Anti_PD_L1_A_tumor, features = "nCount_RNA")

Anti_PD_L1_A_tumor <- RunUMAP(object = Anti_PD_L1_A_tumor, dims = 1:5)
DimPlot(Anti_PD_L1_A_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(Anti_PD_L1_A_tumor, file = "Anti_PD_L1_A_tumor_TSNE_UMAP.rds")
Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_TSNE_UMAP.rds")

Anti_PD_L1_A_tumor <- FindNeighbors(Anti_PD_L1_A_tumor, dims = 1:15)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
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

Idents(Anti_PD_L1_A_tumor) <- Anti_PD_L1_A_tumor$RNA_snn_res.0.7

Anti_PD_L1_A_tumor <- BuildClusterTree(object = Anti_PD_L1_A_tumor, dims = 1:15, reorder = TRUE, 
                                       reorder.numeric = TRUE)
PlotClusterTree(object = Anti_PD_L1_A_tumor)

Anti_PD_L1_A_tumor.markers <- FindAllMarkers(Anti_PD_L1_A_tumor)
m <- Anti_PD_L1_A_tumor.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Anti_PD_L1_A_tumor, features = m$gene)
write.csv(Anti_PD_L1_A_tumor.markers, "Anti_PD_L1_A_tumor_all_markers_before_names.csv")

Anti_PD_L1_A_tumor.table_idents <- as.data.frame(table(Idents(Anti_PD_L1_A_tumor), Anti_PD_L1_A_tumor$cell_types))
write.csv(Anti_PD_L1_A_tumor.table_idents, "Anti_PD_L1_A_tumor_table_idents_before_names.csv")

TSNEPlot(object = Anti_PD_L1_A_tumor, label = TRUE, pt.size = 0.5)

Anti_PD_L1_A_tumor <- RenameIdents(Anti_PD_L1_A_tumor, "1" = "Cytotoxic CD4", "2" = "Tregs", "4" = "Cytotoxic CD4",
                                   "3" = "Tcxcl13")

Anti_PD_L1_A_tumor <- BuildClusterTree(object = Anti_PD_L1_A_tumor, dims = 1:15, reorder = TRUE)
PlotClusterTree(object = Anti_PD_L1_A_tumor)

Anti_PD_L1_A_tumor$my_clusters <- Anti_PD_L1_A_tumor@active.ident

Anti_PD_L1_A_tumor.markers <- FindAllMarkers(Anti_PD_L1_A_tumor)
write.csv(Anti_PD_L1_A_tumor.markers, "Anti_PD_L1_A_tumor_all_markers_names.csv")

saveRDS(Anti_PD_L1_A_tumor, "Anti_PD_L1_A_tumor_clustered_names.rds")
Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_clustered_names.rds")

Anti_PD_L1_A_tumor_cytotoxic <- subset(Anti_PD_L1_A_tumor, subset = my_clusters == "Cytotoxic CD4")

saveRDS(Anti_PD_L1_A_tumor_cytotoxic, "Anti_PD_L1_A_tumor_cytotoxic.rds")
Anti_PD_L1_A_tumor_cytotoxic <- readRDS(file = "Anti_PD_L1_A_tumor_cytotoxic.rds")

Anti_PD_L1_A_tumor_heatmap <- clustify(Anti_PD_L1_A_tumor, ref_hema_microarray, cluster_col = "RNA_snn_res.0.7", 
                                       obj_out = FALSE)
plot_cor_heatmap(Anti_PD_L1_A_tumor_heatmap)
