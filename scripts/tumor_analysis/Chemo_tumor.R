Chemo_tumor <- readRDS(file = "Chemo_tumor_seurat_object.rds")
Chemo_tumor.metadata <- read.csv("GSE149652/Anti-PD-L1_A_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Chemo_tumor.metadata) <- Chemo_tumor.metadata[,1]
Chemo_tumor.index <- rownames(Chemo_tumor.metadata)
rownames(Chemo_tumor.metadata) <- gsub("-", ".", Chemo_tumor.index)

Chemo_tumor <- AddMetaData(Chemo_tumor, metadata = Chemo_tumor.metadata)

Chemo_tumor$nCount_RNA <- Chemo_tumor$n_counts
Chemo_tumor$nFeature_RNA <- Chemo_tumor$n_genes

Chemo_tumor[["percent.mt"]] <- PercentageFeatureSet(Chemo_tumor, pattern = "^MT-")
Chemo_tumor[['percent.ribo']] <- PercentageFeatureSet(Chemo_tumor, pattern = "^RP[SL]")
VlnPlot(Chemo_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Chemo_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Chemo_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Chemo_tumor <- NormalizeData(Chemo_tumor, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(Chemo_tumor, file = "Chemo_tumor_Normalized.rds")
Chemo_tumor <- readRDS(file = "Chemo_tumor_Normalized.rds")

Chemo_tumor <- subset(Chemo_tumor, subset = CD8A == 0)
Chemo_tumor <- subset(Chemo_tumor, subset = CD79A == 0)
Chemo_tumor <- subset(Chemo_tumor, subset = LYZ == 0)
Chemo_tumor <- subset(Chemo_tumor, subset = CST3 == 0)

Chemo_tumor <- FindVariableFeatures(Chemo_tumor, selection.method = "mean.var.plot", 
                                          dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Chemo_tumor), 10)
plot1 <- VariableFeaturePlot(Chemo_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

write.csv(VariableFeatures(Chemo_tumor), "Chemo_tumor_HVG.csv")

saveRDS(Chemo_tumor, file = "Chemo_tumor_HVG.rds")
Chemo_tumor <- readRDS(file = "Chemo_tumor_HVG.rds")

all.genes <- rownames(Chemo_tumor)
Chemo_tumor <- ScaleData(Chemo_tumor, features = all.genes)
Chemo_tumor <- RunPCA(Chemo_tumor, npcs = 15)
DimPlot(Chemo_tumor, reduction = "pca")
ElbowPlot(Chemo_tumor, ndims = 15)

saveRDS(Chemo_tumor, file = "Chemo_tumor_PCA.rds")
Chemo_tumor <- readRDS(file = "Chemo_tumor_PCA.rds")

Chemo_tumor <- RunTSNE(object = Chemo_tumor, dims = 1:15)
TSNEPlot(object = Chemo_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Chemo_tumor, features = "nCount_RNA")

Chemo_tumor <- RunUMAP(object = Chemo_tumor, dims = 1:15)
DimPlot(Chemo_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(Chemo_tumor, file = "Chemo_tumor_TSNE_UMAP.rds")
Chemo_tumor <- readRDS(file = "Chemo_tumor_TSNE_UMAP.rds")

Chemo_tumor <- FindNeighbors(Chemo_tumor, dims = 1:15)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
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

Chemo_tumor <- BuildClusterTree(object = Chemo_tumor, dims = 1:15, reorder = TRUE, 
                                      reorder.numeric = TRUE)
PlotClusterTree(object = Chemo_tumor)

Chemo_tumor.markers <- FindAllMarkers(Chemo_tumor)
m <- Chemo_tumor.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Chemo_tumor, features = m$gene)
write.csv(Chemo_tumor.markers, "Chemo_tumor_all_markers_before_names.csv")

Chemo_tumor.table_idents <- as.data.frame(table(Idents(Chemo_tumor)))
write.csv(Chemo_tumor.table_idents, "Chemo_tumor_table_idents_before_names.csv")

TSNEPlot(object = Chemo_tumor, label = TRUE, pt.size = 0.5)

Chemo_tumor <- RenameIdents(Chemo_tumor, "1" = "Cytotoxic CD4", "2" = "Tregs")

Chemo_tumor <- BuildClusterTree(object = Chemo_tumor, dims = 1:15, reorder = TRUE)
PlotClusterTree(object = Chemo_tumor)

Chemo_tumor$my_clusters <- Chemo_tumor@active.ident

Chemo_tumor.markers <- FindAllMarkers(Chemo_tumor)
write.csv(Chemo_tumor.markers, "Chemo_tumor_all_markers_names.csv")

saveRDS(Chemo_tumor, "Chemo_tumor_clustered_names.rds")
Chemo_tumor <- readRDS(file = "Chemo_tumor_clustered_names.rds")

Chemo_tumor_cytotoxic <- subset(Chemo_tumor, subset = my_clusters == "Cytotoxic CD4")

saveRDS(Chemo_tumor_cytotoxic, "Chemo_tumor_cytotoxic.rds")
Chemo_tumor_cytotoxic <- readRDS(file = "Chemo_tumor_cytotoxic.rds")


Chemo_tumor_heatmap <- clustify(Chemo_tumor, cbmc_ref, cluster_col = "RNA_snn_res.0.7", obj_out = FALSE)
plot_cor_heatmap(Chemo_tumor_heatmap)
