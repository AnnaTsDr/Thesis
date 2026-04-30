Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_seurat_object.rds")
Untreated_A_tumor.metadata <- read.csv("GSE149652/Anti-PD-L1_A_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Untreated_A_tumor.metadata) <- Untreated_A_tumor.metadata[,1]
Untreated_A_tumor.index <- rownames(Untreated_A_tumor.metadata)
rownames(Untreated_A_tumor.metadata) <- gsub("-", ".", Untreated_A_tumor.index)

Untreated_A_tumor <- AddMetaData(Untreated_A_tumor, metadata = Untreated_A_tumor.metadata)

Untreated_A_tumor$nCount_RNA <- Untreated_A_tumor$n_counts
Untreated_A_tumor$nFeature_RNA <- Untreated_A_tumor$n_genes

Untreated_A_tumor[["percent.mt"]] <- PercentageFeatureSet(Untreated_A_tumor, pattern = "^MT-")
Untreated_A_tumor[['percent.ribo']] <- PercentageFeatureSet(Untreated_A_tumor, pattern = "^RP[SL]")
VlnPlot(Untreated_A_tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Untreated_A_tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Untreated_A_tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Untreated_A_tumor <- NormalizeData(Untreated_A_tumor, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(Untreated_A_tumor, file = "Untreated_A_tumor_Normalized.rds")
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_Normalized.rds")

Untreated_A_tumor <- subset(Untreated_A_tumor, subset = CD8A == 0)
Untreated_A_tumor <- subset(Untreated_A_tumor, subset = CD79A == 0)
Untreated_A_tumor <- subset(Untreated_A_tumor, subset = LYZ == 0)
Untreated_A_tumor <- subset(Untreated_A_tumor, subset = CST3 == 0)

Untreated_A_tumor <- FindVariableFeatures(Untreated_A_tumor, selection.method = "mean.var.plot", 
                                          dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Untreated_A_tumor), 10)
plot1 <- VariableFeaturePlot(Untreated_A_tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

write.csv(VariableFeatures(Untreated_A_tumor), "Untreated_A_tumor_HVG.csv")

saveRDS(Untreated_A_tumor, file = "Untreated_A_tumor_HVG.rds")
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_HVG.rds")

all.genes <- rownames(Untreated_A_tumor)
Untreated_A_tumor <- ScaleData(Untreated_A_tumor, features = all.genes)
Untreated_A_tumor <- RunPCA(Untreated_A_tumor, npcs = 15)
DimPlot(Untreated_A_tumor, reduction = "pca")
ElbowPlot(Untreated_A_tumor, ndims = 15)

saveRDS(Untreated_A_tumor, file = "Untreated_A_tumor_PCA.rds")
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_PCA.rds")

Untreated_A_tumor <- RunTSNE(object = Untreated_A_tumor, dims = 1:15)
TSNEPlot(object = Untreated_A_tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(Untreated_A_tumor, features = "nCount_RNA")

Untreated_A_tumor <- RunUMAP(object = Untreated_A_tumor, dims = 1:15)
DimPlot(Untreated_A_tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(Untreated_A_tumor, file = "Untreated_A_tumor_TSNE_UMAP.rds")
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_TSNE_UMAP.rds")

Untreated_A_tumor <- FindNeighbors(Untreated_A_tumor, dims = 1:15)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
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

saveRDS(Untreated_A_tumor, "Untreated_A_tumor_clustered.rds")
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_clustered.rds")

Idents(Untreated_A_tumor) <- Untreated_A_tumor$RNA_snn_res.0.6

Untreated_A_tumor <- BuildClusterTree(object = Untreated_A_tumor, dims = 1:15, reorder = TRUE, 
                                      reorder.numeric = TRUE)
PlotClusterTree(object = Untreated_A_tumor)

Untreated_A_tumor.markers <- FindAllMarkers(Untreated_A_tumor)
m <- Untreated_A_tumor.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Untreated_A_tumor, features = m$gene)
write.csv(Untreated_A_tumor.markers, "Untreated_A_tumor_all_markers_before_names.csv")

Untreated_A_tumor.table_idents <- as.data.frame(table(Idents(Untreated_A_tumor)))
write.csv(Untreated_A_tumor.table_idents, "Untreated_A_tumor_table_idents_before_names.csv")

TSNEPlot(object = Untreated_A_tumor, label = TRUE, pt.size = 0.5)

Untreated_A_tumor <- RenameIdents(Untreated_A_tumor, "1" = "Cytotoxic CD4", "3" = "Tregs")

Untreated_A_tumor <- BuildClusterTree(object = Untreated_A_tumor, dims = 1:15, reorder = TRUE)
PlotClusterTree(object = Untreated_A_tumor)

Untreated_A_tumor$my_clusters <- Untreated_A_tumor@active.ident

Untreated_A_tumor.markers <- FindAllMarkers(Untreated_A_tumor)
write.csv(Untreated_A_tumor.markers, "Untreated_A_tumor_all_markers_names.csv")

Untreated_A_tumor.table_idents <- as.data.frame(table(Idents(Untreated_A_tumor)))
write.csv(Untreated_A_tumor.table_idents, "Untreated_A_tumor_table_idents_before_names.csv")

saveRDS(Untreated_A_tumor, "Untreated_A_tumor_clustered_names.rds")
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_clustered_names.rds")

Untreated_A_tumor_cytotoxic <- subset(Untreated_A_tumor, subset = my_clusters == "Cytotoxic CD4")

saveRDS(Untreated_A_tumor_cytotoxic, "Untreated_A_tumor_cytotoxic.rds")
Untreated_A_tumor_cytotoxic <- readRDS(file = "Untreated_A_tumor_cytotoxic.rds")


Untreated_A_tumor_heatmap <- clustify(Untreated_A_tumor, cbmc_ref, cluster_col = "RNA_snn_res.0.7", obj_out = FALSE)
plot_cor_heatmap(Untreated_A_tumor_heatmap)
