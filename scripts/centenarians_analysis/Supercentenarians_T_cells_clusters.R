
Cluster_1 <- readRDS(file = "Supercentenarians_T_cells_clusters_1.rds")
Cluster_2 <- readRDS(file = "Supercentenarians_T_cells_clusters_2.rds")
Cluster_3 <- readRDS(file = "Supercentenarians_T_cells_clusters_3.rds")
Cluster_4 <- readRDS(file = "Supercentenarians_T_cells_clusters_4.rds")
Cluster_5 <- readRDS(file = "Supercentenarians_T_cells_clusters_5.rds")
Cluster_6 <- readRDS(file = "Supercentenarians_T_cells_clusters_6.rds")
Cluster_7 <- readRDS(file = "Supercentenarians_T_cells_clusters_7.rds")
Cluster_8 <- readRDS(file = "Supercentenarians_T_cells_clusters_8.rds")
Cluster_9 <- readRDS(file = "Supercentenarians_T_cells_clusters_9.rds")
Cluster_10 <- readRDS(file = "Supercentenarians_T_cells_clusters_10.rds")
Cluster_11 <- readRDS(file = "Supercentenarians_T_cells_clusters_11.rds")
Cluster_12 <- readRDS(file = "Supercentenarians_T_cells_clusters_12.rds")
Cluster_13 <- readRDS(file = "Supercentenarians_T_cells_clusters_13.rds")
Cluster_14 <- readRDS(file = "Supercentenarians_T_cells_clusters_14.rds")

#Highly variable expressed genes
Cluster_1 <- FindVariableFeatures(Cluster_1, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.6, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_1),15)
plot1 <- VariableFeaturePlot(Cluster_1)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_1), file = "Supercentenarians_T_Cluster_1_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_1, file = "Supercentenarians_T_Cluster_1_VariableGenes.rds")
Cluster_1 <- readRDS(file = "Supercentenarians_T_Cluster_1_VariableGenes.rds")

all.genes <- rownames(Cluster_1)
Cluster_1 <- ScaleData(Cluster_1, features = all.genes)

saveRDS(Cluster_1, file = "Supercentenarians_T_Cluster_1_ScaledData.rds")
Cluster_1 <- readRDS(file = "Supercentenarians_T_Cluster_1_ScaledData.rds")

Cluster_1 <- RunPCA(Cluster_1, npcs = 20, ndims.print = 1:20, nfeatures.print = 5)
DimPlot(Cluster_1, reduction = "pca")

saveRDS(Cluster_1, file = "Supercentenarians_T_Cluster_1_PCA.rds")
Cluster_1 <- readRDS(file = "Supercentenarians_T_Cluster_1_PCA.rds")

Cluster_1 <- JackStraw(Cluster_1, num.replicate = 100, dims = 20)
Cluster_1 <- ScoreJackStraw(Cluster_1, dims = 1:20)
JackStrawPlot(Cluster_1, dims = 1:20)

ElbowPlot(Cluster_1, ndims = 20)

Cluster_1 <- RunTSNE(Cluster_1, dims = 1:20, perplexity = 10)
DimPlot(Cluster_1, reduction = "tsne")

Cluster_1 <- RunUMAP(Cluster_1, dims = 1:20)
DimPlot(Cluster_1, reduction = "umap")

saveRDS(Cluster_1, file = "Supercentenarians_T_Cluster_1_TSNE_UMAP.rds")
Cluster_1 <- readRDS(file = "Supercentenarians_T_Cluster_1_TSNE_UMAP.rds")

Cluster_1 <- FindNeighbors(Cluster_1, dims = 1:20)

#louvian
Cluster_1 <- FindClusters(Cluster_1, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_1, prefix = "RNA_snn_res.")
clustree(Cluster_1, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_1) <- Cluster_1$RNA_snn_res.0.7

Cluster_1 <- BuildClusterTree(Cluster_1, reorder.numeric = TRUE, reorder = TRUE, dims = 1:3)
PlotClusterTree(object = Cluster_1)

Cluster_1.markers <- FindAllMarkers(Cluster_1)
m <- Cluster_1.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_1, features = m$gene)
write.csv(Cluster_1.markers, "Supercentenarians_T_Cluster_1_merkers.csv")

Cluster_1 <- RenameIdents(Cluster_1, "0" = "DN T")

saveRDS(Cluster_1, file = "Supercentenarians_T_Cluster_1_clustered.rds")
Cluster_1 <- readRDS(file = "Supercentenarians_T_Cluster_1_clustered.rds")