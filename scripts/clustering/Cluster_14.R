Cluster_14 <- readRDS(file = "TMS_merge_T_cells_clusters_14.rds")

#Highly variable expressed genes
Cluster_14 <- FindVariableFeatures(Cluster_14, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_14),15)
plot1 <- VariableFeaturePlot(Cluster_14)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_14), file = "TMS_merge_T_cells_Cluster_14_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_14, file = "TMS_merge_T_cells_Cluster_14_VariableGenes.rds")
Cluster_14 <- readRDS(file = "TMS_merge_T_cells_Cluster_14_VariableGenes.rds")

all.genes <- rownames(Cluster_14)
Cluster_14 <- ScaleData(Cluster_14, features = all.genes)

saveRDS(Cluster_14, file = "TMS_merge_T_cells_Cluster_14_ScaledData.rds")
Cluster_14 <- readRDS(file = "TMS_merge_T_cells_Cluster_14_ScaledData.rds")

Cluster_14 <- RunPCA(Cluster_14, npcs = 7, ndims.print = 1:7, nfeatures.print = 5)
DimPlot(Cluster_14, reduction = "pca")

saveRDS(Cluster_14, file = "TMS_merge_T_cells_Cluster_14_PCA.rds")
Cluster_14 <- readRDS(file = "TMS_merge_T_cells_Cluster_14_PCA.rds")

Cluster_14 <- JackStraw(Cluster_14, num.replicate = 100, dims = 7)
Cluster_14 <- ScoreJackStraw(Cluster_14, dims = 1:7)
JackStrawPlot(Cluster_14, dims = 1:7)

ElbowPlot(Cluster_14, ndims = 7)

Cluster_14 <- RunTSNE(Cluster_14, dims = 1:7, perplexity = 10)
DimPlot(Cluster_14, reduction = "tsne")

Cluster_14 <- RunUMAP(Cluster_14, dims = 1:7)
DimPlot(Cluster_14, reduction = "umap")

saveRDS(Cluster_14, file = "TMS_merge_T_cells_Cluster_14_TSNE_UMAP.rds")
Cluster_14 <- readRDS(file = "TMS_merge_T_cells_Cluster_14_TSNE_UMAP.rds")

Cluster_14 <- FindNeighbors(Cluster_14, dims = 1:7)

#louvian
Cluster_14 <- FindClusters(Cluster_14, resolution = seq(0,3, by = 0.1))

clustree(Cluster_14, prefix = "RNA_snn_res.")
clustree(Cluster_14, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_14) <- Cluster_14$RNA_snn_res.3

Cluster_14 <- BuildClusterTree(Cluster_14, reorder.numeric = TRUE, reorder = TRUE, dims = 1:7)
PlotClusterTree(object = Cluster_14)

Cluster_14.markers <- FindAllMarkers(Cluster_14)
m <- Cluster_14.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_14, features = m$gene)
write.csv(Cluster_14.markers, "TMS_T_cells_Cluster_14_merkers.csv")

Cluster_14 <- RenameIdents(Cluster_14, "1" = "CD8 T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD8 T",
                           "6" = "CD8 T", "7" = "CD8 T", "8" = "DP T", "9" = "DP T", "10" = "CD8 T", "11" = "CD8 T",
                           "12" = "CD8 T", "13" = "CD8 T")

saveRDS(Cluster_14, file = "TMS_merge_T_cells_Cluster_14_clustered.rds")
Cluster_14 <- readRDS(file = "TMS_merge_T_cells_Cluster_14_clustered.rds")

