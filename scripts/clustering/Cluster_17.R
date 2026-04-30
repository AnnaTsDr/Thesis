Cluster_17 <- readRDS(file = "TMS_merge_T_cells_clusters_17.rds")

#Highly variable expressed genes
Cluster_17 <- FindVariableFeatures(Cluster_17, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_17),15)
plot1 <- VariableFeaturePlot(Cluster_17)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_17), file = "TMS_merge_T_cells_Cluster_17_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_17, file = "TMS_merge_T_cells_Cluster_17_VariableGenes.rds")
Cluster_17 <- readRDS(file = "TMS_merge_T_cells_Cluster_17_VariableGenes.rds")

all.genes <- rownames(Cluster_17)
Cluster_17 <- ScaleData(Cluster_17, features = all.genes)

saveRDS(Cluster_17, file = "TMS_merge_T_cells_Cluster_17_ScaledData.rds")
Cluster_17 <- readRDS(file = "TMS_merge_T_cells_Cluster_17_ScaledData.rds")

Cluster_17 <- RunPCA(Cluster_17, npcs = 2, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_17, reduction = "pca")

saveRDS(Cluster_17, file = "TMS_merge_T_cells_Cluster_17_PCA.rds")
Cluster_17 <- readRDS(file = "TMS_merge_T_cells_Cluster_17_PCA.rds")

Cluster_17 <- JackStraw(Cluster_17, num.replicate = 100, dims = 2)
Cluster_17 <- ScoreJackStraw(Cluster_17, dims = 1:2)
JackStrawPlot(Cluster_17, dims = 1:2)

ElbowPlot(Cluster_17, ndims = 2)

Cluster_17 <- RunTSNE(Cluster_17, dims = 1:2, perplexity = 10)
DimPlot(Cluster_17, reduction = "tsne")

Cluster_17 <- RunUMAP(Cluster_17, dims = 1:2)
DimPlot(Cluster_17, reduction = "umap")

saveRDS(Cluster_17, file = "TMS_merge_T_cells_Cluster_17_TSNE_UMAP.rds")
Cluster_17 <- readRDS(file = "TMS_merge_T_cells_Cluster_17_TSNE_UMAP.rds")

Cluster_17 <- FindNeighbors(Cluster_17, dims = 1:2)

#louvian
Cluster_17 <- FindClusters(Cluster_17, resolution = seq(0,3, by = 0.1))

clustree(Cluster_17, prefix = "RNA_snn_res.")
clustree(Cluster_17, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_17) <- Cluster_17$RNA_snn_res.2.7

Cluster_17 <- BuildClusterTree(Cluster_17, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
PlotClusterTree(object = Cluster_17)

Cluster_17.markers <- FindAllMarkers(Cluster_17)
m <- Cluster_17.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_17, features = m$gene)
write.csv(Cluster_17.markers, "TMS_T_cells_Cluster_17_merkers.csv")

Cluster_17 <- RenameIdents(Cluster_17, "1" = "CD8 T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD4 T",
                           "6" = "CD8 T", "7" = "CD8 T", "8" = "CD8 T", "9" = "DP T")

saveRDS(Cluster_17, file = "TMS_merge_T_cells_Cluster_17_clustered.rds")
Cluster_17 <- readRDS(file = "TMS_merge_T_cells_Cluster_17_clustered.rds")
