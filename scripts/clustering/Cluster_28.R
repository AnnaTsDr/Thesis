Cluster_28 <- readRDS(file = "TMS_merge_T_cells_clusters_28.rds")

#Highly variable expressed genes
Cluster_28 <- FindVariableFeatures(Cluster_28, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_28),15)
plot1 <- VariableFeaturePlot(Cluster_28)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_28), file = "TMS_merge_T_cells_Cluster_28_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_28, file = "TMS_merge_T_cells_Cluster_28_VariableGenes.rds")
Cluster_28 <- readRDS(file = "TMS_merge_T_cells_Cluster_28_VariableGenes.rds")

all.genes <- rownames(Cluster_28)
Cluster_28 <- ScaleData(Cluster_28, features = all.genes)

saveRDS(Cluster_28, file = "TMS_merge_T_cells_Cluster_28_ScaledData.rds")
Cluster_28 <- readRDS(file = "TMS_merge_T_cells_Cluster_28_ScaledData.rds")

Cluster_28 <- RunPCA(Cluster_28, npcs = 2, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_28, reduction = "pca")

saveRDS(Cluster_28, file = "TMS_merge_T_cells_Cluster_28_PCA.rds")
Cluster_28 <- readRDS(file = "TMS_merge_T_cells_Cluster_28_PCA.rds")

Cluster_28 <- JackStraw(Cluster_28, num.replicate = 100, dims = 2)
Cluster_28 <- ScoreJackStraw(Cluster_28, dims = 1:2)
JackStrawPlot(Cluster_28, dims = 1:2)

ElbowPlot(Cluster_28, ndims = 2)

Cluster_28 <- RunTSNE(Cluster_28, dims = 1:2, perplexity = 10)
DimPlot(Cluster_28, reduction = "tsne")

Cluster_28 <- RunUMAP(Cluster_28, dims = 1:2)
DimPlot(Cluster_28, reduction = "umap")

saveRDS(Cluster_28, file = "TMS_merge_T_cells_Cluster_28_TSNE_UMAP.rds")
Cluster_28 <- readRDS(file = "TMS_merge_T_cells_Cluster_28_TSNE_UMAP.rds")

Cluster_28 <- FindNeighbors(Cluster_28, dims = 1:2)

#louvian
Cluster_28 <- FindClusters(Cluster_28, resolution = seq(0,3, by = 0.1))

clustree(Cluster_28, prefix = "RNA_snn_res.")
clustree(Cluster_28, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_28) <- Cluster_28$RNA_snn_res.3

Cluster_28 <- BuildClusterTree(Cluster_28, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
PlotClusterTree(object = Cluster_28)

Cluster_28.markers <- FindAllMarkers(Cluster_28)
m <- Cluster_28.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_28, features = m$gene)
write.csv(Cluster_28.markers, "TMS_T_cells_Cluster_28_merkers.csv")

Cluster_28 <- RenameIdents(Cluster_28, "1" = "DN T", "2" = "DP T", "3" = "DN T", "4" = "CD8 T", "5" = "DP T",
                           "6" = "DP T", "7" = "DP T", "8" = "DP T", "9" = "DP T", "10" = "DP T")

saveRDS(Cluster_28, file = "TMS_merge_T_cells_Cluster_28_clustered.rds")
Cluster_28 <- readRDS(file = "TMS_merge_T_cells_Cluster_28_clustered.rds")

