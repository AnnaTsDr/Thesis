Cluster_29 <- readRDS(file = "TMS_merge_T_cells_clusters_29.rds")

#Highly variable expressed genes
Cluster_29 <- FindVariableFeatures(Cluster_29, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_29),15)
plot1 <- VariableFeaturePlot(Cluster_29)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_29), file = "TMS_merge_T_cells_Cluster_29_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_29, file = "TMS_merge_T_cells_Cluster_29_VariableGenes.rds")
Cluster_29 <- readRDS(file = "TMS_merge_T_cells_Cluster_29_VariableGenes.rds")

all.genes <- rownames(Cluster_29)
Cluster_29 <- ScaleData(Cluster_29, features = all.genes)

saveRDS(Cluster_29, file = "TMS_merge_T_cells_Cluster_29_ScaledData.rds")
Cluster_29 <- readRDS(file = "TMS_merge_T_cells_Cluster_29_ScaledData.rds")

Cluster_29 <- RunPCA(Cluster_29, npcs = 10, ndims.print = 1:10, nfeatures.print = 5)
DimPlot(Cluster_29, reduction = "pca")

saveRDS(Cluster_29, file = "TMS_merge_T_cells_Cluster_29_PCA.rds")
Cluster_29 <- readRDS(file = "TMS_merge_T_cells_Cluster_29_PCA.rds")

Cluster_29 <- JackStraw(Cluster_29, num.replicate = 100, dims = 10)
Cluster_29 <- ScoreJackStraw(Cluster_29, dims = 1:10)
JackStrawPlot(Cluster_29, dims = 1:10)

ElbowPlot(Cluster_29, ndims = 10)

Cluster_29 <- RunTSNE(Cluster_29, dims = 1:10, perplexity = 10)
DimPlot(Cluster_29, reduction = "tsne")

Cluster_29 <- RunUMAP(Cluster_29, dims = 1:10)
DimPlot(Cluster_29, reduction = "umap")

saveRDS(Cluster_29, file = "TMS_merge_T_cells_Cluster_29_TSNE_UMAP.rds")
Cluster_29 <- readRDS(file = "TMS_merge_T_cells_Cluster_29_TSNE_UMAP.rds")

Cluster_29 <- FindNeighbors(Cluster_29, dims = 1:10)

#louvian
Cluster_29 <- FindClusters(Cluster_29, resolution = seq(0,2, by = 0.1))

clustree(Cluster_29, prefix = "RNA_snn_res.")
clustree(Cluster_29, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_29) <- Cluster_29$RNA_snn_res.2

Cluster_29 <- BuildClusterTree(Cluster_29, reorder.numeric = TRUE, reorder = TRUE, dims = 1:10)
PlotClusterTree(object = Cluster_29)

Cluster_29.markers <- FindAllMarkers(Cluster_29)
m <- Cluster_29.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_29, features = m$gene)
write.csv(Cluster_29.markers, "TMS_T_cells_Cluster_29_merkers.csv")

Cluster_29 <- RenameIdents(Cluster_29, "1" = "CD8 T", "2" = "DP T", "3" = "CD8 T", "4" = "DP T", "5" = "DP T", 
                           "6" = "CD8 T", "7" = "DP T", "8" = "DP T", "9" = "DP T", "10" = "DP T", "11" = "DP T")

saveRDS(Cluster_29, file = "TMS_merge_T_cells_Cluster_29_clustered.rds")
Cluster_29 <- readRDS(file = "TMS_merge_T_cells_Cluster_29_clustered.rds")

