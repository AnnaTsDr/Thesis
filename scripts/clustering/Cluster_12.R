Cluster_12 <- readRDS(file = "TMS_merge_T_cells_clusters_12.rds")

#Highly variable expressed genes
Cluster_12 <- FindVariableFeatures(Cluster_12, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_12),15)
plot1 <- VariableFeaturePlot(Cluster_12)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_12), file = "TMS_merge_T_cells_Cluster_12_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_12, file = "TMS_merge_T_cells_Cluster_12_VariableGenes.rds")
Cluster_12 <- readRDS(file = "TMS_merge_T_cells_Cluster_12_VariableGenes.rds")

all.genes <- rownames(Cluster_12)
Cluster_12 <- ScaleData(Cluster_12, features = all.genes)

saveRDS(Cluster_12, file = "TMS_merge_T_cells_Cluster_12_ScaledData.rds")
Cluster_12 <- readRDS(file = "TMS_merge_T_cells_Cluster_12_ScaledData.rds")

Cluster_12 <- RunPCA(Cluster_12, npcs = 10, ndims.print = 1:10, nfeatures.print = 5)
DimPlot(Cluster_12, reduction = "pca")

saveRDS(Cluster_12, file = "TMS_merge_T_cells_Cluster_12_PCA.rds")
Cluster_12 <- readRDS(file = "TMS_merge_T_cells_Cluster_12_PCA.rds")

Cluster_12 <- JackStraw(Cluster_12, num.replicate = 100, dims = 10)
Cluster_12 <- ScoreJackStraw(Cluster_12, dims = 1:10)
JackStrawPlot(Cluster_12, dims = 1:10)

ElbowPlot(Cluster_12, ndims = 10)

Cluster_12 <- RunTSNE(Cluster_12, dims = 1:10, perplexity = 10)
DimPlot(Cluster_12, reduction = "tsne")

Cluster_12 <- RunUMAP(Cluster_12, dims = 1:10)
DimPlot(Cluster_12, reduction = "umap")

saveRDS(Cluster_12, file = "TMS_merge_T_cells_Cluster_12_TSNE_UMAP.rds")
Cluster_12 <- readRDS(file = "TMS_merge_T_cells_Cluster_12_TSNE_UMAP.rds")

Cluster_12 <- FindNeighbors(Cluster_12, dims = 1:10)

#louvian
Cluster_12 <- FindClusters(Cluster_12, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_12, prefix = "RNA_snn_res.")
clustree(Cluster_12, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_12) <- Cluster_12$RNA_snn_res.1.2

Cluster_12 <- BuildClusterTree(Cluster_12, reorder.numeric = TRUE, reorder = TRUE, dims = 1:10)
PlotClusterTree(object = Cluster_12)

Cluster_12.markers <- FindAllMarkers(Cluster_12)
m <- Cluster_12.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_12, features = m$gene)
write.csv(Cluster_12.markers, "TMS_T_cells_Cluster_12_merkers.csv")

Cluster_12 <- RenameIdents(Cluster_12, "1" = "DP T", "2" = "DP T", "3" = "DP T", "4" = "DP T", "5" = "CD8 T", 
                           "6" = "DP T", "7" = "DP T", "8" = "DP T")

saveRDS(Cluster_12, file = "TMS_merge_T_cells_Cluster_12_clustered.rds")
Cluster_12 <- readRDS(file = "TMS_merge_T_cells_Cluster_12_clustered.rds")

