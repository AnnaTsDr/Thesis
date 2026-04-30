Cluster_39 <- readRDS(file = "TMS_merge_T_cells_clusters_39.rds")

#Highly variable expressed genes
Cluster_39 <- FindVariableFeatures(Cluster_39, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_39),15)
plot1 <- VariableFeaturePlot(Cluster_39)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_39), file = "TMS_merge_T_cells_Cluster_39_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_39, file = "TMS_merge_T_cells_Cluster_39_VariableGenes.rds")
Cluster_39 <- readRDS(file = "TMS_merge_T_cells_Cluster_39_VariableGenes.rds")

all.genes <- rownames(Cluster_39)
Cluster_39 <- ScaleData(Cluster_39, features = all.genes)

saveRDS(Cluster_39, file = "TMS_merge_T_cells_Cluster_39_ScaledData.rds")
Cluster_39 <- readRDS(file = "TMS_merge_T_cells_Cluster_39_ScaledData.rds")

Cluster_39 <- RunPCA(Cluster_39, npcs = 4, ndims.print = 1:4, nfeatures.print = 5)
DimPlot(Cluster_39, reduction = "pca")

saveRDS(Cluster_39, file = "TMS_merge_T_cells_Cluster_39_PCA.rds")
Cluster_39 <- readRDS(file = "TMS_merge_T_cells_Cluster_39_PCA.rds")

Cluster_39 <- JackStraw(Cluster_39, num.replicate = 100, dims = 4)
Cluster_39 <- ScoreJackStraw(Cluster_39, dims = 1:4)
JackStrawPlot(Cluster_39, dims = 1:4)

ElbowPlot(Cluster_39, ndims = 4)

Cluster_39 <- RunTSNE(Cluster_39, dims = 1:4, perplexity = 10)
DimPlot(Cluster_39, reduction = "tsne")

Cluster_39 <- RunUMAP(Cluster_39, dims = 1:4)
DimPlot(Cluster_39, reduction = "umap")

saveRDS(Cluster_39, file = "TMS_merge_T_cells_Cluster_39_TSNE_UMAP.rds")
Cluster_39 <- readRDS(file = "TMS_merge_T_cells_Cluster_39_TSNE_UMAP.rds")

Cluster_39 <- FindNeighbors(Cluster_39, dims = 1:4)

#louvian
Cluster_39 <- FindClusters(Cluster_39, resolution = seq(0,3, by = 0.1))

clustree(Cluster_39, prefix = "RNA_snn_res.")
clustree(Cluster_39, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_39) <- Cluster_39$RNA_snn_res.3

Cluster_39 <- BuildClusterTree(Cluster_39, reorder.numeric = TRUE, reorder = TRUE, dims = 1:4)
PlotClusterTree(object = Cluster_39)

Cluster_39.markers <- FindAllMarkers(Cluster_39)
m <- Cluster_39.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_39, features = m$gene)
write.csv(Cluster_39.markers, "TMS_T_cells_Cluster_39_merkers.csv")

Cluster_39 <- RenameIdents(Cluster_39, "1" = "DP T", "2" = "CD4 T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD4 T", 
                           "6" = "CD8 T", "7" = "CD8 T", "8" = "CD8 T", "9" = "CD8 T", "10" = "CD8 T", "11" = "CD8 T", 
                           "12" = "CD8 T", "13" = "CD8 T", "14" = "CD8 T", "15" = "CD8 T")

saveRDS(Cluster_39, file = "TMS_merge_T_cells_Cluster_39_clustered.rds")
Cluster_39 <- readRDS(file = "TMS_merge_T_cells_Cluster_39_clustered.rds")
