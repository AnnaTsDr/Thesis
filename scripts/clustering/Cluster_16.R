Cluster_16 <- readRDS(file = "TMS_merge_T_cells_clusters_16.rds")

#Highly variable expressed genes
Cluster_16 <- FindVariableFeatures(Cluster_16, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_16),15)
plot1 <- VariableFeaturePlot(Cluster_16)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_16), file = "TMS_merge_T_cells_Cluster_16_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_16, file = "TMS_merge_T_cells_Cluster_16_VariableGenes.rds")
Cluster_16 <- readRDS(file = "TMS_merge_T_cells_Cluster_16_VariableGenes.rds")

all.genes <- rownames(Cluster_16)
Cluster_16 <- ScaleData(Cluster_16, features = all.genes)

saveRDS(Cluster_16, file = "TMS_merge_T_cells_Cluster_16_ScaledData.rds")
Cluster_16 <- readRDS(file = "TMS_merge_T_cells_Cluster_16_ScaledData.rds")

Cluster_16 <- RunPCA(Cluster_16, npcs = 7, ndims.print = 1:7, nfeatures.print = 5)
DimPlot(Cluster_16, reduction = "pca")

saveRDS(Cluster_16, file = "TMS_merge_T_cells_Cluster_16_PCA.rds")
Cluster_16 <- readRDS(file = "TMS_merge_T_cells_Cluster_16_PCA.rds")

Cluster_16 <- JackStraw(Cluster_16, num.replicate = 100, dims = 7)
Cluster_16 <- ScoreJackStraw(Cluster_16, dims = 1:7)
JackStrawPlot(Cluster_16, dims = 1:7)

ElbowPlot(Cluster_16, ndims = 7)

Cluster_16 <- RunTSNE(Cluster_16, dims = 1:7, perplexity = 10)
DimPlot(Cluster_16, reduction = "tsne")

Cluster_16 <- RunUMAP(Cluster_16, dims = 1:7)
DimPlot(Cluster_16, reduction = "umap")

saveRDS(Cluster_16, file = "TMS_merge_T_cells_Cluster_16_TSNE_UMAP.rds")
Cluster_16 <- readRDS(file = "TMS_merge_T_cells_Cluster_16_TSNE_UMAP.rds")

Cluster_16 <- FindNeighbors(Cluster_16, dims = 1:7)

#louvian
Cluster_16 <- FindClusters(Cluster_16, resolution = seq(0,3, by = 0.1))

clustree(Cluster_16, prefix = "RNA_snn_res.")
clustree(Cluster_16, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_16) <- Cluster_16$RNA_snn_res.3

Cluster_16 <- BuildClusterTree(Cluster_16, reorder.numeric = TRUE, reorder = TRUE, dims = 1:7)
PlotClusterTree(object = Cluster_16)

Cluster_16.markers <- FindAllMarkers(Cluster_16)
m <- Cluster_16.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_16, features = m$gene)
write.csv(Cluster_16.markers, "TMS_T_cells_Cluster_16_merkers.csv")

Cluster_16 <- RenameIdents(Cluster_16, "1" = "DP T", "2" = "DP T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD8 T",
                           "6" = "CD8 T", "7" = "CD8 T", "8" = "DN T", "9" = "CD4 T", "10" = "CD8 T", "11" = "DN T")

saveRDS(Cluster_16, file = "TMS_merge_T_cells_Cluster_16_clustered.rds")
Cluster_16 <- readRDS(file = "TMS_merge_T_cells_Cluster_16_clustered.rds")
