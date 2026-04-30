Cluster_13 <- readRDS(file = "TMS_merge_T_cells_clusters_13.rds")

#Highly variable expressed genes
Cluster_13 <- FindVariableFeatures(Cluster_13, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_13),15)
plot1 <- VariableFeaturePlot(Cluster_13)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_13), file = "TMS_merge_T_cells_Cluster_13_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_13, file = "TMS_merge_T_cells_Cluster_13_VariableGenes.rds")
Cluster_13 <- readRDS(file = "TMS_merge_T_cells_Cluster_13_VariableGenes.rds")

all.genes <- rownames(Cluster_13)
Cluster_13 <- ScaleData(Cluster_13, features = all.genes)

saveRDS(Cluster_13, file = "TMS_merge_T_cells_Cluster_13_ScaledData.rds")
Cluster_13 <- readRDS(file = "TMS_merge_T_cells_Cluster_13_ScaledData.rds")

Cluster_13 <- RunPCA(Cluster_13, npcs = 13, ndims.print = 1:13, nfeatures.print = 5)
DimPlot(Cluster_13, reduction = "pca")

saveRDS(Cluster_13, file = "TMS_merge_T_cells_Cluster_13_PCA.rds")
Cluster_13 <- readRDS(file = "TMS_merge_T_cells_Cluster_13_PCA.rds")

Cluster_13 <- JackStraw(Cluster_13, num.replicate = 100, dims = 13)
Cluster_13 <- ScoreJackStraw(Cluster_13, dims = 1:13)
JackStrawPlot(Cluster_13, dims = 1:13)

ElbowPlot(Cluster_13, ndims = 13)

Cluster_13 <- RunTSNE(Cluster_13, dims = 1:13, perplexity = 10)
DimPlot(Cluster_13, reduction = "tsne")

Cluster_13 <- RunUMAP(Cluster_13, dims = 1:13)
DimPlot(Cluster_13, reduction = "umap")

saveRDS(Cluster_13, file = "TMS_merge_T_cells_Cluster_13_TSNE_UMAP.rds")
Cluster_13 <- readRDS(file = "TMS_merge_T_cells_Cluster_13_TSNE_UMAP.rds")

Cluster_13 <- FindNeighbors(Cluster_13, dims = 1:13)

#louvian
Cluster_13 <- FindClusters(Cluster_13, resolution = seq(0,3, by = 0.1))

clustree(Cluster_13, prefix = "RNA_snn_res.")
clustree(Cluster_13, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_13) <- Cluster_13$RNA_snn_res.3

Cluster_13 <- BuildClusterTree(Cluster_13, reorder.numeric = TRUE, reorder = TRUE, dims = 1:13)
PlotClusterTree(object = Cluster_13)

Cluster_13.markers <- FindAllMarkers(Cluster_13)
m <- Cluster_13.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_13, features = m$gene)
write.csv(Cluster_13.markers, "TMS_T_cells_Cluster_13_merkers.csv")

Cluster_13 <- RenameIdents(Cluster_13, "1" = "DP T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD8 T", 
                           "6" = "CD8 T", "7" = "CD8 T", "8" = "CD8 T", "9" = "CD8 T", "10" = "CD8 T", "11" = "CD8 T",
                           "12" = "DP T", "13" = "DP T", "14" = "CD8 T", "15" = "DP T", "16" = "DP T", "17" = "DP T",
                           "18" = "CD4 T")

saveRDS(Cluster_13, file = "TMS_merge_T_cells_Cluster_13_clustered.rds")
Cluster_13 <- readRDS(file = "TMS_merge_T_cells_Cluster_13_clustered.rds")

