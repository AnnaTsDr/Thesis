Cluster_46 <- readRDS(file = "TMS_merge_T_cells_clusters_46.rds")

#Highly variable expressed genes
Cluster_46 <- FindVariableFeatures(Cluster_46, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_46),15)
plot1 <- VariableFeaturePlot(Cluster_46)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_46), file = "TMS_merge_T_cells_Cluster_46_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_46, file = "TMS_merge_T_cells_Cluster_46_VariableGenes.rds")
Cluster_46 <- readRDS(file = "TMS_merge_T_cells_Cluster_46_VariableGenes.rds")

all.genes <- rownames(Cluster_46)
Cluster_46 <- ScaleData(Cluster_46, features = all.genes)

saveRDS(Cluster_46, file = "TMS_merge_T_cells_Cluster_46_ScaledData.rds")
Cluster_46 <- readRDS(file = "TMS_merge_T_cells_Cluster_46_ScaledData.rds")

Cluster_46 <- RunPCA(Cluster_46, npcs = 45, ndims.print = 1:45, nfeatures.print = 5)
DimPlot(Cluster_46, reduction = "pca")

saveRDS(Cluster_46, file = "TMS_merge_T_cells_Cluster_46_PCA.rds")
Cluster_46 <- readRDS(file = "TMS_merge_T_cells_Cluster_46_PCA.rds")

Cluster_46 <- JackStraw(Cluster_46, num.replicate = 100, dims = 45)
Cluster_46 <- ScoreJackStraw(Cluster_46, dims = 1:45)
JackStrawPlot(Cluster_46, dims = 1:45)

ElbowPlot(Cluster_46, ndims = 45)

Cluster_46 <- RunTSNE(Cluster_46, dims = 1:45, perplexity = 10)
DimPlot(Cluster_46, reduction = "tsne")

Cluster_46 <- RunUMAP(Cluster_46, dims = 1:45)
DimPlot(Cluster_46, reduction = "umap")

saveRDS(Cluster_46, file = "TMS_merge_T_cells_Cluster_46_TSNE_UMAP.rds")
Cluster_46 <- readRDS(file = "TMS_merge_T_cells_Cluster_46_TSNE_UMAP.rds")

Cluster_46 <- FindNeighbors(Cluster_46, dims = 1:45)

#louvian
Cluster_46 <- FindClusters(Cluster_46, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_46, prefix = "RNA_snn_res.")
clustree(Cluster_46, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_46) <- Cluster_46$RNA_snn_res.1.2

Cluster_46 <- BuildClusterTree(Cluster_46, reorder.numeric = TRUE, reorder = TRUE, dims = 1:45)
PlotClusterTree(object = Cluster_46)

Cluster_46.markers <- FindAllMarkers(Cluster_46)
m <- Cluster_46.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_46, features = m$gene)
write.csv(Cluster_46.markers, "TMS_T_cells_Cluster_46_merkers.csv")

Cluster_46 <- RenameIdents(Cluster_46, "1" = "DN T", "2" = "DP T", "3" = "DP T", "4" = "DP T", "5" = "DP T", 
                           "6" = "DP T", "7" = "DP T", "8" = "DP T", "9" = "DP T")

saveRDS(Cluster_46, file = "TMS_merge_T_cells_Cluster_46_clustered.rds")
Cluster_46 <- readRDS(file = "TMS_merge_T_cells_Cluster_46_clustered.rds")

