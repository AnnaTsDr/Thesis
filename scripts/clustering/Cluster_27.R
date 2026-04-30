Cluster_27 <- readRDS(file = "TMS_merge_T_cells_clusters_27.rds")

#Highly variable expressed genes
Cluster_27 <- FindVariableFeatures(Cluster_27, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.8, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_27),15)
plot1 <- VariableFeaturePlot(Cluster_27)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_27), file = "TMS_merge_T_cells_Cluster_27_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_27, file = "TMS_merge_T_cells_Cluster_27_VariableGenes.rds")
Cluster_27 <- readRDS(file = "TMS_merge_T_cells_Cluster_27_VariableGenes.rds")

all.genes <- rownames(Cluster_27)
Cluster_27 <- ScaleData(Cluster_27, features = all.genes)

saveRDS(Cluster_27, file = "TMS_merge_T_cells_Cluster_27_ScaledData.rds")
Cluster_27 <- readRDS(file = "TMS_merge_T_cells_Cluster_27_ScaledData.rds")

Cluster_27 <- RunPCA(Cluster_27, npcs = 8, ndims.print = 1:8, nfeatures.print = 5)
DimPlot(Cluster_27, reduction = "pca")

saveRDS(Cluster_27, file = "TMS_merge_T_cells_Cluster_27_PCA.rds")
Cluster_27 <- readRDS(file = "TMS_merge_T_cells_Cluster_27_PCA.rds")

Cluster_27 <- JackStraw(Cluster_27, num.replicate = 100, dims = 8)
Cluster_27 <- ScoreJackStraw(Cluster_27, dims = 1:8)
JackStrawPlot(Cluster_27, dims = 1:8)

ElbowPlot(Cluster_27, ndims = 8)

Cluster_27 <- RunTSNE(Cluster_27, dims = 1:8, perplexity = 10)
DimPlot(Cluster_27, reduction = "tsne")

Cluster_27 <- RunUMAP(Cluster_27, dims = 1:8)
DimPlot(Cluster_27, reduction = "umap")

saveRDS(Cluster_27, file = "TMS_merge_T_cells_Cluster_27_TSNE_UMAP.rds")
Cluster_27 <- readRDS(file = "TMS_merge_T_cells_Cluster_27_TSNE_UMAP.rds")

Cluster_27 <- FindNeighbors(Cluster_27, dims = 1:8)

#louvian
Cluster_27 <- FindClusters(Cluster_27, resolution = seq(0,2, by = 0.1))

clustree(Cluster_27, prefix = "RNA_snn_res.")
clustree(Cluster_27, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_27) <- Cluster_27$RNA_snn_res.2

Cluster_27 <- BuildClusterTree(Cluster_27, reorder.numeric = TRUE, reorder = TRUE, dims = 1:8)
PlotClusterTree(object = Cluster_27)

Cluster_27.markers <- FindAllMarkers(Cluster_27)
m <- Cluster_27.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_27, features = m$gene)
write.csv(Cluster_27.markers, "TMS_T_cells_Cluster_27_merkers.csv")

Cluster_27 <- RenameIdents(Cluster_27, "1" = "CD8 T", "2" = "DN T", "3" = "CD4 T", "4" = "CD8 T", "5" = "DN T", 
                           "6" = "CD8 T", "7" = "CD4 T", "8" = "CD8 T", "9" = "DP T", "10" = "CD4 T", "11" = "CD8 T")

saveRDS(Cluster_27, file = "TMS_merge_T_cells_Cluster_27_clustered.rds")
Cluster_27 <- readRDS(file = "TMS_merge_T_cells_Cluster_27_clustered.rds")

