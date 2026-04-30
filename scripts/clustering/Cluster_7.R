Cluster_7 <- readRDS(file = "TMS_merge_T_cells_clusters_7.rds")

#Highly variable expressed genes
Cluster_7 <- FindVariableFeatures(Cluster_7, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_7),15)
plot1 <- VariableFeaturePlot(Cluster_7)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_7), file = "TMS_merge_T_cells_Cluster_7_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_7, file = "TMS_merge_T_cells_Cluster_7_VariableGenes.rds")
Cluster_7 <- readRDS(file = "TMS_merge_T_cells_Cluster_7_VariableGenes.rds")

all.genes <- rownames(Cluster_7)
Cluster_7 <- ScaleData(Cluster_7, features = all.genes)

saveRDS(Cluster_7, file = "TMS_merge_T_cells_Cluster_7_ScaledData.rds")
Cluster_7 <- readRDS(file = "TMS_merge_T_cells_Cluster_7_ScaledData.rds")

Cluster_7 <- RunPCA(Cluster_7, npcs = 4, ndims.print = 1:4, nfeatures.print = 5)
DimPlot(Cluster_7, reduction = "pca")

saveRDS(Cluster_7, file = "TMS_merge_T_cells_Cluster_7_PCA.rds")
Cluster_7 <- readRDS(file = "TMS_merge_T_cells_Cluster_7_PCA.rds")

Cluster_7 <- JackStraw(Cluster_7, num.replicate = 100, dims = 4)
Cluster_7 <- ScoreJackStraw(Cluster_7, dims = 1:4)
JackStrawPlot(Cluster_7, dims = 1:4)

ElbowPlot(Cluster_7, ndims = 4)

Cluster_7 <- RunTSNE(Cluster_7, dims = 1:4, perplexity = 10)
DimPlot(Cluster_7, reduction = "tsne")

Cluster_7 <- RunUMAP(Cluster_7, dims = 1:4)
DimPlot(Cluster_7, reduction = "umap")

saveRDS(Cluster_7, file = "TMS_merge_T_cells_Cluster_7_TSNE_UMAP.rds")
Cluster_7 <- readRDS(file = "TMS_merge_T_cells_Cluster_7_TSNE_UMAP.rds")

Cluster_7 <- FindNeighbors(Cluster_7, dims = 1:4)

#louvian
Cluster_7 <- FindClusters(Cluster_7, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_7, prefix = "RNA_snn_res.")
clustree(Cluster_7, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_7) <- Cluster_7$RNA_snn_res.0.9

Cluster_7 <- BuildClusterTree(Cluster_7, reorder.numeric = TRUE, reorder = TRUE, dims = 1:4)
PlotClusterTree(object = Cluster_7)

Cluster_7.markers <- FindAllMarkers(Cluster_7)
m <- Cluster_7.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_7, features = m$gene)
write.csv(Cluster_7.markers, "TMS_T_cells_Cluster_7_merkers.csv")

Cluster_7 <- RenameIdents(Cluster_7, "1" = "CD8 T", "2" = "CD8 T", "3" = "DP T", "4" = "CD8 T", "5" = "CD8 T", 
                          "6" = "DN T")

saveRDS(Cluster_7, file = "TMS_merge_T_cells_Cluster_7_clustered.rds")
Cluster_7 <- readRDS(file = "TMS_merge_T_cells_Cluster_7_clustered.rds")
