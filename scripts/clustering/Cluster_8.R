Cluster_8 <- readRDS(file = "TMS_merge_T_cells_clusters_8.rds")

#Highly variable expressed genes
Cluster_8 <- FindVariableFeatures(Cluster_8, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_8),15)
plot1 <- VariableFeaturePlot(Cluster_8)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_8), file = "TMS_merge_T_cells_Cluster_8_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_8, file = "TMS_merge_T_cells_Cluster_8_VariableGenes.rds")
Cluster_8 <- readRDS(file = "TMS_merge_T_cells_Cluster_8_VariableGenes.rds")

all.genes <- rownames(Cluster_8)
Cluster_8 <- ScaleData(Cluster_8, features = all.genes)

saveRDS(Cluster_8, file = "TMS_merge_T_cells_Cluster_8_ScaledData.rds")
Cluster_8 <- readRDS(file = "TMS_merge_T_cells_Cluster_8_ScaledData.rds")

Cluster_8 <- RunPCA(Cluster_8, npcs = 13, ndims.print = 1:13, nfeatures.print = 5)
DimPlot(Cluster_8, reduction = "pca")

saveRDS(Cluster_8, file = "TMS_merge_T_cells_Cluster_8_PCA.rds")
Cluster_8 <- readRDS(file = "TMS_merge_T_cells_Cluster_8_PCA.rds")

Cluster_8 <- JackStraw(Cluster_8, num.replicate = 100, dims = 13)
Cluster_8 <- ScoreJackStraw(Cluster_8, dims = 1:13)
JackStrawPlot(Cluster_8, dims = 1:13)

ElbowPlot(Cluster_8, ndims = 13)

Cluster_8 <- RunTSNE(Cluster_8, dims = 1:13, perplexity = 10)
DimPlot(Cluster_8, reduction = "tsne")

Cluster_8 <- RunUMAP(Cluster_8, dims = 1:13)
DimPlot(Cluster_8, reduction = "umap")

saveRDS(Cluster_8, file = "TMS_merge_T_cells_Cluster_8_TSNE_UMAP.rds")
Cluster_8 <- readRDS(file = "TMS_merge_T_cells_Cluster_8_TSNE_UMAP.rds")

Cluster_8 <- FindNeighbors(Cluster_8, dims = 1:13)

#louvian
Cluster_8 <- FindClusters(Cluster_8, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_8, prefix = "RNA_snn_res.")
clustree(Cluster_8, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_8) <- Cluster_8$RNA_snn_res.1.2

Cluster_8 <- BuildClusterTree(Cluster_8, reorder.numeric = TRUE, reorder = TRUE, dims = 1:13)
PlotClusterTree(object = Cluster_8)

Cluster_8.markers <- FindAllMarkers(Cluster_8)
m <- Cluster_8.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_8, features = m$gene)
write.csv(Cluster_8.markers, "TMS_T_cells_Cluster_8_merkers.csv")

Cluster_8 <- RenameIdents(Cluster_8, "1" = "DN T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T")

saveRDS(Cluster_8, file = "TMS_merge_T_cells_Cluster_8_clustered.rds")
Cluster_8 <- readRDS(file = "TMS_merge_T_cells_Cluster_8_clustered.rds")

