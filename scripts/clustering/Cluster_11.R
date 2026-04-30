Cluster_11 <- readRDS(file = "TMS_merge_T_cells_clusters_11.rds")

#Highly variable expressed genes
Cluster_11 <- FindVariableFeatures(Cluster_11, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_11),15)
plot1 <- VariableFeaturePlot(Cluster_11)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_11), file = "TMS_merge_T_cells_Cluster_11_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_11, file = "TMS_merge_T_cells_Cluster_11_VariableGenes.rds")
Cluster_11 <- readRDS(file = "TMS_merge_T_cells_Cluster_11_VariableGenes.rds")

all.genes <- rownames(Cluster_11)
Cluster_11 <- ScaleData(Cluster_11, features = all.genes)

saveRDS(Cluster_11, file = "TMS_merge_T_cells_Cluster_11_ScaledData.rds")
Cluster_11 <- readRDS(file = "TMS_merge_T_cells_Cluster_11_ScaledData.rds")

Cluster_11 <- RunPCA(Cluster_11, npcs = 7, ndims.print = 1:7, nfeatures.print = 5)
DimPlot(Cluster_11, reduction = "pca")

saveRDS(Cluster_11, file = "TMS_merge_T_cells_Cluster_11_PCA.rds")
Cluster_11 <- readRDS(file = "TMS_merge_T_cells_Cluster_11_PCA.rds")

Cluster_11 <- JackStraw(Cluster_11, num.replicate = 100, dims = 7)
Cluster_11 <- ScoreJackStraw(Cluster_11, dims = 1:7)
JackStrawPlot(Cluster_11, dims = 1:7)

ElbowPlot(Cluster_11, ndims = 7)

Cluster_11 <- RunTSNE(Cluster_11, dims = 1:7, perplexity = 10)
DimPlot(Cluster_11, reduction = "tsne")

Cluster_11 <- RunUMAP(Cluster_11, dims = 1:7)
DimPlot(Cluster_11, reduction = "umap")

saveRDS(Cluster_11, file = "TMS_merge_T_cells_Cluster_11_TSNE_UMAP.rds")
Cluster_11 <- readRDS(file = "TMS_merge_T_cells_Cluster_11_TSNE_UMAP.rds")

Cluster_11 <- FindNeighbors(Cluster_11, dims = 1:7)

#louvian
Cluster_11 <- FindClusters(Cluster_11, resolution = seq(0,2, by = 0.1))

clustree(Cluster_11, prefix = "RNA_snn_res.")
clustree(Cluster_11, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_11) <- Cluster_11$RNA_snn_res.2

Cluster_11 <- BuildClusterTree(Cluster_11, reorder.numeric = TRUE, reorder = TRUE, dims = 1:7)
PlotClusterTree(object = Cluster_11)

Cluster_11.markers <- FindAllMarkers(Cluster_11)
m <- Cluster_11.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_11, features = m$gene)
write.csv(Cluster_11.markers, "TMS_T_cells_Cluster_11_merkers.csv")

Cluster_11 <- RenameIdents(Cluster_11, "1" = "CD8 T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD8 T",
                           "6" = "CD8 T", "7" = "DN T", "8" = "CD4 T")

saveRDS(Cluster_11, file = "TMS_merge_T_cells_Cluster_11_clustered.rds")
Cluster_11 <- readRDS(file = "TMS_merge_T_cells_Cluster_11_clustered.rds")
