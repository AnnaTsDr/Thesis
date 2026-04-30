Cluster_25 <- readRDS(file = "TMS_merge_T_cells_clusters_25.rds")

#Highly variable expressed genes
Cluster_25 <- FindVariableFeatures(Cluster_25, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_25),15)
plot1 <- VariableFeaturePlot(Cluster_25)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_25), file = "TMS_merge_T_cells_Cluster_25_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_25, file = "TMS_merge_T_cells_Cluster_25_VariableGenes.rds")
Cluster_25 <- readRDS(file = "TMS_merge_T_cells_Cluster_25_VariableGenes.rds")

all.genes <- rownames(Cluster_25)
Cluster_25 <- ScaleData(Cluster_25, features = all.genes)

saveRDS(Cluster_25, file = "TMS_merge_T_cells_Cluster_25_ScaledData.rds")
Cluster_25 <- readRDS(file = "TMS_merge_T_cells_Cluster_25_ScaledData.rds")

Cluster_25 <- RunPCA(Cluster_25, npcs = 9, ndims.print = 1:9, nfeatures.print = 5)
DimPlot(Cluster_25, reduction = "pca")

saveRDS(Cluster_25, file = "TMS_merge_T_cells_Cluster_25_PCA.rds")
Cluster_25 <- readRDS(file = "TMS_merge_T_cells_Cluster_25_PCA.rds")

Cluster_25 <- JackStraw(Cluster_25, num.replicate = 100, dims = 9)
Cluster_25 <- ScoreJackStraw(Cluster_25, dims = 1:9)
JackStrawPlot(Cluster_25, dims = 1:9)

ElbowPlot(Cluster_25, ndims = 9)

Cluster_25 <- RunTSNE(Cluster_25, dims = 1:9, perplexity = 10)
DimPlot(Cluster_25, reduction = "tsne")

Cluster_25 <- RunUMAP(Cluster_25, dims = 1:9)
DimPlot(Cluster_25, reduction = "umap")

saveRDS(Cluster_25, file = "TMS_merge_T_cells_Cluster_25_TSNE_UMAP.rds")
Cluster_25 <- readRDS(file = "TMS_merge_T_cells_Cluster_25_TSNE_UMAP.rds")

Cluster_25 <- FindNeighbors(Cluster_25, dims = 1:9)

#louvian
Cluster_25 <- FindClusters(Cluster_25, resolution = seq(0,3, by = 0.1))

clustree(Cluster_25, prefix = "RNA_snn_res.")
clustree(Cluster_25, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_25) <- Cluster_25$RNA_snn_res.3

Cluster_25 <- BuildClusterTree(Cluster_25, reorder.numeric = TRUE, reorder = TRUE, dims = 1:9)
PlotClusterTree(object = Cluster_25)

Cluster_25.markers <- FindAllMarkers(Cluster_25)
m <- Cluster_25.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_25, features = m$gene)
write.csv(Cluster_25.markers, "TMS_T_cells_Cluster_25_merkers.csv")

Cluster_25 <- RenameIdents(Cluster_25, "1" = "CD8 T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD8 T", 
                           "6" = "CD8 T", "7" ="CD8 T", "8" = "CD4 T", "9" = "CD8 T", "10" = "CD8 T", "11" = "CD8 T", 
                           "12" = "CD8 T", "13" = "CD8 T", "14" = "CD8 T", "15" = "CD8 T", "16" = "CD8 T")

saveRDS(Cluster_25, file = "TMS_merge_T_cells_Cluster_25_clustered.rds")
Cluster_25 <- readRDS(file = "TMS_merge_T_cells_Cluster_25_clustered.rds")

