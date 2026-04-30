Cluster_20 <- readRDS(file = "TMS_merge_T_cells_clusters_20.rds")

#Highly variable expressed genes
Cluster_20 <- FindVariableFeatures(Cluster_20, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_20),15)
plot1 <- VariableFeaturePlot(Cluster_20)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_20), file = "TMS_merge_T_cells_Cluster_20_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_20, file = "TMS_merge_T_cells_Cluster_20_VariableGenes.rds")
Cluster_20 <- readRDS(file = "TMS_merge_T_cells_Cluster_20_VariableGenes.rds")

all.genes <- rownames(Cluster_20)
Cluster_20 <- ScaleData(Cluster_20, features = all.genes)

saveRDS(Cluster_20, file = "TMS_merge_T_cells_Cluster_20_ScaledData.rds")
Cluster_20 <- readRDS(file = "TMS_merge_T_cells_Cluster_20_ScaledData.rds")

Cluster_20 <- RunPCA(Cluster_20, npcs = 9, ndims.print = 1:9, nfeatures.print = 5)
DimPlot(Cluster_20, reduction = "pca")

saveRDS(Cluster_20, file = "TMS_merge_T_cells_Cluster_20_PCA.rds")
Cluster_20 <- readRDS(file = "TMS_merge_T_cells_Cluster_20_PCA.rds")

Cluster_20 <- JackStraw(Cluster_20, num.replicate = 100, dims = 9)
Cluster_20 <- ScoreJackStraw(Cluster_20, dims = 1:9)
JackStrawPlot(Cluster_20, dims = 1:9)

ElbowPlot(Cluster_20, ndims = 9)

Cluster_20 <- RunTSNE(Cluster_20, dims = 1:9, perplexity = 10)
DimPlot(Cluster_20, reduction = "tsne")

Cluster_20 <- RunUMAP(Cluster_20, dims = 1:9)
DimPlot(Cluster_20, reduction = "umap")

saveRDS(Cluster_20, file = "TMS_merge_T_cells_Cluster_20_TSNE_UMAP.rds")
Cluster_20 <- readRDS(file = "TMS_merge_T_cells_Cluster_20_TSNE_UMAP.rds")

Cluster_20 <- FindNeighbors(Cluster_20, dims = 1:9)

#louvian
Cluster_20 <- FindClusters(Cluster_20, resolution = seq(0,2, by = 0.1))

clustree(Cluster_20, prefix = "RNA_snn_res.")
clustree(Cluster_20, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_20) <- Cluster_20$RNA_snn_res.2

Cluster_20 <- BuildClusterTree(Cluster_20, reorder.numeric = TRUE, reorder = TRUE, dims = 1:9)
PlotClusterTree(object = Cluster_20)

Cluster_20.markers <- FindAllMarkers(Cluster_20)
m <- Cluster_20.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_20, features = m$gene)
write.csv(Cluster_20.markers, "TMS_T_cells_Cluster_20_merkers.csv")

Cluster_20 <- RenameIdents(Cluster_20, "1" = "DN T", "2" = "DN T", "3" = "DN T", "4" = "CD8 T", "5" = "CD8 T", 
                           "6" = "DN T", "7" = "DN T", "8" = "DN T", "9" = "CD8 T", "10" = "CD8 T")

saveRDS(Cluster_20, file = "TMS_merge_T_cells_Cluster_20_clustered.rds")
Cluster_20 <- readRDS(file = "TMS_merge_T_cells_Cluster_20_clustered.rds")

