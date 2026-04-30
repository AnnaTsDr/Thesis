Cluster_26 <- readRDS(file = "TMS_merge_T_cells_clusters_26.rds")

#Highly variable expressed genes
Cluster_26 <- FindVariableFeatures(Cluster_26, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0126, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_26),15)
plot1 <- VariableFeaturePlot(Cluster_26)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_26), file = "TMS_merge_T_cells_Cluster_26_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_26, file = "TMS_merge_T_cells_Cluster_26_VariableGenes.rds")
Cluster_26 <- readRDS(file = "TMS_merge_T_cells_Cluster_26_VariableGenes.rds")

all.genes <- rownames(Cluster_26)
Cluster_26 <- ScaleData(Cluster_26, features = all.genes)

saveRDS(Cluster_26, file = "TMS_merge_T_cells_Cluster_26_ScaledData.rds")
Cluster_26 <- readRDS(file = "TMS_merge_T_cells_Cluster_26_ScaledData.rds")

Cluster_26 <- RunPCA(Cluster_26, npcs = 3, ndims.print = 1:3, nfeatures.print = 5)
DimPlot(Cluster_26, reduction = "pca")

saveRDS(Cluster_26, file = "TMS_merge_T_cells_Cluster_26_PCA.rds")
Cluster_26 <- readRDS(file = "TMS_merge_T_cells_Cluster_26_PCA.rds")

Cluster_26 <- JackStraw(Cluster_26, num.replicate = 100, dims = 3)
Cluster_26 <- ScoreJackStraw(Cluster_26, dims = 1:3)
JackStrawPlot(Cluster_26, dims = 1:3)

ElbowPlot(Cluster_26, ndims = 3)

Cluster_26 <- RunTSNE(Cluster_26, dims = 1:3, perplexity = 10)
DimPlot(Cluster_26, reduction = "tsne")

Cluster_26 <- RunUMAP(Cluster_26, dims = 1:3)
DimPlot(Cluster_26, reduction = "umap")

saveRDS(Cluster_26, file = "TMS_merge_T_cells_Cluster_26_TSNE_UMAP.rds")
Cluster_26 <- readRDS(file = "TMS_merge_T_cells_Cluster_26_TSNE_UMAP.rds")

Cluster_26 <- FindNeighbors(Cluster_26, dims = 1:3)

#louvian
Cluster_26 <- FindClusters(Cluster_26, resolution = seq(0,3, by = 0.1))

clustree(Cluster_26, prefix = "RNA_snn_res.")
clustree(Cluster_26, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_26) <- Cluster_26$RNA_snn_res.3

Cluster_26 <- BuildClusterTree(Cluster_26, reorder.numeric = TRUE, reorder = TRUE, dims = 1:3)
PlotClusterTree(object = Cluster_26)

Cluster_26.markers <- FindAllMarkers(Cluster_26)
m <- Cluster_26.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_26, features = m$gene)
write.csv(Cluster_26.markers, "TMS_T_cells_Cluster_26_merkers.csv")

Cluster_26 <- RenameIdents(Cluster_26, "1" = "CD8 T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD8 T",
                           "6" = "CD8 T", "7" = "CD4 T", "8" = "CD8 T", "9" = "CD8 T", "10" = "CD8 T", "11" = "CD8 T", 
                           "12" = "CD8 T", "13" = "CD8 T", "14" = "CD8 T")

saveRDS(Cluster_26, file = "TMS_merge_T_cells_Cluster_26_clustered.rds")
Cluster_26 <- readRDS(file = "TMS_merge_T_cells_Cluster_26_clustered.rds")

