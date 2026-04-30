Cluster_21 <- readRDS(file = "TMS_merge_T_cells_clusters_21.rds")

#Highly variable expressed genes
Cluster_21 <- FindVariableFeatures(Cluster_21, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_21),15)
plot1 <- VariableFeaturePlot(Cluster_21)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_21), file = "TMS_merge_T_cells_Cluster_21_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_21, file = "TMS_merge_T_cells_Cluster_21_VariableGenes.rds")
Cluster_21 <- readRDS(file = "TMS_merge_T_cells_Cluster_21_VariableGenes.rds")

all.genes <- rownames(Cluster_21)
Cluster_21 <- ScaleData(Cluster_21, features = all.genes)

saveRDS(Cluster_21, file = "TMS_merge_T_cells_Cluster_21_ScaledData.rds")
Cluster_21 <- readRDS(file = "TMS_merge_T_cells_Cluster_21_ScaledData.rds")

Cluster_21 <- RunPCA(Cluster_21, npcs = 5, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(Cluster_21, reduction = "pca")

saveRDS(Cluster_21, file = "TMS_merge_T_cells_Cluster_21_PCA.rds")
Cluster_21 <- readRDS(file = "TMS_merge_T_cells_Cluster_21_PCA.rds")

Cluster_21 <- JackStraw(Cluster_21, num.replicate = 100, dims = 5)
Cluster_21 <- ScoreJackStraw(Cluster_21, dims = 1:5)
JackStrawPlot(Cluster_21, dims = 1:5)

ElbowPlot(Cluster_21, ndims = 5)

Cluster_21 <- RunTSNE(Cluster_21, dims = 1:5, perplexity = 10)
DimPlot(Cluster_21, reduction = "tsne")

Cluster_21 <- RunUMAP(Cluster_21, dims = 1:5)
DimPlot(Cluster_21, reduction = "umap")

saveRDS(Cluster_21, file = "TMS_merge_T_cells_Cluster_21_TSNE_UMAP.rds")
Cluster_21 <- readRDS(file = "TMS_merge_T_cells_Cluster_21_TSNE_UMAP.rds")

Cluster_21 <- FindNeighbors(Cluster_21, dims = 1:5)

#louvian
Cluster_21 <- FindClusters(Cluster_21, resolution = seq(0,3, by = 0.1))

clustree(Cluster_21, prefix = "RNA_snn_res.")
clustree(Cluster_21, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_21) <- Cluster_21$RNA_snn_res.3

Cluster_21 <- BuildClusterTree(Cluster_21, reorder.numeric = TRUE, reorder = TRUE, dims = 1:5)
PlotClusterTree(object = Cluster_21)

Cluster_21.markers <- FindAllMarkers(Cluster_21)
m <- Cluster_21.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_21, features = m$gene)
write.csv(Cluster_21.markers, "TMS_T_cells_Cluster_21_merkers.csv")

Cluster_21 <- RenameIdents(Cluster_21, "1" = "DN T", "2" = "CD8 T", "3" = "CD8 T", "4" = "DN T", "5" = "CD8 T", 
                           "6" = "CD8 T", "7" = "CD8 T", "8" = "CD8 T", "9" = "CD8 T", "10" = "CD8 T", "11" = "CD8 T", 
                           "12" = "CD8 T", "13" = "CD8 T", "14" = "CD8 T")

saveRDS(Cluster_21, file = "TMS_merge_T_cells_Cluster_21_clustered.rds")
Cluster_21 <- readRDS(file = "TMS_merge_T_cells_Cluster_21_clustered.rds")

