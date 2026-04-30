Cluster_40 <- readRDS(file = "TMS_merge_T_cells_clusters_40.rds")

#Highly variable expressed genes
Cluster_40 <- FindVariableFeatures(Cluster_40, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_40),15)
plot1 <- VariableFeaturePlot(Cluster_40)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_40), file = "TMS_merge_T_cells_Cluster_40_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_40, file = "TMS_merge_T_cells_Cluster_40_VariableGenes.rds")
Cluster_40 <- readRDS(file = "TMS_merge_T_cells_Cluster_40_VariableGenes.rds")

all.genes <- rownames(Cluster_40)
Cluster_40 <- ScaleData(Cluster_40, features = all.genes)

saveRDS(Cluster_40, file = "TMS_merge_T_cells_Cluster_40_ScaledData.rds")
Cluster_40 <- readRDS(file = "TMS_merge_T_cells_Cluster_40_ScaledData.rds")

Cluster_40 <- RunPCA(Cluster_40, npcs = 14, ndims.print = 1:14, nfeatures.print = 5)
DimPlot(Cluster_40, reduction = "pca")

saveRDS(Cluster_40, file = "TMS_merge_T_cells_Cluster_40_PCA.rds")
Cluster_40 <- readRDS(file = "TMS_merge_T_cells_Cluster_40_PCA.rds")

Cluster_40 <- JackStraw(Cluster_40, num.replicate = 100, dims = 14)
Cluster_40 <- ScoreJackStraw(Cluster_40, dims = 1:14)
JackStrawPlot(Cluster_40, dims = 1:14)

ElbowPlot(Cluster_40, ndims = 14)

Cluster_40 <- RunTSNE(Cluster_40, dims = 1:14, perplexity = 10)
DimPlot(Cluster_40, reduction = "tsne")

Cluster_40 <- RunUMAP(Cluster_40, dims = 1:14)
DimPlot(Cluster_40, reduction = "umap")

saveRDS(Cluster_40, file = "TMS_merge_T_cells_Cluster_40_TSNE_UMAP.rds")
Cluster_40 <- readRDS(file = "TMS_merge_T_cells_Cluster_40_TSNE_UMAP.rds")

Cluster_40 <- FindNeighbors(Cluster_40, dims = 1:14)

#louvian
Cluster_40 <- FindClusters(Cluster_40, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_40, prefix = "RNA_snn_res.")
clustree(Cluster_40, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_40) <- Cluster_40$RNA_snn_res.1.2

Cluster_40 <- BuildClusterTree(Cluster_40, reorder.numeric = TRUE, reorder = TRUE, dims = 1:14)
PlotClusterTree(object = Cluster_40)

Cluster_40.markers <- FindAllMarkers(Cluster_40)
m <- Cluster_40.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_40, features = m$gene)
write.csv(Cluster_40.markers, "TMS_T_cells_Cluster_40_merkers.csv")

Cluster_40 <- RenameIdents(Cluster_40, "1" = "CD4 T", "2" = "CD4 T", "3" = "CD4 T", "4" = "CD4 T", "5" = "CD4 T", 
                           "6" = "CD4 T", "7" = "CD4 T", "8" = "CD4 T")

saveRDS(Cluster_40, file = "TMS_merge_T_cells_Cluster_40_clustered.rds")
Cluster_40 <- readRDS(file = "TMS_merge_T_cells_Cluster_40_clustered.rds")
