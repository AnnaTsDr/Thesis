Cluster_37 <- readRDS(file = "TMS_merge_T_cells_clusters_37.rds")

#Highly variable expressed genes
Cluster_37 <- FindVariableFeatures(Cluster_37, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_37),15)
plot1 <- VariableFeaturePlot(Cluster_37)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_37), file = "TMS_merge_T_cells_Cluster_37_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_37, file = "TMS_merge_T_cells_Cluster_37_VariableGenes.rds")
Cluster_37 <- readRDS(file = "TMS_merge_T_cells_Cluster_37_VariableGenes.rds")

all.genes <- rownames(Cluster_37)
Cluster_37 <- ScaleData(Cluster_37, features = all.genes)

saveRDS(Cluster_37, file = "TMS_merge_T_cells_Cluster_37_ScaledData.rds")
Cluster_37 <- readRDS(file = "TMS_merge_T_cells_Cluster_37_ScaledData.rds")

Cluster_37 <- RunPCA(Cluster_37, npcs = 2, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_37, reduction = "pca")

saveRDS(Cluster_37, file = "TMS_merge_T_cells_Cluster_37_PCA.rds")
Cluster_37 <- readRDS(file = "TMS_merge_T_cells_Cluster_37_PCA.rds")

Cluster_37 <- JackStraw(Cluster_37, num.replicate = 100, dims = 2)
Cluster_37 <- ScoreJackStraw(Cluster_37, dims = 1:2)
JackStrawPlot(Cluster_37, dims = 1:2)

ElbowPlot(Cluster_37, ndims = 2)

Cluster_37 <- RunTSNE(Cluster_37, dims = 1:2, perplexity = 10)
DimPlot(Cluster_37, reduction = "tsne")

Cluster_37 <- RunUMAP(Cluster_37, dims = 1:2)
DimPlot(Cluster_37, reduction = "umap")

saveRDS(Cluster_37, file = "TMS_merge_T_cells_Cluster_37_TSNE_UMAP.rds")
Cluster_37 <- readRDS(file = "TMS_merge_T_cells_Cluster_37_TSNE_UMAP.rds")

Cluster_37 <- FindNeighbors(Cluster_37, dims = 1:2)

#louvian
Cluster_37 <- FindClusters(Cluster_37, resolution = seq(0,3, by = 0.1))

clustree(Cluster_37, prefix = "RNA_snn_res.")
clustree(Cluster_37, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_37) <- Cluster_37$RNA_snn_res.3

Cluster_37 <- BuildClusterTree(Cluster_37, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
PlotClusterTree(object = Cluster_37)

Cluster_37.markers <- FindAllMarkers(Cluster_37)
m <- Cluster_37.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_37, features = m$gene)
write.csv(Cluster_37.markers, "TMS_T_cells_Cluster_37_merkers.csv")

Cluster_37 <- RenameIdents(Cluster_37, "1" = "CD8 T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD8 T", 
                           "6" = "CD8 T", "7" = "CD8 T", "8" = "CD8 T", "9" = "CD8 T", "10" = "CD8 T", "11" = "CD4 T",
                           "12" = "CD4 T")

saveRDS(Cluster_37, file = "TMS_merge_T_cells_Cluster_37_clustered.rds")
Cluster_37 <- readRDS(file = "TMS_merge_T_cells_Cluster_37_clustered.rds")

