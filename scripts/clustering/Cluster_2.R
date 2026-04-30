Cluster_2 <- readRDS(file = "TMS_merge_T_cells_Clusters_2.rds")

#Highly variable expressed genes
Cluster_2 <- FindVariableFeatures(Cluster_2, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.5, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_2),15)
plot1 <- VariableFeaturePlot(Cluster_2)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_2), file = "TMS_merge_T_cells_Cluster_2_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_2, file = "TMS_merge_T_cells_Cluster_2_VariableGenes.rds")
Cluster_2 <- readRDS(file = "TMS_merge_T_cells_Cluster_2_VariableGenes.rds")

all.genes <- rownames(Cluster_2)
Cluster_2 <- ScaleData(Cluster_2, features = all.genes)

saveRDS(Cluster_2, file = "TMS_merge_T_cells_Cluster_2_ScaledData.rds")
Cluster_2 <- readRDS(file = "TMS_merge_T_cells_Cluster_2_ScaledData.rds")

Cluster_2 <- RunPCA(Cluster_2, npcs = 2, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_2, reduction = "pca")

saveRDS(Cluster_2, file = "TMS_merge_T_cells_Cluster_2_PCA.rds")
Cluster_2 <- readRDS(file = "TMS_merge_T_cells_Cluster_2_PCA.rds")

Cluster_2 <- JackStraw(Cluster_2, num.replicate = 100, dims = 2)
Cluster_2 <- ScoreJackStraw(Cluster_2, dims = 1:2)
JackStrawPlot(Cluster_2, dims = 1:2)

ElbowPlot(Cluster_2, ndims = 2)

Cluster_2 <- RunTSNE(Cluster_2, dims = 1:2, perplexity = 10)
DimPlot(Cluster_2, reduction = "tsne")

Cluster_2 <- RunUMAP(Cluster_2, dims = 1:2)
DimPlot(Cluster_2, reduction = "umap")

saveRDS(Cluster_2, file = "TMS_merge_T_cells_Cluster_2_TSNE_UMAP.rds")
Cluster_2 <- readRDS(file = "TMS_merge_T_cells_Cluster_2_TSNE_UMAP.rds")

Cluster_2 <- FindNeighbors(Cluster_2, dims = 1:2)

#louvian
Cluster_2 <- FindClusters(Cluster_2, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_2, prefix = "RNA_snn_res.")
clustree(Cluster_2, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_2) <- Cluster_2$RNA_snn_res.0.9

Cluster_2 <- BuildClusterTree(Cluster_2, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
PlotClusterTree(object = Cluster_2)

Cluster_2.markers <- FindAllMarkers(Cluster_2)
m <- Cluster_2.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_2, features = m$gene)
write.csv(Cluster_2.markers, "TMS_T_cells_Cluster_2_merkers.csv")

Cluster_2 <- RenameIdents(Cluster_2, "1" = "CD8 T", "2" = "DN T", "3" = "DN T")

saveRDS(Cluster_2, file = "TMS_merge_T_cells_Cluster_2_clustered.rds")
Cluster_2 <- readRDS(file = "TMS_merge_T_cells_Cluster_2_clustered.rds")

#homogenize cluster