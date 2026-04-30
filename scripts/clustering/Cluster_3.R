Cluster_3 <- readRDS(file = "TMS_merge_T_cells_Clusters_3.rds")

#Highly variable expressed genes
Cluster_3 <- FindVariableFeatures(Cluster_3, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_3),15)
plot1 <- VariableFeaturePlot(Cluster_3)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_3), file = "TMS_merge_T_cells_Cluster_3_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_3, file = "TMS_merge_T_cells_Cluster_3_VariableGenes.rds")
Cluster_3 <- readRDS(file = "TMS_merge_T_cells_Cluster_3_VariableGenes.rds")

all.genes <- rownames(Cluster_3)
Cluster_3 <- ScaleData(Cluster_3, features = all.genes)

saveRDS(Cluster_3, file = "TMS_merge_T_cells_Cluster_3_ScaledData.rds")
Cluster_3 <- readRDS(file = "TMS_merge_T_cells_Cluster_3_ScaledData.rds")

Cluster_3 <- RunPCA(Cluster_3, npcs = 2, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_3, reduction = "pca")

saveRDS(Cluster_3, file = "TMS_merge_T_cells_Cluster_3_PCA.rds")
Cluster_3 <- readRDS(file = "TMS_merge_T_cells_Cluster_3_PCA.rds")

Cluster_3 <- JackStraw(Cluster_3, num.replicate = 100, dims = 2)
Cluster_3 <- ScoreJackStraw(Cluster_3, dims = 1:2)
JackStrawPlot(Cluster_3, dims = 1:2)

ElbowPlot(Cluster_3, ndims = 2)

Cluster_3 <- RunTSNE(Cluster_3, dims = 1:2, perplexity = 10)
DimPlot(Cluster_3, reduction = "tsne")

Cluster_3 <- RunUMAP(Cluster_3, dims = 1:2)
DimPlot(Cluster_3, reduction = "umap")

saveRDS(Cluster_3, file = "TMS_merge_T_cells_Cluster_3_TSNE_UMAP.rds")
Cluster_3 <- readRDS(file = "TMS_merge_T_cells_Cluster_3_TSNE_UMAP.rds")

Cluster_3 <- FindNeighbors(Cluster_3, dims = 1:2)

#louvian
Cluster_3 <- FindClusters(Cluster_3, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_3, prefix = "RNA_snn_res.")
clustree(Cluster_3, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_3) <- Cluster_3$RNA_snn_res.1

Cluster_3 <- BuildClusterTree(Cluster_3, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
PlotClusterTree(object = Cluster_3)

Cluster_3.markers <- FindAllMarkers(Cluster_3)
m <- Cluster_3.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_3, features = m$gene)
write.csv(Cluster_3.markers, "TMS_T_cells_Cluster_3_merkers.csv")

Cluster_3 <- RenameIdents(Cluster_3, "1" = "Not T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T")

saveRDS(Cluster_3, file = "TMS_merge_T_cells_Cluster_3_clustered.rds")
Cluster_3 <- readRDS(file = "TMS_merge_T_cells_Cluster_3_clustered.rds")

#maybe not T