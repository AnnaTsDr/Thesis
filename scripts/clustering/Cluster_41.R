Cluster_41 <- readRDS(file = "TMS_merge_T_cells_clusters_41.rds")

#Highly variable expressed genes
Cluster_41 <- FindVariableFeatures(Cluster_41, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_41),15)
plot1 <- VariableFeaturePlot(Cluster_41)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_41), file = "TMS_merge_T_cells_Cluster_41_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_41, file = "TMS_merge_T_cells_Cluster_41_VariableGenes.rds")
Cluster_41 <- readRDS(file = "TMS_merge_T_cells_Cluster_41_VariableGenes.rds")

all.genes <- rownames(Cluster_41)
Cluster_41 <- ScaleData(Cluster_41, features = all.genes)

saveRDS(Cluster_41, file = "TMS_merge_T_cells_Cluster_41_ScaledData.rds")
Cluster_41 <- readRDS(file = "TMS_merge_T_cells_Cluster_41_ScaledData.rds")

Cluster_41 <- RunPCA(Cluster_41, npcs = 27, ndims.print = 1:27, nfeatures.print = 5)
DimPlot(Cluster_41, reduction = "pca")

saveRDS(Cluster_41, file = "TMS_merge_T_cells_Cluster_41_PCA.rds")
Cluster_41 <- readRDS(file = "TMS_merge_T_cells_Cluster_41_PCA.rds")

Cluster_41 <- JackStraw(Cluster_41, num.replicate = 100, dims = 27)
Cluster_41 <- ScoreJackStraw(Cluster_41, dims = 1:27)
JackStrawPlot(Cluster_41, dims = 1:27)

ElbowPlot(Cluster_41, ndims = 27)

Cluster_41 <- RunTSNE(Cluster_41, dims = 1:27, perplexity = 10)
DimPlot(Cluster_41, reduction = "tsne")

Cluster_41 <- RunUMAP(Cluster_41, dims = 1:27)
DimPlot(Cluster_41, reduction = "umap")

saveRDS(Cluster_41, file = "TMS_merge_T_cells_Cluster_41_TSNE_UMAP.rds")
Cluster_41 <- readRDS(file = "TMS_merge_T_cells_Cluster_41_TSNE_UMAP.rds")

Cluster_41 <- FindNeighbors(Cluster_41, dims = 1:27)

#louvian
Cluster_41 <- FindClusters(Cluster_41, resolution = seq(0,3, by = 0.1))

clustree(Cluster_41, prefix = "RNA_snn_res.")
clustree(Cluster_41, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_41) <- Cluster_41$RNA_snn_res.0

#Cluster_41 <- BuildClusterTree(Cluster_41, reorder.numeric = TRUE, reorder = TRUE, dims = 1:27)
#PlotClusterTree(object = Cluster_41)

#Cluster_41.markers <- FindAllMarkers(Cluster_41)
#m <- Cluster_41.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_41, features = m$gene)
#write.csv(Cluster_41.markers, "TMS_T_cells_Cluster_41_merkers.csv")

Cluster_41 <- RenameIdents(Cluster_41, "0" = "CD8 T")

saveRDS(Cluster_41, file = "TMS_merge_T_cells_Cluster_41_clustered.rds")
Cluster_41 <- readRDS(file = "TMS_merge_T_cells_Cluster_41_clustered.rds")


