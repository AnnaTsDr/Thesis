Cluster_10 <- readRDS(file = "TMS_merge_T_cells_clusters_10.rds")

#Highly variable expressed genes
Cluster_10 <- FindVariableFeatures(Cluster_10, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_10),15)
plot1 <- VariableFeaturePlot(Cluster_10)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_10), file = "TMS_merge_T_cells_Cluster_10_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_10, file = "TMS_merge_T_cells_Cluster_10_VariableGenes.rds")
Cluster_10 <- readRDS(file = "TMS_merge_T_cells_Cluster_10_VariableGenes.rds")

all.genes <- rownames(Cluster_10)
Cluster_10 <- ScaleData(Cluster_10, features = all.genes)

saveRDS(Cluster_10, file = "TMS_merge_T_cells_Cluster_10_ScaledData.rds")
Cluster_10 <- readRDS(file = "TMS_merge_T_cells_Cluster_10_ScaledData.rds")

Cluster_10 <- RunPCA(Cluster_10, npcs = 5, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(Cluster_10, reduction = "pca")

saveRDS(Cluster_10, file = "TMS_merge_T_cells_Cluster_10_PCA.rds")
Cluster_10 <- readRDS(file = "TMS_merge_T_cells_Cluster_10_PCA.rds")

Cluster_10 <- JackStraw(Cluster_10, num.replicate = 100, dims = 5)
Cluster_10 <- ScoreJackStraw(Cluster_10, dims = 1:5)
JackStrawPlot(Cluster_10, dims = 1:5)

ElbowPlot(Cluster_10, ndims = 5)

Cluster_10 <- RunTSNE(Cluster_10, dims = 1:5, perplexity = 10)
DimPlot(Cluster_10, reduction = "tsne")

Cluster_10 <- RunUMAP(Cluster_10, dims = 1:5)
DimPlot(Cluster_10, reduction = "umap")

saveRDS(Cluster_10, file = "TMS_merge_T_cells_Cluster_10_TSNE_UMAP.rds")
Cluster_10 <- readRDS(file = "TMS_merge_T_cells_Cluster_10_TSNE_UMAP.rds")

Cluster_10 <- FindNeighbors(Cluster_10, dims = 1:5)

#louvian
Cluster_10 <- FindClusters(Cluster_10, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_10, prefix = "RNA_snn_res.")
clustree(Cluster_10, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_10) <- Cluster_10$RNA_snn_res.0

#Cluster_10 <- BuildClusterTree(Cluster_10, reorder.numeric = TRUE, reorder = TRUE, dims = 1:5)
#PlotClusterTree(object = Cluster_10)

#Cluster_10.markers <- FindAllMarkers(Cluster_10)
#m <- Cluster_10.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_10, features = m$gene)
#write.csv(Cluster_10.markers, "TMS_T_cells_Cluster_10_merkers.csv")

Cluster_10 <- RenameIdents(Cluster_10, "0" = "DN T")

saveRDS(Cluster_10, file = "TMS_merge_T_cells_Cluster_10_clustered.rds")
Cluster_10 <- readRDS(file = "TMS_merge_T_cells_Cluster_10_clustered.rds")

