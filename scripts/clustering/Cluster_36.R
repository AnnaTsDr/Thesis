Cluster_36 <- readRDS(file = "TMS_merge_T_cells_clusters_36.rds")

#Highly variable expressed genes
Cluster_36 <- FindVariableFeatures(Cluster_36, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.3, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_36),15)
plot1 <- VariableFeaturePlot(Cluster_36)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_36), file = "TMS_merge_T_cells_Cluster_36_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_36, file = "TMS_merge_T_cells_Cluster_36_VariableGenes.rds")
Cluster_36 <- readRDS(file = "TMS_merge_T_cells_Cluster_36_VariableGenes.rds")

all.genes <- rownames(Cluster_36)
Cluster_36 <- ScaleData(Cluster_36, features = all.genes)

saveRDS(Cluster_36, file = "TMS_merge_T_cells_Cluster_36_ScaledData.rds")
Cluster_36 <- readRDS(file = "TMS_merge_T_cells_Cluster_36_ScaledData.rds")

Cluster_36 <- RunPCA(Cluster_36, npcs = 9, ndims.print = 1:9, nfeatures.print = 5)
DimPlot(Cluster_36, reduction = "pca")

saveRDS(Cluster_36, file = "TMS_merge_T_cells_Cluster_36_PCA.rds")
Cluster_36 <- readRDS(file = "TMS_merge_T_cells_Cluster_36_PCA.rds")

Cluster_36 <- JackStraw(Cluster_36, num.replicate = 100, dims = 9)
Cluster_36 <- ScoreJackStraw(Cluster_36, dims = 1:9)
JackStrawPlot(Cluster_36, dims = 1:9)

ElbowPlot(Cluster_36, ndims = 9)

Cluster_36 <- RunTSNE(Cluster_36, dims = 1:9, perplexity = 10)
DimPlot(Cluster_36, reduction = "tsne")

Cluster_36 <- RunUMAP(Cluster_36, dims = 1:9)
DimPlot(Cluster_36, reduction = "umap")

saveRDS(Cluster_36, file = "TMS_merge_T_cells_Cluster_36_TSNE_UMAP.rds")
Cluster_36 <- readRDS(file = "TMS_merge_T_cells_Cluster_36_TSNE_UMAP.rds")

Cluster_36 <- FindNeighbors(Cluster_36, dims = 1:9)

#louvian
Cluster_36 <- FindClusters(Cluster_36, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_36, prefix = "RNA_snn_res.")
clustree(Cluster_36, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_36) <- Cluster_36$RNA_snn_res.0

#Cluster_36 <- BuildClusterTree(Cluster_36, reorder.numeric = TRUE, reorder = TRUE, dims = 1:9)
#PlotClusterTree(object = Cluster_36)

#Cluster_36.markers <- FindAllMarkers(Cluster_36)
#m <- Cluster_36.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_36, features = m$gene)
#write.csv(Cluster_36.markers, "TMS_T_cells_Cluster_36_merkers.csv")

Cluster_36 <- RenameIdents(Cluster_36, "0" = "CD8 T")

saveRDS(Cluster_36, file = "TMS_merge_T_cells_Cluster_36_clustered.rds")
Cluster_36 <- readRDS(file = "TMS_merge_T_cells_Cluster_36_clustered.rds")

