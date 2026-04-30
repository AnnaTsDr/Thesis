Cluster_9 <- readRDS(file = "TMS_merge_T_cells_clusters_9.rds")

#Highly variable expressed genes
Cluster_9 <- FindVariableFeatures(Cluster_9, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.8, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_9),15)
plot1 <- VariableFeaturePlot(Cluster_9)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_9), file = "TMS_merge_T_cells_Cluster_9_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_9, file = "TMS_merge_T_cells_Cluster_9_VariableGenes.rds")
Cluster_9 <- readRDS(file = "TMS_merge_T_cells_Cluster_9_VariableGenes.rds")

all.genes <- rownames(Cluster_9)
Cluster_9 <- ScaleData(Cluster_9, features = all.genes)

saveRDS(Cluster_9, file = "TMS_merge_T_cells_Cluster_9_ScaledData.rds")
Cluster_9 <- readRDS(file = "TMS_merge_T_cells_Cluster_9_ScaledData.rds")

Cluster_9 <- RunPCA(Cluster_9, npcs = 13, ndims.print = 1:13, nfeatures.print = 5)
DimPlot(Cluster_9, reduction = "pca")

saveRDS(Cluster_9, file = "TMS_merge_T_cells_Cluster_9_PCA.rds")
Cluster_9 <- readRDS(file = "TMS_merge_T_cells_Cluster_9_PCA.rds")

Cluster_9 <- JackStraw(Cluster_9, num.replicate = 100, dims = 13)
Cluster_9 <- ScoreJackStraw(Cluster_9, dims = 1:13)
JackStrawPlot(Cluster_9, dims = 1:13)

ElbowPlot(Cluster_9, ndims = 13)

Cluster_9 <- RunTSNE(Cluster_9, dims = 1:13, perplexity = 10)
DimPlot(Cluster_9, reduction = "tsne")

Cluster_9 <- RunUMAP(Cluster_9, dims = 1:13)
DimPlot(Cluster_9, reduction = "umap")

saveRDS(Cluster_9, file = "TMS_merge_T_cells_Cluster_9_TSNE_UMAP.rds")
Cluster_9 <- readRDS(file = "TMS_merge_T_cells_Cluster_9_TSNE_UMAP.rds")

Cluster_9 <- FindNeighbors(Cluster_9, dims = 1:13)

#louvian
Cluster_9 <- FindClusters(Cluster_9, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_9, prefix = "RNA_snn_res.")
clustree(Cluster_9, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_9) <- Cluster_9$RNA_snn_res.0

#Cluster_9 <- BuildClusterTree(Cluster_9, reorder.numeric = TRUE, reorder = TRUE, dims = 1:13)
#PlotClusterTree(object = Cluster_9)

#Cluster_9.markers <- FindAllMarkers(Cluster_9)
#m <- Cluster_9.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_9, features = m$gene)
#write.csv(Cluster_9.markers, "TMS_T_cells_Cluster_9_merkers.csv")

Cluster_9 <- RenameIdents(Cluster_9, "0" = "DN T")

saveRDS(Cluster_9, file = "TMS_merge_T_cells_Cluster_9_clustered.rds")
Cluster_9 <- readRDS(file = "TMS_merge_T_cells_Cluster_9_clustered.rds")

