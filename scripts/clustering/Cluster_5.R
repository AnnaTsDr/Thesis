Cluster_5 <- readRDS(file = "TMS_merge_T_cells_clusters_5.rds")

#Highly variable expressed genes
Cluster_5 <- FindVariableFeatures(Cluster_5, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_5),15)
plot1 <- VariableFeaturePlot(Cluster_5)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_5), file = "TMS_merge_T_cells_Cluster_5_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_5, file = "TMS_merge_T_cells_Cluster_5_VariableGenes.rds")
Cluster_5 <- readRDS(file = "TMS_merge_T_cells_Cluster_5_VariableGenes.rds")

all.genes <- rownames(Cluster_5)
Cluster_5 <- ScaleData(Cluster_5, features = all.genes)

saveRDS(Cluster_5, file = "TMS_merge_T_cells_Cluster_5_ScaledData.rds")
Cluster_5 <- readRDS(file = "TMS_merge_T_cells_Cluster_5_ScaledData.rds")

Cluster_5 <- RunPCA(Cluster_5, npcs = 20, ndims.print = 1:20, nfeatures.print = 5)
DimPlot(Cluster_5, reduction = "pca")

saveRDS(Cluster_5, file = "TMS_merge_T_cells_Cluster_5_PCA.rds")
Cluster_5 <- readRDS(file = "TMS_merge_T_cells_Cluster_5_PCA.rds")

Cluster_5 <- JackStraw(Cluster_5, num.replicate = 100, dims = 20)
Cluster_5 <- ScoreJackStraw(Cluster_5, dims = 1:20)
JackStrawPlot(Cluster_5, dims = 1:20)

ElbowPlot(Cluster_5, ndims = 20)

Cluster_5 <- RunTSNE(Cluster_5, dims = 1:20, perplexity = 10)
DimPlot(Cluster_5, reduction = "tsne")

Cluster_5 <- RunUMAP(Cluster_5, dims = 1:20)
DimPlot(Cluster_5, reduction = "umap")

saveRDS(Cluster_5, file = "TMS_merge_T_cells_Cluster_5_TSNE_UMAP.rds")
Cluster_5 <- readRDS(file = "TMS_merge_T_cells_Cluster_5_TSNE_UMAP.rds")

Cluster_5 <- FindNeighbors(Cluster_5, dims = 1:20)

#louvian
Cluster_5 <- FindClusters(Cluster_5, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_5, prefix = "RNA_snn_res.")
clustree(Cluster_5, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_5) <- Cluster_5$RNA_snn_res.0

#Cluster_5 <- BuildClusterTree(Cluster_5, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
#PlotClusterTree(object = Cluster_5)

#Cluster_5.markers <- FindAllMarkers(Cluster_5)
#m <- Cluster_5.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_5, features = m$gene)
#write.csv(Cluster_5.markers, "TMS_T_cells_Cluster_5_merkers.csv")

Cluster_5 <- RenameIdents(Cluster_5, "0" = "DN T")

saveRDS(Cluster_5, file = "TMS_merge_T_cells_Cluster_5_clustered.rds")
Cluster_5 <- readRDS(file = "TMS_merge_T_cells_Cluster_5_clustered.rds")
