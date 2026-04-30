Cluster_47 <- readRDS(file = "TMS_merge_T_cells_clusters_47.rds")

#Highly variable expressed genes
Cluster_47 <- FindVariableFeatures(Cluster_47, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_47),15)
plot1 <- VariableFeaturePlot(Cluster_47)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_47), file = "TMS_merge_T_cells_Cluster_47_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_47, file = "TMS_merge_T_cells_Cluster_47_VariableGenes.rds")
Cluster_47 <- readRDS(file = "TMS_merge_T_cells_Cluster_47_VariableGenes.rds")

all.genes <- rownames(Cluster_47)
Cluster_47 <- ScaleData(Cluster_47, features = all.genes)

saveRDS(Cluster_47, file = "TMS_merge_T_cells_Cluster_47_ScaledData.rds")
Cluster_47 <- readRDS(file = "TMS_merge_T_cells_Cluster_47_ScaledData.rds")

Cluster_47 <- RunPCA(Cluster_47, npcs = 3, ndims.print = 1:3, nfeatures.print = 5)
DimPlot(Cluster_47, reduction = "pca")

saveRDS(Cluster_47, file = "TMS_merge_T_cells_Cluster_47_PCA.rds")
Cluster_47 <- readRDS(file = "TMS_merge_T_cells_Cluster_47_PCA.rds")

Cluster_47 <- JackStraw(Cluster_47, num.replicate = 100, dims = 3)
Cluster_47 <- ScoreJackStraw(Cluster_47, dims = 1:3)
JackStrawPlot(Cluster_47, dims = 1:3)

ElbowPlot(Cluster_47, ndims = 3)

Cluster_47 <- RunTSNE(Cluster_47, dims = 1:3, perplexity = 10)
DimPlot(Cluster_47, reduction = "tsne")

Cluster_47 <- RunUMAP(Cluster_47, dims = 1:3)
DimPlot(Cluster_47, reduction = "umap")

saveRDS(Cluster_47, file = "TMS_merge_T_cells_Cluster_47_TSNE_UMAP.rds")
Cluster_47 <- readRDS(file = "TMS_merge_T_cells_Cluster_47_TSNE_UMAP.rds")

Cluster_47 <- FindNeighbors(Cluster_47, dims = 1:3)

#louvian
Cluster_47 <- FindClusters(Cluster_47, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_47, prefix = "RNA_snn_res.")
clustree(Cluster_47, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_47) <- Cluster_47$RNA_snn_res.0

#Cluster_47 <- BuildClusterTree(Cluster_47, reorder.numeric = TRUE, reorder = TRUE, dims = 1:3)
#PlotClusterTree(object = Cluster_47)

#Cluster_47.markers <- FindAllMarkers(Cluster_47)
#m <- Cluster_47.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_47, features = m$gene)
#write.csv(Cluster_47.markers, "TMS_T_cells_Cluster_47_merkers.csv")

Cluster_47 <- RenameIdents(Cluster_47, "0" = "DP T")

saveRDS(Cluster_47, file = "TMS_merge_T_cells_Cluster_47_clustered.rds")
Cluster_47 <- readRDS(file = "TMS_merge_T_cells_Cluster_47_clustered.rds")

