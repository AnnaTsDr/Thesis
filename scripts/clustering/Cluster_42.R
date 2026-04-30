Cluster_42 <- readRDS(file = "TMS_merge_T_cells_clusters_42.rds")

#Highly variable expressed genes
Cluster_42 <- FindVariableFeatures(Cluster_42, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_42),15)
plot1 <- VariableFeaturePlot(Cluster_42)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_42), file = "TMS_merge_T_cells_Cluster_42_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_42, file = "TMS_merge_T_cells_Cluster_42_VariableGenes.rds")
Cluster_42 <- readRDS(file = "TMS_merge_T_cells_Cluster_42_VariableGenes.rds")

all.genes <- rownames(Cluster_42)
Cluster_42 <- ScaleData(Cluster_42, features = all.genes)

saveRDS(Cluster_42, file = "TMS_merge_T_cells_Cluster_42_ScaledData.rds")
Cluster_42 <- readRDS(file = "TMS_merge_T_cells_Cluster_42_ScaledData.rds")

Cluster_42 <- RunPCA(Cluster_42, npcs = 2, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_42, reduction = "pca")

saveRDS(Cluster_42, file = "TMS_merge_T_cells_Cluster_42_PCA.rds")
Cluster_42 <- readRDS(file = "TMS_merge_T_cells_Cluster_42_PCA.rds")

Cluster_42 <- JackStraw(Cluster_42, num.replicate = 100, dims = 2)
Cluster_42 <- ScoreJackStraw(Cluster_42, dims = 1:2)
JackStrawPlot(Cluster_42, dims = 1:2)

ElbowPlot(Cluster_42, ndims = 2)

Cluster_42 <- RunTSNE(Cluster_42, dims = 1:2, perplexity = 10)
DimPlot(Cluster_42, reduction = "tsne")

Cluster_42 <- RunUMAP(Cluster_42, dims = 1:2)
DimPlot(Cluster_42, reduction = "umap")

saveRDS(Cluster_42, file = "TMS_merge_T_cells_Cluster_42_TSNE_UMAP.rds")
Cluster_42 <- readRDS(file = "TMS_merge_T_cells_Cluster_42_TSNE_UMAP.rds")

Cluster_42 <- FindNeighbors(Cluster_42, dims = 1:2)

#louvian
Cluster_42 <- FindClusters(Cluster_42, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_42, prefix = "RNA_snn_res.")
clustree(Cluster_42, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_42) <- Cluster_42$RNA_snn_res.0

#Cluster_42 <- BuildClusterTree(Cluster_42, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
#PlotClusterTree(object = Cluster_42)

#Cluster_42.markers <- FindAllMarkers(Cluster_42)
#m <- Cluster_42.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_42, features = m$gene)
#write.csv(Cluster_42.markers, "TMS_T_cells_Cluster_42_merkers.csv")

Cluster_42 <- RenameIdents(Cluster_42, "0" = "CD8 T")

saveRDS(Cluster_42, file = "TMS_merge_T_cells_Cluster_42_clustered.rds")
Cluster_42 <- readRDS(file = "TMS_merge_T_cells_Cluster_42_clustered.rds")
