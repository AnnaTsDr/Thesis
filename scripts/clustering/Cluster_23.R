Cluster_23 <- readRDS(file = "TMS_merge_T_cells_clusters_23.rds")

#Highly variable expressed genes
Cluster_23 <- FindVariableFeatures(Cluster_23, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.3, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_23),15)
plot1 <- VariableFeaturePlot(Cluster_23)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_23), file = "TMS_merge_T_cells_Cluster_23_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_23, file = "TMS_merge_T_cells_Cluster_23_VariableGenes.rds")
Cluster_23 <- readRDS(file = "TMS_merge_T_cells_Cluster_23_VariableGenes.rds")

all.genes <- rownames(Cluster_23)
Cluster_23 <- ScaleData(Cluster_23, features = all.genes)

saveRDS(Cluster_23, file = "TMS_merge_T_cells_Cluster_23_ScaledData.rds")
Cluster_23 <- readRDS(file = "TMS_merge_T_cells_Cluster_23_ScaledData.rds")

Cluster_23 <- RunPCA(Cluster_23, npcs = 3, ndims.print = 1:3, nfeatures.print = 5)
DimPlot(Cluster_23, reduction = "pca")

saveRDS(Cluster_23, file = "TMS_merge_T_cells_Cluster_23_PCA.rds")
Cluster_23 <- readRDS(file = "TMS_merge_T_cells_Cluster_23_PCA.rds")

Cluster_23 <- JackStraw(Cluster_23, num.replicate = 100, dims = 3)
Cluster_23 <- ScoreJackStraw(Cluster_23, dims = 1:3)
JackStrawPlot(Cluster_23, dims = 1:3)

ElbowPlot(Cluster_23, ndims = 3)

Cluster_23 <- RunTSNE(Cluster_23, dims = 1:3, perplexity = 10)
DimPlot(Cluster_23, reduction = "tsne")

Cluster_23 <- RunUMAP(Cluster_23, dims = 1:3)
DimPlot(Cluster_23, reduction = "umap")

saveRDS(Cluster_23, file = "TMS_merge_T_cells_Cluster_23_TSNE_UMAP.rds")
Cluster_23 <- readRDS(file = "TMS_merge_T_cells_Cluster_23_TSNE_UMAP.rds")

Cluster_23 <- FindNeighbors(Cluster_23, dims = 1:3)

#louvian
Cluster_23 <- FindClusters(Cluster_23, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_23, prefix = "RNA_snn_res.")
clustree(Cluster_23, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_23) <- Cluster_23$RNA_snn_res.0

#Cluster_23 <- BuildClusterTree(Cluster_23, reorder.numeric = TRUE, reorder = TRUE, dims = 1:3)
#PlotClusterTree(object = Cluster_23)

#Cluster_23.markers <- FindAllMarkers(Cluster_23)
#m <- Cluster_23.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_23, features = m$gene)
#write.csv(Cluster_23.markers, "TMS_T_cells_Cluster_23_merkers.csv")

Cluster_23 <- RenameIdents(Cluster_23, "0" = "CD8 T")

saveRDS(Cluster_23, file = "TMS_merge_T_cells_Cluster_23_clustered.rds")
Cluster_23 <- readRDS(file = "TMS_merge_T_cells_Cluster_23_clustered.rds")

