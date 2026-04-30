Cluster_30 <- readRDS(file = "TMS_merge_T_cells_clusters_30.rds")

#Highly variable expressed genes
Cluster_30 <- FindVariableFeatures(Cluster_30, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.8, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_30),15)
plot1 <- VariableFeaturePlot(Cluster_30)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_30), file = "TMS_merge_T_cells_Cluster_30_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_30, file = "TMS_merge_T_cells_Cluster_30_VariableGenes.rds")
Cluster_30 <- readRDS(file = "TMS_merge_T_cells_Cluster_30_VariableGenes.rds")

all.genes <- rownames(Cluster_30)
Cluster_30 <- ScaleData(Cluster_30, features = all.genes)

saveRDS(Cluster_30, file = "TMS_merge_T_cells_Cluster_30_ScaledData.rds")
Cluster_30 <- readRDS(file = "TMS_merge_T_cells_Cluster_30_ScaledData.rds")

Cluster_30 <- RunPCA(Cluster_30, npcs = 22, ndims.print = 1:22, nfeatures.print = 5)
DimPlot(Cluster_30, reduction = "pca")

saveRDS(Cluster_30, file = "TMS_merge_T_cells_Cluster_30_PCA.rds")
Cluster_30 <- readRDS(file = "TMS_merge_T_cells_Cluster_30_PCA.rds")

Cluster_30 <- JackStraw(Cluster_30, num.replicate = 100, dims = 22)
Cluster_30 <- ScoreJackStraw(Cluster_30, dims = 1:22)
JackStrawPlot(Cluster_30, dims = 1:22)

ElbowPlot(Cluster_30, ndims = 22)

Cluster_30 <- RunTSNE(Cluster_30, dims = 1:22, perplexity = 10)
DimPlot(Cluster_30, reduction = "tsne")

Cluster_30 <- RunUMAP(Cluster_30, dims = 1:22)
DimPlot(Cluster_30, reduction = "umap")

saveRDS(Cluster_30, file = "TMS_merge_T_cells_Cluster_30_TSNE_UMAP.rds")
Cluster_30 <- readRDS(file = "TMS_merge_T_cells_Cluster_30_TSNE_UMAP.rds")

Cluster_30 <- FindNeighbors(Cluster_30, dims = 1:22)

#louvian
Cluster_30 <- FindClusters(Cluster_30, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_30, prefix = "RNA_snn_res.")
clustree(Cluster_30, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_30) <- Cluster_30$RNA_snn_res.0

#Cluster_30 <- BuildClusterTree(Cluster_30, reorder.numeric = TRUE, reorder = TRUE, dims = 1:22)
#PlotClusterTree(object = Cluster_30)

#Cluster_30.markers <- FindAllMarkers(Cluster_30)
#m <- Cluster_30.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_30, features = m$gene)
#write.csv(Cluster_30.markers, "TMS_T_cells_Cluster_30_merkers.csv")

Cluster_30 <- RenameIdents(Cluster_30, "0" = "CD4 T")


saveRDS(Cluster_30, file = "TMS_merge_T_cells_Cluster_30_clustered.rds")
Cluster_30 <- readRDS(file = "TMS_merge_T_cells_Cluster_30_clustered.rds")

