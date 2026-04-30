Cluster_32 <- readRDS(file = "TMS_merge_T_cells_clusters_32.rds")

#Highly variable expressed genes
Cluster_32 <- FindVariableFeatures(Cluster_32, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_32),15)
plot1 <- VariableFeaturePlot(Cluster_32)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_32), file = "TMS_merge_T_cells_Cluster_32_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_32, file = "TMS_merge_T_cells_Cluster_32_VariableGenes.rds")
Cluster_32 <- readRDS(file = "TMS_merge_T_cells_Cluster_32_VariableGenes.rds")

all.genes <- rownames(Cluster_32)
Cluster_32 <- ScaleData(Cluster_32, features = all.genes)

saveRDS(Cluster_32, file = "TMS_merge_T_cells_Cluster_32_ScaledData.rds")
Cluster_32 <- readRDS(file = "TMS_merge_T_cells_Cluster_32_ScaledData.rds")

Cluster_32 <- RunPCA(Cluster_32, npcs = 11, ndims.print = 1:11, nfeatures.print = 5)
DimPlot(Cluster_32, reduction = "pca")

saveRDS(Cluster_32, file = "TMS_merge_T_cells_Cluster_32_PCA.rds")
Cluster_32 <- readRDS(file = "TMS_merge_T_cells_Cluster_32_PCA.rds")

Cluster_32 <- JackStraw(Cluster_32, num.replicate = 100, dims = 11)
Cluster_32 <- ScoreJackStraw(Cluster_32, dims = 1:11)
JackStrawPlot(Cluster_32, dims = 1:11)

ElbowPlot(Cluster_32, ndims = 11)

Cluster_32 <- RunTSNE(Cluster_32, dims = 1:11, perplexity = 10)
DimPlot(Cluster_32, reduction = "tsne")

Cluster_32 <- RunUMAP(Cluster_32, dims = 1:11)
DimPlot(Cluster_32, reduction = "umap")

saveRDS(Cluster_32, file = "TMS_merge_T_cells_Cluster_32_TSNE_UMAP.rds")
Cluster_32 <- readRDS(file = "TMS_merge_T_cells_Cluster_32_TSNE_UMAP.rds")

Cluster_32 <- FindNeighbors(Cluster_32, dims = 1:11)

#louvian
Cluster_32 <- FindClusters(Cluster_32, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_32, prefix = "RNA_snn_res.")
clustree(Cluster_32, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_32) <- Cluster_32$RNA_snn_res.0

#Cluster_32 <- BuildClusterTree(Cluster_32, reorder.numeric = TRUE, reorder = TRUE, dims = 1:11)
#PlotClusterTree(object = Cluster_32)

#Cluster_32.markers <- FindAllMarkers(Cluster_32)
#m <- Cluster_32.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_32, features = m$gene)
#write.csv(Cluster_32.markers, "TMS_T_cells_Cluster_32_merkers.csv")

Cluster_32 <- RenameIdents(Cluster_32, "0" = "DP T")

saveRDS(Cluster_32, file = "TMS_merge_T_cells_Cluster_32_clustered.rds")
Cluster_32 <- readRDS(file = "TMS_merge_T_cells_Cluster_32_clustered.rds")

