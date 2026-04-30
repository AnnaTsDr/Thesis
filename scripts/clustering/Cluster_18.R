Cluster_18 <- readRDS(file = "TMS_merge_T_cells_clusters_18.rds")

#Highly variable expressed genes
Cluster_18 <- FindVariableFeatures(Cluster_18, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.6, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_18),15)
plot1 <- VariableFeaturePlot(Cluster_18)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_18), file = "TMS_merge_T_cells_Cluster_18_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_18, file = "TMS_merge_T_cells_Cluster_18_VariableGenes.rds")
Cluster_18 <- readRDS(file = "TMS_merge_T_cells_Cluster_18_VariableGenes.rds")

all.genes <- rownames(Cluster_18)
Cluster_18 <- ScaleData(Cluster_18, features = all.genes)

saveRDS(Cluster_18, file = "TMS_merge_T_cells_Cluster_18_ScaledData.rds")
Cluster_18 <- readRDS(file = "TMS_merge_T_cells_Cluster_18_ScaledData.rds")

Cluster_18 <- RunPCA(Cluster_18, npcs = 28, ndims.print = 1:28, nfeatures.print = 5)
DimPlot(Cluster_18, reduction = "pca")

saveRDS(Cluster_18, file = "TMS_merge_T_cells_Cluster_18_PCA.rds")
Cluster_18 <- readRDS(file = "TMS_merge_T_cells_Cluster_18_PCA.rds")

Cluster_18 <- JackStraw(Cluster_18, num.replicate = 100, dims = 28)
Cluster_18 <- ScoreJackStraw(Cluster_18, dims = 1:28)
JackStrawPlot(Cluster_18, dims = 1:28)

ElbowPlot(Cluster_18, ndims = 28)

Cluster_18 <- RunTSNE(Cluster_18, dims = 1:28, perplexity = 10)
DimPlot(Cluster_18, reduction = "tsne")

Cluster_18 <- RunUMAP(Cluster_18, dims = 1:28)
DimPlot(Cluster_18, reduction = "umap")

saveRDS(Cluster_18, file = "TMS_merge_T_cells_Cluster_18_TSNE_UMAP.rds")
Cluster_18 <- readRDS(file = "TMS_merge_T_cells_Cluster_18_TSNE_UMAP.rds")

Cluster_18 <- FindNeighbors(Cluster_18, dims = 1:28)

#louvian
Cluster_18 <- FindClusters(Cluster_18, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_18, prefix = "RNA_snn_res.")
clustree(Cluster_18, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_18) <- Cluster_18$RNA_snn_res.1

Cluster_18 <- BuildClusterTree(Cluster_18, reorder.numeric = TRUE, reorder = TRUE, dims = 1:28)
PlotClusterTree(object = Cluster_18)

Cluster_18.markers <- FindAllMarkers(Cluster_18)
m <- Cluster_18.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_18, features = m$gene)
write.csv(Cluster_18.markers, "TMS_T_cells_Cluster_18_merkers.csv")

Cluster_18 <- RenameIdents(Cluster_18, "1" = "CD8 T", "2" = "DP T", "3" = "DP T", "4" = "DP T", "5" = "DP T", 
                           "6" = "DP T")

saveRDS(Cluster_18, file = "TMS_merge_T_cells_Cluster_18_clustered.rds")
Cluster_18 <- readRDS(file = "TMS_merge_T_cells_Cluster_18_clustered.rds")
