Cluster_35 <- readRDS(file = "TMS_merge_T_cells_clusters_35.rds")

#Highly variable expressed genes
Cluster_35 <- FindVariableFeatures(Cluster_35, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.3, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_35),15)
plot1 <- VariableFeaturePlot(Cluster_35)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_35), file = "TMS_merge_T_cells_Cluster_35_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_35, file = "TMS_merge_T_cells_Cluster_35_VariableGenes.rds")
Cluster_35 <- readRDS(file = "TMS_merge_T_cells_Cluster_35_VariableGenes.rds")

all.genes <- rownames(Cluster_35)
Cluster_35 <- ScaleData(Cluster_35, features = all.genes)

saveRDS(Cluster_35, file = "TMS_merge_T_cells_Cluster_35_ScaledData.rds")
Cluster_35 <- readRDS(file = "TMS_merge_T_cells_Cluster_35_ScaledData.rds")

Cluster_35 <- RunPCA(Cluster_35, npcs = 2, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_35, reduction = "pca")

saveRDS(Cluster_35, file = "TMS_merge_T_cells_Cluster_35_PCA.rds")
Cluster_35 <- readRDS(file = "TMS_merge_T_cells_Cluster_35_PCA.rds")

Cluster_35 <- JackStraw(Cluster_35, num.replicate = 100, dims = 2)
Cluster_35 <- ScoreJackStraw(Cluster_35, dims = 1:2)
JackStrawPlot(Cluster_35, dims = 1:2)

ElbowPlot(Cluster_35, ndims = 2)

Cluster_35 <- RunTSNE(Cluster_35, dims = 1:2, perplexity = 10)
DimPlot(Cluster_35, reduction = "tsne")

Cluster_35 <- RunUMAP(Cluster_35, dims = 1:2)
DimPlot(Cluster_35, reduction = "umap")

saveRDS(Cluster_35, file = "TMS_merge_T_cells_Cluster_35_TSNE_UMAP.rds")
Cluster_35 <- readRDS(file = "TMS_merge_T_cells_Cluster_35_TSNE_UMAP.rds")

Cluster_35 <- FindNeighbors(Cluster_35, dims = 1:2)

#louvian
Cluster_35 <- FindClusters(Cluster_35, resolution = seq(0,3, by = 0.1))

clustree(Cluster_35, prefix = "RNA_snn_res.")
clustree(Cluster_35, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_35) <- Cluster_35$RNA_snn_res.3

Cluster_35 <- BuildClusterTree(Cluster_35, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
PlotClusterTree(object = Cluster_35)

Cluster_35.markers <- FindAllMarkers(Cluster_35)
m <- Cluster_35.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_35, features = m$gene)
write.csv(Cluster_35.markers, "TMS_T_cells_Cluster_35_merkers.csv")

Cluster_35 <- RenameIdents(Cluster_35, "1" = "CD4 T", "2" = "CD4 T", "3" = "CD4 T", "4" = "CD8 T", "5" = "CD8 T", 
                           "6" = "CD8 T", "7" = "DP T", "8" = "CD8 T", "9" = "CD8 T", "10" = "CD8 T", "11" = "DP T")

saveRDS(Cluster_35, file = "TMS_merge_T_cells_Cluster_35_clustered.rds")
Cluster_35 <- readRDS(file = "TMS_merge_T_cells_Cluster_35_clustered.rds")

