Cluster_38 <- readRDS(file = "TMS_merge_T_cells_clusters_38.rds")

#Highly variable expressed genes
Cluster_38 <- FindVariableFeatures(Cluster_38, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.8, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_38),15)
plot1 <- VariableFeaturePlot(Cluster_38)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_38), file = "TMS_merge_T_cells_Cluster_38_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_38, file = "TMS_merge_T_cells_Cluster_38_VariableGenes.rds")
Cluster_38 <- readRDS(file = "TMS_merge_T_cells_Cluster_38_VariableGenes.rds")

all.genes <- rownames(Cluster_38)
Cluster_38 <- ScaleData(Cluster_38, features = all.genes)

saveRDS(Cluster_38, file = "TMS_merge_T_cells_Cluster_38_ScaledData.rds")
Cluster_38 <- readRDS(file = "TMS_merge_T_cells_Cluster_38_ScaledData.rds")

Cluster_38 <- RunPCA(Cluster_38, npcs = 13, ndims.print = 1:13, nfeatures.print = 5)
DimPlot(Cluster_38, reduction = "pca")

saveRDS(Cluster_38, file = "TMS_merge_T_cells_Cluster_38_PCA.rds")
Cluster_38 <- readRDS(file = "TMS_merge_T_cells_Cluster_38_PCA.rds")

Cluster_38 <- JackStraw(Cluster_38, num.replicate = 100, dims = 13)
Cluster_38 <- ScoreJackStraw(Cluster_38, dims = 1:13)
JackStrawPlot(Cluster_38, dims = 1:13)

ElbowPlot(Cluster_38, ndims = 13)

Cluster_38 <- RunTSNE(Cluster_38, dims = 1:13, perplexity = 10)
DimPlot(Cluster_38, reduction = "tsne")

Cluster_38 <- RunUMAP(Cluster_38, dims = 1:13)
DimPlot(Cluster_38, reduction = "umap")

saveRDS(Cluster_38, file = "TMS_merge_T_cells_Cluster_38_TSNE_UMAP.rds")
Cluster_38 <- readRDS(file = "TMS_merge_T_cells_Cluster_38_TSNE_UMAP.rds")

Cluster_38 <- FindNeighbors(Cluster_38, dims = 1:13)

#louvian
Cluster_38 <- FindClusters(Cluster_38, resolution = seq(0,3, by = 0.1))

clustree(Cluster_38, prefix = "RNA_snn_res.")
clustree(Cluster_38, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_38) <- Cluster_38$RNA_snn_res.3

Cluster_38 <- BuildClusterTree(Cluster_38, reorder.numeric = TRUE, reorder = TRUE, dims = 1:13)
PlotClusterTree(object = Cluster_38)

Cluster_38.markers <- FindAllMarkers(Cluster_38)
m <- Cluster_38.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_38, features = m$gene)
write.csv(Cluster_38.markers, "TMS_T_cells_Cluster_38_merkers.csv")

Cluster_38 <- RenameIdents(Cluster_38, "1" = "CD8 T", "2" = "CD8 T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD8 T", 
                           "6" = "CD8 T", "7" = "CD4 T", "8" = "CD4 T", "9" = "CD4 T", "10" = "DP T", "11" = "CD4 T",
                           "12" = "CD8 T", "13" = "CD4 T", "14" = "CD8 T", "15" = "CD8 T")

saveRDS(Cluster_38, file = "TMS_merge_T_cells_Cluster_38_clustered.rds")
Cluster_38 <- readRDS(file = "TMS_merge_T_cells_Cluster_38_clustered.rds")
