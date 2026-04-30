Cluster_44 <- readRDS(file = "TMS_merge_T_cells_clusters_44.rds")

#Highly variable expressed genes
Cluster_44 <- FindVariableFeatures(Cluster_44, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_44),15)
plot1 <- VariableFeaturePlot(Cluster_44)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_44), file = "TMS_merge_T_cells_Cluster_44_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_44, file = "TMS_merge_T_cells_Cluster_44_VariableGenes.rds")
Cluster_44 <- readRDS(file = "TMS_merge_T_cells_Cluster_44_VariableGenes.rds")

all.genes <- rownames(Cluster_44)
Cluster_44 <- ScaleData(Cluster_44, features = all.genes)

saveRDS(Cluster_44, file = "TMS_merge_T_cells_Cluster_44_ScaledData.rds")
Cluster_44 <- readRDS(file = "TMS_merge_T_cells_Cluster_44_ScaledData.rds")

Cluster_44 <- RunPCA(Cluster_44, npcs = 5, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(Cluster_44, reduction = "pca")

saveRDS(Cluster_44, file = "TMS_merge_T_cells_Cluster_44_PCA.rds")
Cluster_44 <- readRDS(file = "TMS_merge_T_cells_Cluster_44_PCA.rds")

Cluster_44 <- JackStraw(Cluster_44, num.replicate = 100, dims = 5)
Cluster_44 <- ScoreJackStraw(Cluster_44, dims = 1:5)
JackStrawPlot(Cluster_44, dims = 1:5)

ElbowPlot(Cluster_44, ndims = 5)

Cluster_44 <- RunTSNE(Cluster_44, dims = 1:5, perplexity = 10)
DimPlot(Cluster_44, reduction = "tsne")

Cluster_44 <- RunUMAP(Cluster_44, dims = 1:5)
DimPlot(Cluster_44, reduction = "umap")

saveRDS(Cluster_44, file = "TMS_merge_T_cells_Cluster_44_TSNE_UMAP.rds")
Cluster_44 <- readRDS(file = "TMS_merge_T_cells_Cluster_44_TSNE_UMAP.rds")

Cluster_44 <- FindNeighbors(Cluster_44, dims = 1:5)

#louvian
Cluster_44 <- FindClusters(Cluster_44, resolution = seq(0,3, by = 0.1))

clustree(Cluster_44, prefix = "RNA_snn_res.")
clustree(Cluster_44, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_44) <- Cluster_44$RNA_snn_res.3

Cluster_44 <- BuildClusterTree(Cluster_44, reorder.numeric = TRUE, reorder = TRUE, dims = 1:5)
PlotClusterTree(object = Cluster_44)

Cluster_44.markers <- FindAllMarkers(Cluster_44)
m <- Cluster_44.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_44, features = m$gene)
write.csv(Cluster_44.markers, "TMS_T_cells_Cluster_44_merkers.csv")

Cluster_44 <- RenameIdents(Cluster_44, "1" = "CD8 T", "2" = "B cell", "3" = "B cell", "4" = "B cell", "5" = "B cell",
                           "6" = "B cell", "7" = "CD4 T", "8" = "CD8 T", "9" = "CD4 T", "10" = "CD4 T", "11" = "DP T",
                           "12" = "CD8 T", "13" = "CD8 T", "14" = "CD8 T", "15" = "CD8 T", "16" = "CD8 T")

saveRDS(Cluster_44, file = "TMS_merge_T_cells_Cluster_44_clustered.rds")
Cluster_44 <- readRDS(file = "TMS_merge_T_cells_Cluster_44_clustered.rds")

