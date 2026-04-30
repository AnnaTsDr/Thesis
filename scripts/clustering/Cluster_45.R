Cluster_45 <- readRDS(file = "TMS_merge_T_cells_clusters_45.rds")

#Highly variable expressed genes
Cluster_45 <- FindVariableFeatures(Cluster_45, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_45),15)
plot1 <- VariableFeaturePlot(Cluster_45)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_45), file = "TMS_merge_T_cells_Cluster_45_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_45, file = "TMS_merge_T_cells_Cluster_45_VariableGenes.rds")
Cluster_45 <- readRDS(file = "TMS_merge_T_cells_Cluster_45_VariableGenes.rds")

all.genes <- rownames(Cluster_45)
Cluster_45 <- ScaleData(Cluster_45, features = all.genes)

saveRDS(Cluster_45, file = "TMS_merge_T_cells_Cluster_45_ScaledData.rds")
Cluster_45 <- readRDS(file = "TMS_merge_T_cells_Cluster_45_ScaledData.rds")

Cluster_45 <- RunPCA(Cluster_45, npcs = 3, ndims.print = 1:3, nfeatures.print = 5)
DimPlot(Cluster_45, reduction = "pca")

saveRDS(Cluster_45, file = "TMS_merge_T_cells_Cluster_45_PCA.rds")
Cluster_45 <- readRDS(file = "TMS_merge_T_cells_Cluster_45_PCA.rds")

Cluster_45 <- JackStraw(Cluster_45, num.replicate = 100, dims = 3)
Cluster_45 <- ScoreJackStraw(Cluster_45, dims = 1:3)
JackStrawPlot(Cluster_45, dims = 1:3)

ElbowPlot(Cluster_45, ndims = 3)

Cluster_45 <- RunTSNE(Cluster_45, dims = 1:3, perplexity = 10)
DimPlot(Cluster_45, reduction = "tsne")

Cluster_45 <- RunUMAP(Cluster_45, dims = 1:3)
DimPlot(Cluster_45, reduction = "umap")

saveRDS(Cluster_45, file = "TMS_merge_T_cells_Cluster_45_TSNE_UMAP.rds")
Cluster_45 <- readRDS(file = "TMS_merge_T_cells_Cluster_45_TSNE_UMAP.rds")

Cluster_45 <- FindNeighbors(Cluster_45, dims = 1:3)

#louvian
Cluster_45 <- FindClusters(Cluster_45, resolution = seq(0,3, by = 0.1))

clustree(Cluster_45, prefix = "RNA_snn_res.")
clustree(Cluster_45, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_45) <- Cluster_45$RNA_snn_res.3

Cluster_45 <- BuildClusterTree(Cluster_45, reorder.numeric = TRUE, reorder = TRUE, dims = 1:3)
PlotClusterTree(object = Cluster_45)

Cluster_45.markers <- FindAllMarkers(Cluster_45)
m <- Cluster_45.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_45, features = m$gene)
write.csv(Cluster_45.markers, "TMS_T_cells_Cluster_45_merkers.csv")

Cluster_45 <- RenameIdents(Cluster_45, "1" = "CD8 T", "2" = "DN T", "3" = "DN T", "4" = "DN T", "5" = "CD8 T",
                           "6" = "CD4 T", "7" = "CD8 T", "8" = "CD4 T", "9" = "CD8 T", "10" = "DN T", "11" = "DN T")

saveRDS(Cluster_45, file = "TMS_merge_T_cells_Cluster_45_clustered.rds")
Cluster_45 <- readRDS(file = "TMS_merge_T_cells_Cluster_45_clustered.rds")
