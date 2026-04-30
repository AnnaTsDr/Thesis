Cluster_4 <- readRDS(file = "TMS_merge_T_cells_clusters_4.rds")

#Highly variable expressed genes
Cluster_4 <- FindVariableFeatures(Cluster_4, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_4),15)
plot1 <- VariableFeaturePlot(Cluster_4)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_4), file = "TMS_merge_T_cells_Cluster_4_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_4, file = "TMS_merge_T_cells_Cluster_4_VariableGenes.rds")
Cluster_4 <- readRDS(file = "TMS_merge_T_cells_Cluster_4_VariableGenes.rds")

all.genes <- rownames(Cluster_4)
Cluster_4 <- ScaleData(Cluster_4, features = all.genes)

saveRDS(Cluster_4, file = "TMS_merge_T_cells_Cluster_4_ScaledData.rds")
Cluster_4 <- readRDS(file = "TMS_merge_T_cells_Cluster_4_ScaledData.rds")

Cluster_4 <- RunPCA(Cluster_4, npcs = 6, ndims.print = 1:6, nfeatures.print = 5)
DimPlot(Cluster_4, reduction = "pca")

saveRDS(Cluster_4, file = "TMS_merge_T_cells_Cluster_4_PCA.rds")
Cluster_4 <- readRDS(file = "TMS_merge_T_cells_Cluster_4_PCA.rds")

Cluster_4 <- JackStraw(Cluster_4, num.replicate = 100, dims = 6)
Cluster_4 <- ScoreJackStraw(Cluster_4, dims = 1:6)
JackStrawPlot(Cluster_4, dims = 1:6)

ElbowPlot(Cluster_4, ndims = 6)

Cluster_4 <- RunTSNE(Cluster_4, dims = 1:6, perplexity = 10)
DimPlot(Cluster_4, reduction = "tsne")

Cluster_4 <- RunUMAP(Cluster_4, dims = 1:6)
DimPlot(Cluster_4, reduction = "umap")

saveRDS(Cluster_4, file = "TMS_merge_T_cells_Cluster_4_TSNE_UMAP.rds")
Cluster_4 <- readRDS(file = "TMS_merge_T_cells_Cluster_4_TSNE_UMAP.rds")

Cluster_4 <- FindNeighbors(Cluster_4, dims = 1:6)

#louvian
Cluster_4 <- FindClusters(Cluster_4, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_4, prefix = "RNA_snn_res.")
clustree(Cluster_4, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_4) <- Cluster_4$RNA_snn_res.1.2

Cluster_4 <- BuildClusterTree(Cluster_4, reorder.numeric = TRUE, reorder = TRUE, dims = 1:6)
PlotClusterTree(object = Cluster_4)

Cluster_4.markers <- FindAllMarkers(Cluster_4)
m <- Cluster_4.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_4, features = m$gene)
write.csv(Cluster_4.markers, "TMS_T_cells_Cluster_4_merkers.csv")

Cluster_4 <- RenameIdents(Cluster_4, "1" = "Not T", "2" = "Not T", "3" = "Not T", "4" = "CD4 T", "5" = "CD4 T", 
                          "6" = "CD4 T")

saveRDS(Cluster_4, file = "TMS_merge_T_cells_Cluster_4_clustered.rds")
Cluster_4 <- readRDS(file = "TMS_merge_T_cells_Cluster_4_clustered.rds")

