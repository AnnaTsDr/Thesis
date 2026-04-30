Cluster_22 <- readRDS(file = "TMS_merge_T_cells_clusters_22.rds")

#Highly variable expressed genes
Cluster_22 <- FindVariableFeatures(Cluster_22, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_22),15)
plot1 <- VariableFeaturePlot(Cluster_22)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_22), file = "TMS_merge_T_cells_Cluster_22_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_22, file = "TMS_merge_T_cells_Cluster_22_VariableGenes.rds")
Cluster_22 <- readRDS(file = "TMS_merge_T_cells_Cluster_22_VariableGenes.rds")

all.genes <- rownames(Cluster_22)
Cluster_22 <- ScaleData(Cluster_22, features = all.genes)

saveRDS(Cluster_22, file = "TMS_merge_T_cells_Cluster_22_ScaledData.rds")
Cluster_22 <- readRDS(file = "TMS_merge_T_cells_Cluster_22_ScaledData.rds")

Cluster_22 <- RunPCA(Cluster_22, npcs = 16, ndims.print = 1:16, nfeatures.print = 5)
DimPlot(Cluster_22, reduction = "pca")

saveRDS(Cluster_22, file = "TMS_merge_T_cells_Cluster_22_PCA.rds")
Cluster_22 <- readRDS(file = "TMS_merge_T_cells_Cluster_22_PCA.rds")

Cluster_22 <- JackStraw(Cluster_22, num.replicate = 100, dims = 16)
Cluster_22 <- ScoreJackStraw(Cluster_22, dims = 1:16)
JackStrawPlot(Cluster_22, dims = 1:16)

ElbowPlot(Cluster_22, ndims = 16)

Cluster_22 <- RunTSNE(Cluster_22, dims = 1:16, perplexity = 10)
DimPlot(Cluster_22, reduction = "tsne")

Cluster_22 <- RunUMAP(Cluster_22, dims = 1:16)
DimPlot(Cluster_22, reduction = "umap")

saveRDS(Cluster_22, file = "TMS_merge_T_cells_Cluster_22_TSNE_UMAP.rds")
Cluster_22 <- readRDS(file = "TMS_merge_T_cells_Cluster_22_TSNE_UMAP.rds")

Cluster_22 <- FindNeighbors(Cluster_22, dims = 1:16)

#louvian
Cluster_22 <- FindClusters(Cluster_22, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_22, prefix = "RNA_snn_res.")
clustree(Cluster_22, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_22) <- Cluster_22$RNA_snn_res.0

#Cluster_22 <- BuildClusterTree(Cluster_22, reorder.numeric = TRUE, reorder = TRUE, dims = 1:16)
#PlotClusterTree(object = Cluster_22)

#Cluster_22.markers <- FindAllMarkers(Cluster_22)
#m <- Cluster_22.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_22, features = m$gene)
#write.csv(Cluster_22.markers, "TMS_T_cells_Cluster_22_merkers.csv")

Cluster_22 <- RenameIdents(Cluster_22, "0" = "CD8 T")

saveRDS(Cluster_22, file = "TMS_merge_T_cells_Cluster_22_clustered.rds")
Cluster_22 <- readRDS(file = "TMS_merge_T_cells_Cluster_22_clustered.rds")

