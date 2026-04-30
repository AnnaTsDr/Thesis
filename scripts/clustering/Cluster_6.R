Cluster_6 <- readRDS(file = "TMS_merge_T_cells_clusters_6.rds")

#Highly variable expressed genes
Cluster_6 <- FindVariableFeatures(Cluster_6, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.5, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_6),15)
plot1 <- VariableFeaturePlot(Cluster_6)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_6), file = "TMS_merge_T_cells_Cluster_6_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_6, file = "TMS_merge_T_cells_Cluster_6_VariableGenes.rds")
Cluster_6 <- readRDS(file = "TMS_merge_T_cells_Cluster_6_VariableGenes.rds")

all.genes <- rownames(Cluster_6)
Cluster_6 <- ScaleData(Cluster_6, features = all.genes)

saveRDS(Cluster_6, file = "TMS_merge_T_cells_Cluster_6_ScaledData.rds")
Cluster_6 <- readRDS(file = "TMS_merge_T_cells_Cluster_6_ScaledData.rds")

Cluster_6 <- RunPCA(Cluster_6, npcs = 2, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_6, reduction = "pca")

saveRDS(Cluster_6, file = "TMS_merge_T_cells_Cluster_6_PCA.rds")
Cluster_6 <- readRDS(file = "TMS_merge_T_cells_Cluster_6_PCA.rds")

Cluster_6 <- JackStraw(Cluster_6, num.replicate = 100, dims = 2)
Cluster_6 <- ScoreJackStraw(Cluster_6, dims = 1:2)
JackStrawPlot(Cluster_6, dims = 1:2)

ElbowPlot(Cluster_6, ndims = 2)

Cluster_6 <- RunTSNE(Cluster_6, dims = 1:2, perplexity = 10)
DimPlot(Cluster_6, reduction = "tsne")

Cluster_6 <- RunUMAP(Cluster_6, dims = 1:2)
DimPlot(Cluster_6, reduction = "umap")

saveRDS(Cluster_6, file = "TMS_merge_T_cells_Cluster_6_TSNE_UMAP.rds")
Cluster_6 <- readRDS(file = "TMS_merge_T_cells_Cluster_6_TSNE_UMAP.rds")

Cluster_6 <- FindNeighbors(Cluster_6, dims = 1:2)

#louvian
Cluster_6 <- FindClusters(Cluster_6, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_6, prefix = "RNA_snn_res.")
clustree(Cluster_6, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_6) <- Cluster_6$RNA_snn_res.0

#Cluster_6 <- BuildClusterTree(Cluster_6, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
#PlotClusterTree(object = Cluster_6)

#Cluster_6.markers <- FindAllMarkers(Cluster_6)
#m <- Cluster_6.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_6, features = m$gene)
#write.csv(Cluster_6.markers, "TMS_T_cells_Cluster_6_merkers.csv")

Cluster_6 <- RenameIdents(Cluster_6, "0" = "CD4 T")

saveRDS(Cluster_6, file = "TMS_merge_T_cells_Cluster_6_clustered.rds")
Cluster_6 <- readRDS(file = "TMS_merge_T_cells_Cluster_6_clustered.rds")

#homogenize cluster