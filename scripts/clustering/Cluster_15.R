Cluster_15 <- readRDS(file = "TMS_merge_T_cells_clusters_15.rds")

#Highly variable expressed genes
Cluster_15 <- FindVariableFeatures(Cluster_15, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_15),15)
plot1 <- VariableFeaturePlot(Cluster_15)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_15), file = "TMS_merge_T_cells_Cluster_15_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_15, file = "TMS_merge_T_cells_Cluster_15_VariableGenes.rds")
Cluster_15 <- readRDS(file = "TMS_merge_T_cells_Cluster_15_VariableGenes.rds")

all.genes <- rownames(Cluster_15)
Cluster_15 <- ScaleData(Cluster_15, features = all.genes)

saveRDS(Cluster_15, file = "TMS_merge_T_cells_Cluster_15_ScaledData.rds")
Cluster_15 <- readRDS(file = "TMS_merge_T_cells_Cluster_15_ScaledData.rds")

Cluster_15 <- RunPCA(Cluster_15, npcs = 3, ndims.print = 1:3, nfeatures.print = 5)
DimPlot(Cluster_15, reduction = "pca")

saveRDS(Cluster_15, file = "TMS_merge_T_cells_Cluster_15_PCA.rds")
Cluster_15 <- readRDS(file = "TMS_merge_T_cells_Cluster_15_PCA.rds")

Cluster_15 <- JackStraw(Cluster_15, num.replicate = 100, dims = 3)
Cluster_15 <- ScoreJackStraw(Cluster_15, dims = 1:3)
JackStrawPlot(Cluster_15, dims = 1:3)

ElbowPlot(Cluster_15, ndims = 3)

Cluster_15 <- RunTSNE(Cluster_15, dims = 1:3, perplexity = 10)
DimPlot(Cluster_15, reduction = "tsne")

Cluster_15 <- RunUMAP(Cluster_15, dims = 1:3)
DimPlot(Cluster_15, reduction = "umap")

saveRDS(Cluster_15, file = "TMS_merge_T_cells_Cluster_15_TSNE_UMAP.rds")
Cluster_15 <- readRDS(file = "TMS_merge_T_cells_Cluster_15_TSNE_UMAP.rds")

Cluster_15 <- FindNeighbors(Cluster_15, dims = 1:3)

#louvian
Cluster_15 <- FindClusters(Cluster_15, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_15, prefix = "RNA_snn_res.")
clustree(Cluster_15, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_15) <- Cluster_15$RNA_snn_res.0

#Cluster_15 <- BuildClusterTree(Cluster_15, reorder.numeric = TRUE, reorder = TRUE, dims = 1:3)
#PlotClusterTree(object = Cluster_15)

#Cluster_15.markers <- FindAllMarkers(Cluster_15)
#m <- Cluster_15.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_15, features = m$gene)
#write.csv(Cluster_15.markers, "TMS_T_cells_Cluster_15_merkers.csv")

Cluster_15 <- RenameIdents(Cluster_15, "0" = "CD8 T")

saveRDS(Cluster_15, file = "TMS_merge_T_cells_Cluster_15_clustered.rds")
Cluster_15 <- readRDS(file = "TMS_merge_T_cells_Cluster_15_clustered.rds")

