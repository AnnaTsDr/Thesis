Cluster_24 <- readRDS(file = "TMS_merge_T_cells_clusters_24.rds")

#Highly variable expressed genes
Cluster_24 <- FindVariableFeatures(Cluster_24, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.4, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_24),15)
plot1 <- VariableFeaturePlot(Cluster_24)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_24), file = "TMS_merge_T_cells_Cluster_24_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_24, file = "TMS_merge_T_cells_Cluster_24_VariableGenes.rds")
Cluster_24 <- readRDS(file = "TMS_merge_T_cells_Cluster_24_VariableGenes.rds")

all.genes <- rownames(Cluster_24)
Cluster_24 <- ScaleData(Cluster_24, features = all.genes)

saveRDS(Cluster_24, file = "TMS_merge_T_cells_Cluster_24_ScaledData.rds")
Cluster_24 <- readRDS(file = "TMS_merge_T_cells_Cluster_24_ScaledData.rds")

Cluster_24 <- RunPCA(Cluster_24, npcs = 2, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_24, reduction = "pca")

saveRDS(Cluster_24, file = "TMS_merge_T_cells_Cluster_24_PCA.rds")
Cluster_24 <- readRDS(file = "TMS_merge_T_cells_Cluster_24_PCA.rds")

Cluster_24 <- JackStraw(Cluster_24, num.replicate = 100, dims = 2)
Cluster_24 <- ScoreJackStraw(Cluster_24, dims = 1:2)
JackStrawPlot(Cluster_24, dims = 1:2)

ElbowPlot(Cluster_24, ndims = 2)

Cluster_24 <- RunTSNE(Cluster_24, dims = 1:2, perplexity = 10)
DimPlot(Cluster_24, reduction = "tsne")

Cluster_24 <- RunUMAP(Cluster_24, dims = 1:2)
DimPlot(Cluster_24, reduction = "umap")

saveRDS(Cluster_24, file = "TMS_merge_T_cells_Cluster_24_TSNE_UMAP.rds")
Cluster_24 <- readRDS(file = "TMS_merge_T_cells_Cluster_24_TSNE_UMAP.rds")

Cluster_24 <- FindNeighbors(Cluster_24, dims = 1:2)

#louvian
Cluster_24 <- FindClusters(Cluster_24, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_24, prefix = "RNA_snn_res.")
clustree(Cluster_24, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_24) <- Cluster_24$RNA_snn_res.0

#Cluster_24 <- BuildClusterTree(Cluster_24, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
#PlotClusterTree(object = Cluster_24)

#Cluster_24.markers <- FindAllMarkers(Cluster_24)
#m <- Cluster_24.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_24, features = m$gene)
#write.csv(Cluster_24.markers, "TMS_T_cells_Cluster_24_merkers.csv")

Cluster_24 <- RenameIdents(Cluster_24, "0" = "DN T")

saveRDS(Cluster_24, file = "TMS_merge_T_cells_Cluster_24_clustered.rds")
Cluster_24 <- readRDS(file = "TMS_merge_T_cells_Cluster_24_clustered.rds")

