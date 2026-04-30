Cluster_33 <- readRDS(file = "TMS_merge_T_cells_clusters_33.rds")

#Highly variable expressed genes
Cluster_33 <- FindVariableFeatures(Cluster_33, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_33),15)
plot1 <- VariableFeaturePlot(Cluster_33)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_33), file = "TMS_merge_T_cells_Cluster_33_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_33, file = "TMS_merge_T_cells_Cluster_33_VariableGenes.rds")
Cluster_33 <- readRDS(file = "TMS_merge_T_cells_Cluster_33_VariableGenes.rds")

all.genes <- rownames(Cluster_33)
Cluster_33 <- ScaleData(Cluster_33, features = all.genes)

saveRDS(Cluster_33, file = "TMS_merge_T_cells_Cluster_33_ScaledData.rds")
Cluster_33 <- readRDS(file = "TMS_merge_T_cells_Cluster_33_ScaledData.rds")

Cluster_33 <- RunPCA(Cluster_33, npcs = 12, ndims.print = 1:12, nfeatures.print = 5)
DimPlot(Cluster_33, reduction = "pca")

saveRDS(Cluster_33, file = "TMS_merge_T_cells_Cluster_33_PCA.rds")
Cluster_33 <- readRDS(file = "TMS_merge_T_cells_Cluster_33_PCA.rds")

Cluster_33 <- JackStraw(Cluster_33, num.replicate = 100, dims = 12)
Cluster_33 <- ScoreJackStraw(Cluster_33, dims = 1:12)
JackStrawPlot(Cluster_33, dims = 1:12)

ElbowPlot(Cluster_33, ndims = 12)

Cluster_33 <- RunTSNE(Cluster_33, dims = 1:12, perplexity = 10)
DimPlot(Cluster_33, reduction = "tsne")

Cluster_33 <- RunUMAP(Cluster_33, dims = 1:12)
DimPlot(Cluster_33, reduction = "umap")

saveRDS(Cluster_33, file = "TMS_merge_T_cells_Cluster_33_TSNE_UMAP.rds")
Cluster_33 <- readRDS(file = "TMS_merge_T_cells_Cluster_33_TSNE_UMAP.rds")

Cluster_33 <- FindNeighbors(Cluster_33, dims = 1:12)

#louvian
Cluster_33 <- FindClusters(Cluster_33, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_33, prefix = "RNA_snn_res.")
clustree(Cluster_33, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_33) <- Cluster_33$RNA_snn_res.0.6

Cluster_33 <- BuildClusterTree(Cluster_33, reorder.numeric = TRUE, reorder = TRUE, dims = 1:12)
PlotClusterTree(object = Cluster_33)

Cluster_33.markers <- FindAllMarkers(Cluster_33)
m <- Cluster_33.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_33, features = m$gene)
write.csv(Cluster_33.markers, "TMS_T_cells_Cluster_33_merkers.csv")

Cluster_33 <- RenameIdents(Cluster_33, "1" = "CD4 T", "2" = "CD4 T", "3" = "CD8 T")

saveRDS(Cluster_33, file = "TMS_merge_T_cells_Cluster_33_clustered.rds")
Cluster_33 <- readRDS(file = "TMS_merge_T_cells_Cluster_33_clustered.rds")

