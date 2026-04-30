Cluster_34 <- readRDS(file = "TMS_merge_T_cells_clusters_34.rds")

#Highly variable expressed genes
Cluster_34 <- FindVariableFeatures(Cluster_34, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_34),15)
plot1 <- VariableFeaturePlot(Cluster_34)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_34), file = "TMS_merge_T_cells_Cluster_34_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_34, file = "TMS_merge_T_cells_Cluster_34_VariableGenes.rds")
Cluster_34 <- readRDS(file = "TMS_merge_T_cells_Cluster_34_VariableGenes.rds")

all.genes <- rownames(Cluster_34)
Cluster_34 <- ScaleData(Cluster_34, features = all.genes)

saveRDS(Cluster_34, file = "TMS_merge_T_cells_Cluster_34_ScaledData.rds")
Cluster_34 <- readRDS(file = "TMS_merge_T_cells_Cluster_34_ScaledData.rds")

Cluster_34 <- RunPCA(Cluster_34, npcs = 6, ndims.print = 1:6, nfeatures.print = 5)
DimPlot(Cluster_34, reduction = "pca")

saveRDS(Cluster_34, file = "TMS_merge_T_cells_Cluster_34_PCA.rds")
Cluster_34 <- readRDS(file = "TMS_merge_T_cells_Cluster_34_PCA.rds")

Cluster_34 <- JackStraw(Cluster_34, num.replicate = 100, dims = 6)
Cluster_34 <- ScoreJackStraw(Cluster_34, dims = 1:6)
JackStrawPlot(Cluster_34, dims = 1:6)

ElbowPlot(Cluster_34, ndims = 6)

Cluster_34 <- RunTSNE(Cluster_34, dims = 1:6, perplexity = 10)
DimPlot(Cluster_34, reduction = "tsne")

Cluster_34 <- RunUMAP(Cluster_34, dims = 1:6)
DimPlot(Cluster_34, reduction = "umap")

saveRDS(Cluster_34, file = "TMS_merge_T_cells_Cluster_34_TSNE_UMAP.rds")
Cluster_34 <- readRDS(file = "TMS_merge_T_cells_Cluster_34_TSNE_UMAP.rds")

Cluster_34 <- FindNeighbors(Cluster_34, dims = 1:6)

#louvian
Cluster_34 <- FindClusters(Cluster_34, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_34, prefix = "RNA_snn_res.")
clustree(Cluster_34, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_34) <- Cluster_34$RNA_snn_res.0

#Cluster_34 <- BuildClusterTree(Cluster_34, reorder.numeric = TRUE, reorder = TRUE, dims = 1:6)
#PlotClusterTree(object = Cluster_34)

#Cluster_34.markers <- FindAllMarkers(Cluster_34)
#m <- Cluster_34.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_34, features = m$gene)
#write.csv(Cluster_34.markers, "TMS_T_cells_Cluster_34_merkers.csv")

Cluster_34 <- RenameIdents(Cluster_34, "0" = "CD8 T")

saveRDS(Cluster_34, file = "TMS_merge_T_cells_Cluster_34_clustered.rds")
Cluster_34 <- readRDS(file = "TMS_merge_T_cells_Cluster_34_clustered.rds")

