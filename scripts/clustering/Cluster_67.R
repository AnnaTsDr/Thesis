Cluster_67 <- readRDS(file = "TMS_merge_T_cells_clusters_67.rds")

#Highly variable expressed genes
Cluster_67 <- FindVariableFeatures(Cluster_67, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.7, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_67),15)
plot1 <- VariableFeaturePlot(Cluster_67)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_67), file = "TMS_merge_T_cells_Cluster_67_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_67, file = "TMS_merge_T_cells_Cluster_67_VariableGenes.rds")
Cluster_67 <- readRDS(file = "TMS_merge_T_cells_Cluster_67_VariableGenes.rds")

all.genes <- rownames(Cluster_67)
Cluster_67 <- ScaleData(Cluster_67, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_67, file = "TMS_merge_T_cells_Cluster_67_ScaledData.rds")
Cluster_67 <- readRDS(file = "TMS_merge_T_cells_Cluster_67_ScaledData.rds")

Cluster_67 <- RunPCA(Cluster_67,npcs = 63, ndims.print = 1:63, nfeatures.print = 5)
DimPlot(Cluster_67, reduction = "pca")

saveRDS(Cluster_67, file = "TMS_merge_T_cells_Cluster_67_PCA.rds")
Cluster_67 <- readRDS(file = "TMS_merge_T_cells_Cluster_67_PCA.rds")

Cluster_67 <- JackStraw(Cluster_67, num.replicate = 100, dims = 63)
Cluster_67 <- ScoreJackStraw(Cluster_67, dims = 1:63)
JackStrawPlot(Cluster_67, dims = 1:63)

ElbowPlot(Cluster_67, ndims = 63)

Cluster_67 <- RunTSNE(Cluster_67, dims = 1:63, perplexity = 10)
DimPlot(Cluster_67, reduction = "tsne")

Cluster_67 <- RunUMAP(Cluster_67, dims = 1:63)
DimPlot(Cluster_67, reduction = "umap")

saveRDS(Cluster_67, file = "TMS_merge_T_cells_Cluster_67_TSNE_UMAP.rds")
Cluster_67 <- readRDS(file = "TMS_merge_T_cells_Cluster_67_TSNE_UMAP.rds")

Cluster_67 <- FindNeighbors(Cluster_67, dims = 1:63)

#louvian
Cluster_67 <- FindClusters(Cluster_67, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_67, prefix = "RNA_snn_res.")
clustree(Cluster_67, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_67) <- Cluster_67$RNA_snn_res.1

Cluster_67 <- BuildClusterTree(Cluster_67, reorder.numeric = TRUE, reorder = TRUE, dims = 1:63)
PlotClusterTree(object = Cluster_67)

Cluster_67.markers <- FindAllMarkers(Cluster_67)
m <- Cluster_67.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_67, features = m$gene)
write.csv(Cluster_67.markers, "TMS_T_cells_Cluster_67_merkers.csv")

Cluster_67 <- RenameIdents(Cluster_67, "1" = "DP T", "2" = "DN T", "3" = "DP T", "4" = "DP T", "5" = "DP T")

saveRDS(Cluster_67, file = "TMS_merge_T_cells_Cluster_67_clustered.rds")
Cluster_67 <- readRDS(file = "TMS_merge_T_cells_Cluster_67_clustered.rds")

