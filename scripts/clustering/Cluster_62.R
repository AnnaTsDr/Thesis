Cluster_62 <- readRDS(file = "TMS_merge_T_cells_clusters_62.rds")

#Highly variable expressed genes
Cluster_62 <- FindVariableFeatures(Cluster_62, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_62),15)
plot1 <- VariableFeaturePlot(Cluster_62)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_62), file = "TMS_merge_T_cells_Cluster_62_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_62, file = "TMS_merge_T_cells_Cluster_62_VariableGenes.rds")
Cluster_62 <- readRDS(file = "TMS_merge_T_cells_Cluster_62_VariableGenes.rds")

all.genes <- rownames(Cluster_62)
Cluster_62 <- ScaleData(Cluster_62, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_62, file = "TMS_merge_T_cells_Cluster_62_ScaledData.rds")
Cluster_62 <- readRDS(file = "TMS_merge_T_cells_Cluster_62_ScaledData.rds")

Cluster_62 <- RunPCA(Cluster_62, ndims.print = 1:8, nfeatures.print = 5)
DimPlot(Cluster_62, reduction = "pca")

saveRDS(Cluster_62, file = "TMS_merge_T_cells_Cluster_62_PCA.rds")
Cluster_62 <- readRDS(file = "TMS_merge_T_cells_Cluster_62_PCA.rds")

Cluster_62 <- JackStraw(Cluster_62, num.replicate = 100, dims = 8)
Cluster_62 <- ScoreJackStraw(Cluster_62, dims = 1:8)
JackStrawPlot(Cluster_62, dims = 1:8)

ElbowPlot(Cluster_62, ndims = 8)

Cluster_62 <- RunTSNE(Cluster_62, dims = 1:8, perplexity = 10)
DimPlot(Cluster_62, reduction = "tsne")

Cluster_62 <- RunUMAP(Cluster_62, dims = 1:8)
DimPlot(Cluster_62, reduction = "umap")

saveRDS(Cluster_62, file = "TMS_merge_T_cells_Cluster_62_TSNE_UMAP.rds")
Cluster_62 <- readRDS(file = "TMS_merge_T_cells_Cluster_62_TSNE_UMAP.rds")

Cluster_62 <- FindNeighbors(Cluster_62, dims = 1:8)

#louvian
Cluster_62 <- FindClusters(Cluster_62, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_62, prefix = "RNA_snn_res.")
clustree(Cluster_62, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_62) <- Cluster_62$RNA_snn_res.0.5

Cluster_62 <- BuildClusterTree(Cluster_62, reorder.numeric = TRUE, reorder = TRUE, dims = 1:8)
PlotClusterTree(object = Cluster_62)

Cluster_62.markers <- FindAllMarkers(Cluster_62)
m <- Cluster_62.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_62, features = m$gene)
write.csv(Cluster_62.markers, "TMS_T_cells_Cluster_62_merkers.csv")

Cluster_62 <- RenameIdents(Cluster_62, "1" = "DN T", "2" = "CD4 T", "3" = "DP T", "4" = "CD8 T")

saveRDS(Cluster_62, file = "TMS_merge_T_cells_Cluster_62_clustered.rds")
Cluster_62 <- readRDS(file = "TMS_merge_T_cells_Cluster_62_clustered.rds")

