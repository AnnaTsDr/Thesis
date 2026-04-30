Cluster_64 <- readRDS(file = "TMS_merge_T_cells_clusters_64.rds")

#Highly variable expressed genes
Cluster_64 <- FindVariableFeatures(Cluster_64, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_64),15)
plot1 <- VariableFeaturePlot(Cluster_64)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_64), file = "TMS_merge_T_cells_Cluster_64_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_64, file = "TMS_merge_T_cells_Cluster_64_VariableGenes.rds")
Cluster_64 <- readRDS(file = "TMS_merge_T_cells_Cluster_64_VariableGenes.rds")

all.genes <- rownames(Cluster_64)
Cluster_64 <- ScaleData(Cluster_64, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_64, file = "TMS_merge_T_cells_Cluster_64_ScaledData.rds")
Cluster_64 <- readRDS(file = "TMS_merge_T_cells_Cluster_64_ScaledData.rds")

Cluster_64 <- RunPCA(Cluster_64, ndims.print = 1:11, nfeatures.print = 5)
DimPlot(Cluster_64, reduction = "pca")

saveRDS(Cluster_64, file = "TMS_merge_T_cells_Cluster_64_PCA.rds")
Cluster_64 <- readRDS(file = "TMS_merge_T_cells_Cluster_64_PCA.rds")

Cluster_64 <- JackStraw(Cluster_64, num.replicate = 100, dims = 11)
Cluster_64 <- ScoreJackStraw(Cluster_64, dims = 1:11)
JackStrawPlot(Cluster_64, dims = 1:11)

ElbowPlot(Cluster_64, ndims = 11)

Cluster_64 <- RunTSNE(Cluster_64, dims = 1:11, perplexity = 10)
DimPlot(Cluster_64, reduction = "tsne")

Cluster_64 <- RunUMAP(Cluster_64, dims = 1:11)
DimPlot(Cluster_64, reduction = "umap")

saveRDS(Cluster_64, file = "TMS_merge_T_cells_Cluster_64_TSNE_UMAP.rds")
Cluster_64 <- readRDS(file = "TMS_merge_T_cells_Cluster_64_TSNE_UMAP.rds")

Cluster_64 <- FindNeighbors(Cluster_64, dims = 1:11)

#louvian
Cluster_64 <- FindClusters(Cluster_64, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_64, prefix = "RNA_snn_res.")
clustree(Cluster_64, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_64) <- Cluster_64$RNA_snn_res.1

Cluster_64 <- BuildClusterTree(Cluster_64, reorder.numeric = TRUE, reorder = TRUE, dims = 1:11)
PlotClusterTree(object = Cluster_64)

Cluster_64.markers <- FindAllMarkers(Cluster_64)
m <- Cluster_64.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_64, features = m$gene)
write.csv(Cluster_64.markers, "TMS_T_cells_Cluster_64_merkers.csv")

Cluster_64 <- RenameIdents(Cluster_64, "1" = "CD8 T", "2" = "CD4 T", "3" = "DP T", "4" = "CD8 T")

saveRDS(Cluster_64, file = "TMS_merge_T_cells_Cluster_64_clustered.rds")
Cluster_64 <- readRDS(file = "TMS_merge_T_cells_Cluster_64_clustered.rds")

