Cluster_61 <- readRDS(file = "TMS_merge_T_cells_clusters_61.rds")

#Highly variable expressed genes
Cluster_61 <- FindVariableFeatures(Cluster_61, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_61),15)
plot1 <- VariableFeaturePlot(Cluster_61)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_61), file = "TMS_merge_T_cells_Cluster_61_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_61, file = "TMS_merge_T_cells_Cluster_61_VariableGenes.rds")
Cluster_61 <- readRDS(file = "TMS_merge_T_cells_Cluster_61_VariableGenes.rds")

all.genes <- rownames(Cluster_61)
Cluster_61 <- ScaleData(Cluster_61, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_61, file = "TMS_merge_T_cells_Cluster_61_ScaledData.rds")
Cluster_61 <- readRDS(file = "TMS_merge_T_cells_Cluster_61_ScaledData.rds")

Cluster_61 <- RunPCA(Cluster_61, ndims.print = 1:11, nfeatures.print = 5)
DimPlot(Cluster_61, reduction = "pca")

saveRDS(Cluster_61, file = "TMS_merge_T_cells_Cluster_61_PCA.rds")
Cluster_61 <- readRDS(file = "TMS_merge_T_cells_Cluster_61_PCA.rds")

Cluster_61 <- JackStraw(Cluster_61, num.replicate = 100, dims = 11)
Cluster_61 <- ScoreJackStraw(Cluster_61, dims = 1:11)
JackStrawPlot(Cluster_61, dims = 1:11)

ElbowPlot(Cluster_61, ndims = 11)

Cluster_61 <- RunTSNE(Cluster_61, dims = 1:11, perplexity = 10)
DimPlot(Cluster_61, reduction = "tsne")

Cluster_61 <- RunUMAP(Cluster_61, dims = 1:11)
DimPlot(Cluster_61, reduction = "umap")

saveRDS(Cluster_61, file = "TMS_merge_T_cells_Cluster_61_TSNE_UMAP.rds")
Cluster_61 <- readRDS(file = "TMS_merge_T_cells_Cluster_61_TSNE_UMAP.rds")

Cluster_61 <- FindNeighbors(Cluster_61, dims = 1:11)

#louvian
Cluster_61 <- FindClusters(Cluster_61, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_61, prefix = "RNA_snn_res.")
clustree(Cluster_61, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_61) <- Cluster_61$RNA_snn_res.1.2

Cluster_61 <- BuildClusterTree(Cluster_61, reorder.numeric = TRUE, reorder = TRUE, dims = 1:11)
PlotClusterTree(object = Cluster_61)

Cluster_61.markers <- FindAllMarkers(Cluster_61)
m <- Cluster_61.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_61, features = m$gene)
write.csv(Cluster_61.markers, "TMS_T_cells_Cluster_61_merkers.csv")

Cluster_61 <- RenameIdents(Cluster_61, "1" = "DP T", "2" = "CD8 T", "3" = "CD4 T", "4" = "DN T", "5" = "DP T", 
                           "6" = "DP T", "7" = "CD8 T")

saveRDS(Cluster_61, file = "TMS_merge_T_cells_Cluster_61_clustered.rds")
Cluster_61 <- readRDS(file = "TMS_merge_T_cells_Cluster_61_clustered.rds")

