Cluster_65 <- readRDS(file = "TMS_merge_T_cells_clusters_65.rds")

#Highly variable expressed genes
Cluster_65 <- FindVariableFeatures(Cluster_65, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_65),15)
plot1 <- VariableFeaturePlot(Cluster_65)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_65), file = "TMS_merge_T_cells_Cluster_65_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_65, file = "TMS_merge_T_cells_Cluster_65_VariableGenes.rds")
Cluster_65 <- readRDS(file = "TMS_merge_T_cells_Cluster_65_VariableGenes.rds")

all.genes <- rownames(Cluster_65)
Cluster_65 <- ScaleData(Cluster_65, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_65, file = "TMS_merge_T_cells_Cluster_65_ScaledData.rds")
Cluster_65 <- readRDS(file = "TMS_merge_T_cells_Cluster_65_ScaledData.rds")

Cluster_65 <- RunPCA(Cluster_65, ndims.print = 1:3, nfeatures.print = 5)
DimPlot(Cluster_65, reduction = "pca")

saveRDS(Cluster_65, file = "TMS_merge_T_cells_Cluster_65_PCA.rds")
Cluster_65 <- readRDS(file = "TMS_merge_T_cells_Cluster_65_PCA.rds")

Cluster_65 <- JackStraw(Cluster_65, num.replicate = 100, dims = 3)
Cluster_65 <- ScoreJackStraw(Cluster_65, dims = 1:3)
JackStrawPlot(Cluster_65, dims = 1:3)

ElbowPlot(Cluster_65, ndims = 3)

Cluster_65 <- RunTSNE(Cluster_65, dims = 1:3, perplexity = 10)
DimPlot(Cluster_65, reduction = "tsne")

Cluster_65 <- RunUMAP(Cluster_65, dims = 1:3)
DimPlot(Cluster_65, reduction = "umap")

saveRDS(Cluster_65, file = "TMS_merge_T_cells_Cluster_65_TSNE_UMAP.rds")
Cluster_65 <- readRDS(file = "TMS_merge_T_cells_Cluster_65_TSNE_UMAP.rds")

Cluster_65 <- FindNeighbors(Cluster_65, dims = 1:3)

#louvian
Cluster_65 <- FindClusters(Cluster_65, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_65, prefix = "RNA_snn_res.")
clustree(Cluster_65, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_65) <- Cluster_65$RNA_snn_res.0.9

Cluster_65 <- BuildClusterTree(Cluster_65, reorder.numeric = TRUE, reorder = TRUE, dims = 1:3)
PlotClusterTree(object = Cluster_65)

Cluster_65.markers <- FindAllMarkers(Cluster_65)
m <- Cluster_65.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_65, features = m$gene)
write.csv(Cluster_65.markers, "TMS_T_cells_Cluster_65_merkers.csv")

Cluster_65 <- RenameIdents(Cluster_65, "1" = "DN T", "2" = "DN T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD8 T")

saveRDS(Cluster_65, file = "TMS_merge_T_cells_Cluster_65_clustered.rds")
Cluster_65 <- readRDS(file = "TMS_merge_T_cells_Cluster_65_clustered.rds")

