Cluster_19 <- readRDS(file = "TMS_merge_T_cells_clusters_19.rds")

#Highly variable expressed genes
Cluster_19 <- FindVariableFeatures(Cluster_19, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_19),15)
plot1 <- VariableFeaturePlot(Cluster_19)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_19), file = "TMS_merge_T_cells_Cluster_19_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_19, file = "TMS_merge_T_cells_Cluster_19_VariableGenes.rds")
Cluster_19 <- readRDS(file = "TMS_merge_T_cells_Cluster_19_VariableGenes.rds")

all.genes <- rownames(Cluster_19)
Cluster_19 <- ScaleData(Cluster_19, features = all.genes)

saveRDS(Cluster_19, file = "TMS_merge_T_cells_Cluster_19_ScaledData.rds")
Cluster_19 <- readRDS(file = "TMS_merge_T_cells_Cluster_19_ScaledData.rds")

Cluster_19 <- RunPCA(Cluster_19, npcs = 5, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(Cluster_19, reduction = "pca")

saveRDS(Cluster_19, file = "TMS_merge_T_cells_Cluster_19_PCA.rds")
Cluster_19 <- readRDS(file = "TMS_merge_T_cells_Cluster_19_PCA.rds")

Cluster_19 <- JackStraw(Cluster_19, num.replicate = 100, dims = 5)
Cluster_19 <- ScoreJackStraw(Cluster_19, dims = 1:5)
JackStrawPlot(Cluster_19, dims = 1:5)

ElbowPlot(Cluster_19, ndims = 5)

Cluster_19 <- RunTSNE(Cluster_19, dims = 1:5, perplexity = 10)
DimPlot(Cluster_19, reduction = "tsne")

Cluster_19 <- RunUMAP(Cluster_19, dims = 1:5)
DimPlot(Cluster_19, reduction = "umap")

saveRDS(Cluster_19, file = "TMS_merge_T_cells_Cluster_19_TSNE_UMAP.rds")
Cluster_19 <- readRDS(file = "TMS_merge_T_cells_Cluster_19_TSNE_UMAP.rds")

Cluster_19 <- FindNeighbors(Cluster_19, dims = 1:5)

#louvian
Cluster_19 <- FindClusters(Cluster_19, resolution = seq(0,3, by = 0.1))

clustree(Cluster_19, prefix = "RNA_snn_res.")
clustree(Cluster_19, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_19) <- Cluster_19$RNA_snn_res.3

Cluster_19 <- BuildClusterTree(Cluster_19, reorder.numeric = TRUE, reorder = TRUE, dims = 1:5)
PlotClusterTree(object = Cluster_19)

Cluster_19.markers <- FindAllMarkers(Cluster_19)
m <- Cluster_19.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_19, features = m$gene)
write.csv(Cluster_19.markers, "TMS_T_cells_Cluster_19_merkers.csv")

Cluster_19 <- RenameIdents(Cluster_19, "1" = "DN T", "2" = "CD4 T", "3" = "DN T", "4" = "DN T", "5" = "DN T", 
                           "6" = "DN T", "7" = "DN T", "8" = "DN T", "9" = "DN T", "10" = "DN T", "11" = "DN T", 
                           "12" = "CD4 T", "13" = "DN T", "14" = "DN T", "15" = "DN T", "16" = "CD8 T")

saveRDS(Cluster_19, file = "TMS_merge_T_cells_Cluster_19_clustered.rds")
Cluster_19 <- readRDS(file = "TMS_merge_T_cells_Cluster_19_clustered.rds")

