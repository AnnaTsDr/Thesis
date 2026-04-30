Cluster_31 <- readRDS(file = "TMS_merge_T_cells_clusters_31.rds")

#Highly variable expressed genes
Cluster_31 <- FindVariableFeatures(Cluster_31, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_31),15)
plot1 <- VariableFeaturePlot(Cluster_31)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_31), file = "TMS_merge_T_cells_Cluster_31_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_31, file = "TMS_merge_T_cells_Cluster_31_VariableGenes.rds")
Cluster_31 <- readRDS(file = "TMS_merge_T_cells_Cluster_31_VariableGenes.rds")

all.genes <- rownames(Cluster_31)
Cluster_31 <- ScaleData(Cluster_31, features = all.genes)

saveRDS(Cluster_31, file = "TMS_merge_T_cells_Cluster_31_ScaledData.rds")
Cluster_31 <- readRDS(file = "TMS_merge_T_cells_Cluster_31_ScaledData.rds")

Cluster_31 <- RunPCA(Cluster_31, npcs = 8, ndims.print = 1:8, nfeatures.print = 5)
DimPlot(Cluster_31, reduction = "pca")

saveRDS(Cluster_31, file = "TMS_merge_T_cells_Cluster_31_PCA.rds")
Cluster_31 <- readRDS(file = "TMS_merge_T_cells_Cluster_31_PCA.rds")

Cluster_31 <- JackStraw(Cluster_31, num.replicate = 100, dims = 8)
Cluster_31 <- ScoreJackStraw(Cluster_31, dims = 1:8)
JackStrawPlot(Cluster_31, dims = 1:8)

ElbowPlot(Cluster_31, ndims = 8)

Cluster_31 <- RunTSNE(Cluster_31, dims = 1:8, perplexity = 10)
DimPlot(Cluster_31, reduction = "tsne")

Cluster_31 <- RunUMAP(Cluster_31, dims = 1:8)
DimPlot(Cluster_31, reduction = "umap")

saveRDS(Cluster_31, file = "TMS_merge_T_cells_Cluster_31_TSNE_UMAP.rds")
Cluster_31 <- readRDS(file = "TMS_merge_T_cells_Cluster_31_TSNE_UMAP.rds")

Cluster_31 <- FindNeighbors(Cluster_31, dims = 1:8)

#louvian
Cluster_31 <- FindClusters(Cluster_31, resolution = seq(0,3, by = 0.1))

clustree(Cluster_31, prefix = "RNA_snn_res.")
clustree(Cluster_31, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_31) <- Cluster_31$RNA_snn_res.3

Cluster_31 <- BuildClusterTree(Cluster_31, reorder.numeric = TRUE, reorder = TRUE, dims = 1:8)
PlotClusterTree(object = Cluster_31)

Cluster_31.markers <- FindAllMarkers(Cluster_31)
m <- Cluster_31.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_31, features = m$gene)
write.csv(Cluster_31.markers, "TMS_T_cells_Cluster_31_merkers.csv")

Cluster_31 <- RenameIdents(Cluster_31, "1" = "CD4 T", "2" = "CD4 T", "3" = "CD4 T", "4" = "CD4 T", "5" = "CD4 T", 
                           "6" = "CD4 T", "7" = "CD4 T", "8" = "CD4 T", "9" = "CD4 T", "10" = "CD4 T", "11" = "CD4 T", 
                           "12" = "DN T", "13" = "CD4 T", "14" = "CD4 T", "15" = "CD4 T", "16" = "CD4 T")

saveRDS(Cluster_31, file = "TMS_merge_T_cells_Cluster_31_clustered.rds")
Cluster_31 <- readRDS(file = "TMS_merge_T_cells_Cluster_31_clustered.rds")
