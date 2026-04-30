Cluster_43 <- readRDS(file = "TMS_merge_T_cells_clusters_43.rds")

#Highly variable expressed genes
Cluster_43 <- FindVariableFeatures(Cluster_43, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.8, Inf), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_43),15)
plot1 <- VariableFeaturePlot(Cluster_43)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_43), file = "TMS_merge_T_cells_Cluster_43_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_43, file = "TMS_merge_T_cells_Cluster_43_VariableGenes.rds")
Cluster_43 <- readRDS(file = "TMS_merge_T_cells_Cluster_43_VariableGenes.rds")

all.genes <- rownames(Cluster_43)
Cluster_43 <- ScaleData(Cluster_43, features = all.genes)

saveRDS(Cluster_43, file = "TMS_merge_T_cells_Cluster_43_ScaledData.rds")
Cluster_43 <- readRDS(file = "TMS_merge_T_cells_Cluster_43_ScaledData.rds")

Cluster_43 <- RunPCA(Cluster_43, npcs = 7, ndims.print = 1:7, nfeatures.print = 5)
DimPlot(Cluster_43, reduction = "pca")

saveRDS(Cluster_43, file = "TMS_merge_T_cells_Cluster_43_PCA.rds")
Cluster_43 <- readRDS(file = "TMS_merge_T_cells_Cluster_43_PCA.rds")

Cluster_43 <- JackStraw(Cluster_43, num.replicate = 100, dims = 7)
Cluster_43 <- ScoreJackStraw(Cluster_43, dims = 1:7)
JackStrawPlot(Cluster_43, dims = 1:7)

ElbowPlot(Cluster_43, ndims = 7)

Cluster_43 <- RunTSNE(Cluster_43, dims = 1:7, perplexity = 10)
DimPlot(Cluster_43, reduction = "tsne")

Cluster_43 <- RunUMAP(Cluster_43, dims = 1:7)
DimPlot(Cluster_43, reduction = "umap")

saveRDS(Cluster_43, file = "TMS_merge_T_cells_Cluster_43_TSNE_UMAP.rds")
Cluster_43 <- readRDS(file = "TMS_merge_T_cells_Cluster_43_TSNE_UMAP.rds")

Cluster_43 <- FindNeighbors(Cluster_43, dims = 1:7)

#louvian
Cluster_43 <- FindClusters(Cluster_43, resolution = seq(0,3, by = 0.1))

clustree(Cluster_43, prefix = "RNA_snn_res.")
clustree(Cluster_43, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_43) <- Cluster_43$RNA_snn_res.3

Cluster_43 <- BuildClusterTree(Cluster_43, reorder.numeric = TRUE, reorder = TRUE, dims = 1:7)
PlotClusterTree(object = Cluster_43)

Cluster_43.markers <- FindAllMarkers(Cluster_43)
m <- Cluster_43.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_43, features = m$gene)
write.csv(Cluster_43.markers, "TMS_T_cells_Cluster_43_merkers.csv")

Cluster_43 <- RenameIdents(Cluster_43, "1" = "CD8 T", "2" = "CD4 T", "3" = "CD8 T", "4" = "CD4 T", "5" = "CD8 T", 
                           "6" = "DN T", "7" = "CD4 T", "8" = "CD4 T", "9" = "CD4 T", "10" = "CD4 T", "11" = "DN T",
                           "12" = "CD8 T", "13" = "CD8 T", "14" = "CD8 T", "15" = "CD8 T", "16" = "CD8 T",
                           "17" = "CD8 T", "18" = "CD8 T", "19" = "CD8 T")

saveRDS(Cluster_43, file = "TMS_merge_T_cells_Cluster_43_clustered.rds")
Cluster_43 <- readRDS(file = "TMS_merge_T_cells_Cluster_43_clustered.rds")

