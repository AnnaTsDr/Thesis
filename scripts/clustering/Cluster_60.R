Cluster_60 <- readRDS(file = "TMS_merge_T_cells_clusters_60.rds")

#Highly variable expressed genes
Cluster_60 <- FindVariableFeatures(Cluster_60, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_60),15)
plot1 <- VariableFeaturePlot(Cluster_60)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_60), file = "TMS_merge_T_cells_Cluster_60_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_60, file = "TMS_merge_T_cells_Cluster_60_VariableGenes.rds")
Cluster_60 <- readRDS(file = "TMS_merge_T_cells_Cluster_60_VariableGenes.rds")

all.genes <- rownames(Cluster_60)
Cluster_60 <- ScaleData(Cluster_60, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_60, file = "TMS_merge_T_cells_Cluster_60_ScaledData.rds")
Cluster_60 <- readRDS(file = "TMS_merge_T_cells_Cluster_60_ScaledData.rds")

Cluster_60 <- RunPCA(Cluster_60, ndims.print = 1:15, nfeatures.print = 5)
DimPlot(Cluster_60, reduction = "pca")

saveRDS(Cluster_60, file = "TMS_merge_T_cells_Cluster_60_PCA.rds")
Cluster_60 <- readRDS(file = "TMS_merge_T_cells_Cluster_60_PCA.rds")

Cluster_60 <- JackStraw(Cluster_60, num.replicate = 100, dims = 15)
Cluster_60 <- ScoreJackStraw(Cluster_60, dims = 1:15)
JackStrawPlot(Cluster_60, dims = 1:15)

ElbowPlot(Cluster_60, ndims = 15)

Cluster_60 <- RunTSNE(Cluster_60, dims = 1:15, perplexity = 10)
DimPlot(Cluster_60, reduction = "tsne")

Cluster_60 <- RunUMAP(Cluster_60, dims = 1:15)
DimPlot(Cluster_60, reduction = "umap")

saveRDS(Cluster_60, file = "TMS_merge_T_cells_Cluster_60_TSNE_UMAP.rds")
Cluster_60 <- readRDS(file = "TMS_merge_T_cells_Cluster_60_TSNE_UMAP.rds")

Cluster_60 <- FindNeighbors(Cluster_60, dims = 1:15)

#louvian
Cluster_60 <- FindClusters(Cluster_60, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_60, prefix = "RNA_snn_res.")
clustree(Cluster_60, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_60) <- Cluster_60$RNA_snn_res.0

#Cluster_60 <- BuildClusterTree(Cluster_60, reorder.numeric = TRUE, reorder = TRUE, dims = 1:15)
#PlotClusterTree(object = Cluster_60)

#Cluster_60.markers <- FindAllMarkers(Cluster_60)
#m <- Cluster_60.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_60, features = m$gene)
#write.csv(Cluster_60.markers, "TMS_T_cells_Cluster_60_merkers.csv")

Cluster_60 <- RenameIdents(Cluster_60, "0" = "CD4 T")

saveRDS(Cluster_60, file = "TMS_merge_T_cells_Cluster_60_clustered.rds")
Cluster_60 <- readRDS(file = "TMS_merge_T_cells_Cluster_60_clustered.rds")

