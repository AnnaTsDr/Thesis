Cluster_50 <- readRDS(file = "TMS_merge_T_cells_clusters_50.rds")

#Highly variable expressed genes
Cluster_50 <- FindVariableFeatures(Cluster_50, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.1, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_50),15)
plot1 <- VariableFeaturePlot(Cluster_50)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_50), file = "TMS_merge_T_cells_Cluster_50_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_50, file = "TMS_merge_T_cells_Cluster_50_VariableGenes.rds")
Cluster_50 <- readRDS(file = "TMS_merge_T_cells_Cluster_50_VariableGenes.rds")

all.genes <- rownames(Cluster_50)
Cluster_50 <- ScaleData(Cluster_50, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_50, file = "TMS_merge_T_cells_Cluster_50_ScaledData.rds")
Cluster_50 <- readRDS(file = "TMS_merge_T_cells_Cluster_50_ScaledData.rds")

Cluster_50 <- RunPCA(Cluster_50, ndims.print = 1:2, nfeatures.print = 5)
DimPlot(Cluster_50, reduction = "pca")

saveRDS(Cluster_50, file = "TMS_merge_T_cells_Cluster_50_PCA.rds")
Cluster_50 <- readRDS(file = "TMS_merge_T_cells_Cluster_50_PCA.rds")

Cluster_50 <- JackStraw(Cluster_50, num.replicate = 100, dims = 2)
Cluster_50 <- ScoreJackStraw(Cluster_50, dims = 1:2)
JackStrawPlot(Cluster_50, dims = 1:2)

ElbowPlot(Cluster_50, ndims = 2)

Cluster_50 <- RunTSNE(Cluster_50, dims = 1:2, perplexity = 10)
DimPlot(Cluster_50, reduction = "tsne")

Cluster_50 <- RunUMAP(Cluster_50, dims = 1:2)
DimPlot(Cluster_50, reduction = "umap")

saveRDS(Cluster_50, file = "TMS_merge_T_cells_Cluster_50_TSNE_UMAP.rds")
Cluster_50 <- readRDS(file = "TMS_merge_T_cells_Cluster_50_TSNE_UMAP.rds")

Cluster_50 <- FindNeighbors(Cluster_50, dims = 1:2)

#louvian
Cluster_50 <- FindClusters(Cluster_50, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_50, prefix = "RNA_snn_res.")
clustree(Cluster_50, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_50) <- Cluster_50$RNA_snn_res.0

#Cluster_50 <- BuildClusterTree(Cluster_50, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
#PlotClusterTree(object = Cluster_50)

#Cluster_50.markers <- FindAllMarkers(Cluster_50)
#m <- Cluster_50.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_50, features = m$gene)
#write.csv(Cluster_50.markers, "TMS_T_cells_Cluster_50_merkers.csv")

Cluster_50 <- RenameIdents(Cluster_50, "0" = "DN T")

saveRDS(Cluster_50, file = "TMS_merge_T_cells_Cluster_50_clustered.rds")
Cluster_50 <- readRDS(file = "TMS_merge_T_cells_Cluster_50_clustered.rds")

