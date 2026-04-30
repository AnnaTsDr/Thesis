Cluster_51 <- readRDS(file = "TMS_merge_T_cells_clusters_51.rds")

#Highly variable expressed genes
Cluster_51 <- FindVariableFeatures(Cluster_51, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.3, Inf), mean.cutoff = c(0.1, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_51),15)
plot1 <- VariableFeaturePlot(Cluster_51)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_51), file = "TMS_merge_T_cells_Cluster_51_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_51, file = "TMS_merge_T_cells_Cluster_51_VariableGenes.rds")
Cluster_51 <- readRDS(file = "TMS_merge_T_cells_Cluster_51_VariableGenes.rds")

all.genes <- rownames(Cluster_51)
Cluster_51 <- ScaleData(Cluster_51, features = all.genes, vars.to.regress = c("age", "mouse.id"))

saveRDS(Cluster_51, file = "TMS_merge_T_cells_Cluster_51_ScaledData.rds")
Cluster_51 <- readRDS(file = "TMS_merge_T_cells_Cluster_51_ScaledData.rds")

Cluster_51 <- RunPCA(Cluster_51, ndims.print = 1:6, nfeatures.print = 5)
DimPlot(Cluster_51, reduction = "pca")

saveRDS(Cluster_51, file = "TMS_merge_T_cells_Cluster_51_PCA.rds")
Cluster_51 <- readRDS(file = "TMS_merge_T_cells_Cluster_51_PCA.rds")

Cluster_51 <- JackStraw(Cluster_51, num.replicate = 100, dims = 6)
Cluster_51 <- ScoreJackStraw(Cluster_51, dims = 1:6)
JackStrawPlot(Cluster_51, dims = 1:6)

ElbowPlot(Cluster_51, ndims = 6)

Cluster_51 <- RunTSNE(Cluster_51, dims = 1:6, perplexity = 10)
DimPlot(Cluster_51, reduction = "tsne")

Cluster_51 <- RunUMAP(Cluster_51, dims = 1:6)
DimPlot(Cluster_51, reduction = "umap")

saveRDS(Cluster_51, file = "TMS_merge_T_cells_Cluster_51_TSNE_UMAP.rds")
Cluster_51 <- readRDS(file = "TMS_merge_T_cells_Cluster_51_TSNE_UMAP.rds")

Cluster_51 <- FindNeighbors(Cluster_51, dims = 1:6)

#louvian
Cluster_51 <- FindClusters(Cluster_51, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_51, prefix = "RNA_snn_res.")
clustree(Cluster_51, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_51) <- Cluster_51$RNA_snn_res.0

#Cluster_51 <- BuildClusterTree(Cluster_51, reorder.numeric = TRUE, reorder = TRUE, dims = 1:6)
#PlotClusterTree(object = Cluster_51)

#Cluster_51.markers <- FindAllMarkers(Cluster_51)
#m <- Cluster_51.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_51, features = m$gene)
#write.csv(Cluster_51.markers, "TMS_T_cells_Cluster_51_merkers.csv")

Cluster_51 <- RenameIdents(Cluster_51, "0" = "CD8 T")

saveRDS(Cluster_51, file = "TMS_merge_T_cells_Cluster_51_clustered.rds")
Cluster_51 <- readRDS(file = "TMS_merge_T_cells_Cluster_51_clustered.rds")

