Cluster_59 <- readRDS(file = "TMS_merge_T_cells_clusters_59.rds")

#Highly variable expressed genes
Cluster_59 <- FindVariableFeatures(Cluster_59, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.1, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_59),15)
plot1 <- VariableFeaturePlot(Cluster_59)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_59), file = "TMS_merge_T_cells_Cluster_59_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_59, file = "TMS_merge_T_cells_Cluster_59_VariableGenes.rds")
Cluster_59 <- readRDS(file = "TMS_merge_T_cells_Cluster_59_VariableGenes.rds")

all.genes <- rownames(Cluster_59)
Cluster_59 <- ScaleData(Cluster_59, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_59, file = "TMS_merge_T_cells_Cluster_59_ScaledData.rds")
Cluster_59 <- readRDS(file = "TMS_merge_T_cells_Cluster_59_ScaledData.rds")

Cluster_59 <- RunPCA(Cluster_59, ndims.print = 1:8, nfeatures.print = 5)
DimPlot(Cluster_59, reduction = "pca")

saveRDS(Cluster_59, file = "TMS_merge_T_cells_Cluster_59_PCA.rds")
Cluster_59 <- readRDS(file = "TMS_merge_T_cells_Cluster_59_PCA.rds")

Cluster_59 <- JackStraw(Cluster_59, num.replicate = 100, dims = 8)
Cluster_59 <- ScoreJackStraw(Cluster_59, dims = 1:8)
JackStrawPlot(Cluster_59, dims = 1:8)

ElbowPlot(Cluster_59, ndims = 8)

Cluster_59 <- RunTSNE(Cluster_59, dims = 1:8, perplexity = 10)
DimPlot(Cluster_59, reduction = "tsne")

Cluster_59 <- RunUMAP(Cluster_59, dims = 1:8)
DimPlot(Cluster_59, reduction = "umap")

saveRDS(Cluster_59, file = "TMS_merge_T_cells_Cluster_59_TSNE_UMAP.rds")
Cluster_59 <- readRDS(file = "TMS_merge_T_cells_Cluster_59_TSNE_UMAP.rds")

Cluster_59 <- FindNeighbors(Cluster_59, dims = 1:8)

#louvian
Cluster_59 <- FindClusters(Cluster_59, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_59, prefix = "RNA_snn_res.")
clustree(Cluster_59, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_59) <- Cluster_59$RNA_snn_res.0.8

Cluster_59 <- BuildClusterTree(Cluster_59, reorder.numeric = TRUE, reorder = TRUE, dims = 1:8)
PlotClusterTree(object = Cluster_59)

Cluster_59.markers <- FindAllMarkers(Cluster_59)
m <- Cluster_59.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_59, features = m$gene)
write.csv(Cluster_59.markers, "TMS_T_cells_Cluster_59_merkers.csv")

Cluster_59 <- RenameIdents(Cluster_59, "1" = "CD4 T", "2" = "CD4 T", "3" = "CD4 T", "4" = "Not T")

saveRDS(Cluster_59, file = "TMS_merge_T_cells_Cluster_59_clustered.rds")
Cluster_59 <- readRDS(file = "TMS_merge_T_cells_Cluster_59_clustered.rds")

