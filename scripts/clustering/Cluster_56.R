Cluster_56 <- readRDS(file = "TMS_merge_T_cells_clusters_56.rds")

#Highly variable expressed genes
Cluster_56 <- FindVariableFeatures(Cluster_56, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_56),15)
plot1 <- VariableFeaturePlot(Cluster_56)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_56), file = "TMS_merge_T_cells_Cluster_56_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_56, file = "TMS_merge_T_cells_Cluster_56_VariableGenes.rds")
Cluster_56 <- readRDS(file = "TMS_merge_T_cells_Cluster_56_VariableGenes.rds")

all.genes <- rownames(Cluster_56)
Cluster_56 <- ScaleData(Cluster_56, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_56, file = "TMS_merge_T_cells_Cluster_56_ScaledData.rds")
Cluster_56 <- readRDS(file = "TMS_merge_T_cells_Cluster_56_ScaledData.rds")

Cluster_56 <- RunPCA(Cluster_56, ndims.print = 1:3, nfeatures.print = 5)
DimPlot(Cluster_56, reduction = "pca")

saveRDS(Cluster_56, file = "TMS_merge_T_cells_Cluster_56_PCA.rds")
Cluster_56 <- readRDS(file = "TMS_merge_T_cells_Cluster_56_PCA.rds")

Cluster_56 <- JackStraw(Cluster_56, num.replicate = 100, dims = 3)
Cluster_56 <- ScoreJackStraw(Cluster_56, dims = 1:3)
JackStrawPlot(Cluster_56, dims = 1:3)

ElbowPlot(Cluster_56, ndims = 3)

Cluster_56 <- RunTSNE(Cluster_56, dims = 1:3, perplexity = 10)
DimPlot(Cluster_56, reduction = "tsne")

Cluster_56 <- RunUMAP(Cluster_56, dims = 1:3)
DimPlot(Cluster_56, reduction = "umap")

saveRDS(Cluster_56, file = "TMS_merge_T_cells_Cluster_56_TSNE_UMAP.rds")
Cluster_56 <- readRDS(file = "TMS_merge_T_cells_Cluster_56_TSNE_UMAP.rds")

Cluster_56 <- FindNeighbors(Cluster_56, dims = 1:3)

#louvian
Cluster_56 <- FindClusters(Cluster_56, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_56, prefix = "RNA_snn_res.")
clustree(Cluster_56, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_56) <- Cluster_56$RNA_snn_res.1.2

Cluster_56 <- BuildClusterTree(Cluster_56, reorder.numeric = TRUE, reorder = TRUE, dims = 1:3)
PlotClusterTree(object = Cluster_56)

Cluster_56.markers <- FindAllMarkers(Cluster_56)
m <- Cluster_56.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_56, features = m$gene)
write.csv(Cluster_56.markers, "TMS_T_cells_Cluster_56_merkers.csv")

Cluster_56 <- RenameIdents(Cluster_56, "1" = "Not T", "2" = "DN T", "3" = "CD8 T", "4" = "CD8 T", "5" = "CD4 T")

saveRDS(Cluster_56, file = "TMS_merge_T_cells_Cluster_56_clustered.rds")
Cluster_56 <- readRDS(file = "TMS_merge_T_cells_Cluster_56_clustered.rds")

