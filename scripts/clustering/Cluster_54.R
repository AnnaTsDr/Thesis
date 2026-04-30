Cluster_54 <- readRDS(file = "TMS_merge_T_cells_clusters_54.rds")

#Highly variable expressed genes
Cluster_54 <- FindVariableFeatures(Cluster_54, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.8, Inf), mean.cutoff = c(0.1, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_54),15)
plot1 <- VariableFeaturePlot(Cluster_54)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_54), file = "TMS_merge_T_cells_Cluster_54_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_54, file = "TMS_merge_T_cells_Cluster_54_VariableGenes.rds")
Cluster_54 <- readRDS(file = "TMS_merge_T_cells_Cluster_54_VariableGenes.rds")

all.genes <- rownames(Cluster_54)
Cluster_54 <- ScaleData(Cluster_54, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_54, file = "TMS_merge_T_cells_Cluster_54_ScaledData.rds")
Cluster_54 <- readRDS(file = "TMS_merge_T_cells_Cluster_54_ScaledData.rds")

Cluster_54 <- RunPCA(Cluster_54, npcs = 7, ndims.print = 1:7, nfeatures.print = 5)
DimPlot(Cluster_54, reduction = "pca")

saveRDS(Cluster_54, file = "TMS_merge_T_cells_Cluster_54_PCA.rds")
Cluster_54 <- readRDS(file = "TMS_merge_T_cells_Cluster_54_PCA.rds")

Cluster_54 <- JackStraw(Cluster_54, num.replicate = 100, dims = 7)
Cluster_54 <- ScoreJackStraw(Cluster_54, dims = 1:7)
JackStrawPlot(Cluster_54, dims = 1:7)

ElbowPlot(Cluster_54, ndims = 7)

Cluster_54 <- RunTSNE(Cluster_54, dims = 1:7, perplexity = 10)
DimPlot(Cluster_54, reduction = "tsne")

Cluster_54 <- RunUMAP(Cluster_54, dims = 1:7)
DimPlot(Cluster_54, reduction = "umap")

saveRDS(Cluster_54, file = "TMS_merge_T_cells_Cluster_54_TSNE_UMAP.rds")
Cluster_54 <- readRDS(file = "TMS_merge_T_cells_Cluster_54_TSNE_UMAP.rds")

Cluster_54 <- FindNeighbors(Cluster_54, dims = 1:7)

#louvian
Cluster_54 <- FindClusters(Cluster_54, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_54, prefix = "RNA_snn_res.")
clustree(Cluster_54, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_54) <- Cluster_54$RNA_snn_res.1.2

Cluster_54 <- BuildClusterTree(Cluster_54, reorder.numeric = TRUE, reorder = TRUE, dims = 1:7)
PlotClusterTree(object = Cluster_54)

Cluster_54.markers <- FindAllMarkers(Cluster_54)
m <- Cluster_54.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_54, features = m$gene)
write.csv(Cluster_54.markers, "TMS_T_cells_Cluster_54_merkers.csv")

Cluster_54 <- RenameIdents(Cluster_54, "1" = "Not T", "2" = "CD4 T", "3" = "CD4 T", "4" = "CD4 T", "5" = "DN T")

saveRDS(Cluster_54, file = "TMS_merge_T_cells_Cluster_54_clustered.rds")
Cluster_54 <- readRDS(file = "TMS_merge_T_cells_Cluster_54_clustered.rds")

