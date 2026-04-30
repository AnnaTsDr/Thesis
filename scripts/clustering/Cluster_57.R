Cluster_57 <- readRDS(file = "TMS_merge_T_cells_clusters_57.rds")

#Highly variable expressed genes
Cluster_57 <- FindVariableFeatures(Cluster_57, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_57),15)
plot1 <- VariableFeaturePlot(Cluster_57)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_57), file = "TMS_merge_T_cells_Cluster_57_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_57, file = "TMS_merge_T_cells_Cluster_57_VariableGenes.rds")
Cluster_57 <- readRDS(file = "TMS_merge_T_cells_Cluster_57_VariableGenes.rds")

all.genes <- rownames(Cluster_57)
Cluster_57 <- ScaleData(Cluster_57, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_57, file = "TMS_merge_T_cells_Cluster_57_ScaledData.rds")
Cluster_57 <- readRDS(file = "TMS_merge_T_cells_Cluster_57_ScaledData.rds")

Cluster_57 <- RunPCA(Cluster_57, ndims.print = 1:6, nfeatures.print = 5)
DimPlot(Cluster_57, reduction = "pca")

saveRDS(Cluster_57, file = "TMS_merge_T_cells_Cluster_57_PCA.rds")
Cluster_57 <- readRDS(file = "TMS_merge_T_cells_Cluster_57_PCA.rds")

Cluster_57 <- JackStraw(Cluster_57, num.replicate = 100, dims = 6)
Cluster_57 <- ScoreJackStraw(Cluster_57, dims = 1:6)
JackStrawPlot(Cluster_57, dims = 1:6)

ElbowPlot(Cluster_57, ndims = 6)

Cluster_57 <- RunTSNE(Cluster_57, dims = 1:6, perplexity = 10)
DimPlot(Cluster_57, reduction = "tsne")

Cluster_57 <- RunUMAP(Cluster_57, dims = 1:6)
DimPlot(Cluster_57, reduction = "umap")

saveRDS(Cluster_57, file = "TMS_merge_T_cells_Cluster_57_TSNE_UMAP.rds")
Cluster_57 <- readRDS(file = "TMS_merge_T_cells_Cluster_57_TSNE_UMAP.rds")

Cluster_57 <- FindNeighbors(Cluster_57, dims = 1:6)

#louvian
Cluster_57 <- FindClusters(Cluster_57, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_57, prefix = "RNA_snn_res.")
clustree(Cluster_57, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_57) <- Cluster_57$RNA_snn_res.1.2

Cluster_57 <- BuildClusterTree(Cluster_57, reorder.numeric = TRUE, reorder = TRUE, dims = 1:6)
PlotClusterTree(object = Cluster_57)

Cluster_57.markers <- FindAllMarkers(Cluster_57)
m <- Cluster_57.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_57, features = m$gene)
write.csv(Cluster_57.markers, "TMS_T_cells_Cluster_57_merkers.csv")

Cluster_57 <- RenameIdents(Cluster_57, "1" = "CD4 T", "2" = "Not T", "3" = "CD4 T", "4" = "CD4 T", "5" = "CD4 T", 
                           "6" = "CD4 T")

saveRDS(Cluster_57, file = "TMS_merge_T_cells_Cluster_57_clustered.rds")
Cluster_57 <- readRDS(file = "TMS_merge_T_cells_Cluster_57_clustered.rds")

