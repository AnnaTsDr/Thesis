Cluster_63 <- readRDS(file = "TMS_merge_T_cells_clusters_63.rds")

#Highly variable expressed genes
Cluster_63 <- FindVariableFeatures(Cluster_63, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_63),15)
plot1 <- VariableFeaturePlot(Cluster_63)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_63), file = "TMS_merge_T_cells_Cluster_63_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_63, file = "TMS_merge_T_cells_Cluster_63_VariableGenes.rds")
Cluster_63 <- readRDS(file = "TMS_merge_T_cells_Cluster_63_VariableGenes.rds")

all.genes <- rownames(Cluster_63)
Cluster_63 <- ScaleData(Cluster_63, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_63, file = "TMS_merge_T_cells_Cluster_63_ScaledData.rds")
Cluster_63 <- readRDS(file = "TMS_merge_T_cells_Cluster_63_ScaledData.rds")

Cluster_63 <- RunPCA(Cluster_63, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(Cluster_63, reduction = "pca")

saveRDS(Cluster_63, file = "TMS_merge_T_cells_Cluster_63_PCA.rds")
Cluster_63 <- readRDS(file = "TMS_merge_T_cells_Cluster_63_PCA.rds")

Cluster_63 <- JackStraw(Cluster_63, num.replicate = 100, dims = 5)
Cluster_63 <- ScoreJackStraw(Cluster_63, dims = 1:5)
JackStrawPlot(Cluster_63, dims = 1:5)

ElbowPlot(Cluster_63, ndims = 5)

Cluster_63 <- RunTSNE(Cluster_63, dims = 1:5, perplexity = 10)
DimPlot(Cluster_63, reduction = "tsne")

Cluster_63 <- RunUMAP(Cluster_63, dims = 1:5)
DimPlot(Cluster_63, reduction = "umap")

saveRDS(Cluster_63, file = "TMS_merge_T_cells_Cluster_63_TSNE_UMAP.rds")
Cluster_63 <- readRDS(file = "TMS_merge_T_cells_Cluster_63_TSNE_UMAP.rds")

Cluster_63 <- FindNeighbors(Cluster_63, dims = 1:5)

#louvian
Cluster_63 <- FindClusters(Cluster_63, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_63, prefix = "RNA_snn_res.")
clustree(Cluster_63, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_63) <- Cluster_63$RNA_snn_res.1.2

Cluster_63 <- BuildClusterTree(Cluster_63, reorder.numeric = TRUE, reorder = TRUE, dims = 1:5)
PlotClusterTree(object = Cluster_63)

Cluster_63.markers <- FindAllMarkers(Cluster_63)
m <- Cluster_63.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_63, features = m$gene)
write.csv(Cluster_63.markers, "TMS_T_cells_Cluster_63_merkers.csv")

Cluster_63 <- RenameIdents(Cluster_63, "1" = "CD4 T", "2" = "CD8 T", "3" = "CD8 T", "4" = "B cell", "5" = "CD8 T", 
                           "6" = "CD8 T")

saveRDS(Cluster_63, file = "TMS_merge_T_cells_Cluster_63_clustered.rds")
Cluster_63 <- readRDS(file = "TMS_merge_T_cells_Cluster_63_clustered.rds")

