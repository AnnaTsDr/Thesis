Cluster_68 <- readRDS(file = "TMS_merge_T_cells_clusters_68.rds")

#Highly variable expressed genes
Cluster_68 <- FindVariableFeatures(Cluster_68, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_68),15)
plot1 <- VariableFeaturePlot(Cluster_68)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_68), file = "TMS_merge_T_cells_Cluster_68_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_68, file = "TMS_merge_T_cells_Cluster_68_VariableGenes.rds")
Cluster_68 <- readRDS(file = "TMS_merge_T_cells_Cluster_68_VariableGenes.rds")

all.genes <- rownames(Cluster_68)
Cluster_68 <- ScaleData(Cluster_68, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_68, file = "TMS_merge_T_cells_Cluster_68_ScaledData.rds")
Cluster_68 <- readRDS(file = "TMS_merge_T_cells_Cluster_68_ScaledData.rds")

Cluster_68 <- RunPCA(Cluster_68, ndims.print = 1:21, nfeatures.print = 5)
DimPlot(Cluster_68, reduction = "pca")

saveRDS(Cluster_68, file = "TMS_merge_T_cells_Cluster_68_PCA.rds")
Cluster_68 <- readRDS(file = "TMS_merge_T_cells_Cluster_68_PCA.rds")

Cluster_68 <- JackStraw(Cluster_68, num.replicate = 100, dims = 21)
Cluster_68 <- ScoreJackStraw(Cluster_68, dims = 1:21)
JackStrawPlot(Cluster_68, dims = 1:21)

ElbowPlot(Cluster_68, ndims = 21)

Cluster_68 <- RunTSNE(Cluster_68, dims = 1:21, perplexity = 10)
DimPlot(Cluster_68, reduction = "tsne")

Cluster_68 <- RunUMAP(Cluster_68, dims = 1:21)
DimPlot(Cluster_68, reduction = "umap")

saveRDS(Cluster_68, file = "TMS_merge_T_cells_Cluster_68_TSNE_UMAP.rds")
Cluster_68 <- readRDS(file = "TMS_merge_T_cells_Cluster_68_TSNE_UMAP.rds")

Cluster_68 <- FindNeighbors(Cluster_68, dims = 1:21)

#louvian
Cluster_68 <- FindClusters(Cluster_68, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_68, prefix = "RNA_snn_res.")
clustree(Cluster_68, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_68) <- Cluster_68$RNA_snn_res.0

#Cluster_68 <- BuildClusterTree(Cluster_68, reorder.numeric = TRUE, reorder = TRUE, dims = 1:15)
#PlotClusterTree(object = Cluster_68)

#Cluster_68.markers <- FindAllMarkers(Cluster_68)
#m <- Cluster_68.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_68, features = m$gene)
#write.csv(Cluster_68.markers, "TMS_T_cells_Cluster_68_merkers.csv")

Cluster_68 <- RenameIdents(Cluster_68, "0" = "DP T")

saveRDS(Cluster_68, file = "TMS_merge_T_cells_Cluster_68_clustered.rds")
Cluster_68 <- readRDS(file = "TMS_merge_T_cells_Cluster_68_clustered.rds")

