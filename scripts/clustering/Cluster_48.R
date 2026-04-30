Cluster_48 <- readRDS(file = "TMS_merge_T_cells_clusters_48.rds")

#Highly variable expressed genes
Cluster_48 <- FindVariableFeatures(Cluster_48, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.3, Inf), mean.cutoff = c(0.05, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_48),15)
plot1 <- VariableFeaturePlot(Cluster_48)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_48), file = "TMS_merge_T_cells_Cluster_48_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_48, file = "TMS_merge_T_cells_Cluster_48_VariableGenes.rds")
Cluster_48 <- readRDS(file = "TMS_merge_T_cells_Cluster_48_VariableGenes.rds")

all.genes <- rownames(Cluster_48)
Cluster_48 <- ScaleData(Cluster_48, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_48, file = "TMS_merge_T_cells_Cluster_48_ScaledData.rds")
Cluster_48 <- readRDS(file = "TMS_merge_T_cells_Cluster_48_ScaledData.rds")

Cluster_48 <- RunPCA(Cluster_48, ndims.print = 1:10, nfeatures.print = 5)
DimPlot(Cluster_48, reduction = "pca")

saveRDS(Cluster_48, file = "TMS_merge_T_cells_Cluster_48_PCA.rds")
Cluster_48 <- readRDS(file = "TMS_merge_T_cells_Cluster_48_PCA.rds")

Cluster_48 <- JackStraw(Cluster_48, num.replicate = 100, dims = 10)
Cluster_48 <- ScoreJackStraw(Cluster_48, dims = 1:10)
JackStrawPlot(Cluster_48, dims = 1:10)

ElbowPlot(Cluster_48, ndims = 10)

Cluster_48 <- RunTSNE(Cluster_48, dims = 1:10, perplexity = 10)
DimPlot(Cluster_48, reduction = "tsne")

Cluster_48 <- RunUMAP(Cluster_48, dims = 1:10)
DimPlot(Cluster_48, reduction = "umap")

saveRDS(Cluster_48, file = "TMS_merge_T_cells_Cluster_48_TSNE_UMAP.rds")
Cluster_48 <- readRDS(file = "TMS_merge_T_cells_Cluster_48_TSNE_UMAP.rds")

Cluster_48 <- FindNeighbors(Cluster_48, dims = 1:10)

#louvian
Cluster_48 <- FindClusters(Cluster_48, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_48, prefix = "RNA_snn_res.")
clustree(Cluster_48, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_48) <- Cluster_48$RNA_snn_res.0

#Cluster_48 <- BuildClusterTree(Cluster_48, reorder.numeric = TRUE, reorder = TRUE, dims = 1:5)
#PlotClusterTree(object = Cluster_48)

#Cluster_48.markers <- FindAllMarkers(Cluster_48)
#m <- Cluster_48.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_48, features = m$gene)
#write.csv(Cluster_48.markers, "TMS_T_cells_Cluster_48_merkers.csv")

Cluster_48 <- RenameIdents(Cluster_48, "0" = "DN T")

saveRDS(Cluster_48, file = "TMS_merge_T_cells_Cluster_48_clustered.rds")
Cluster_48 <- readRDS(file = "TMS_merge_T_cells_Cluster_48_clustered.rds")

