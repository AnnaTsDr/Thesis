Cluster_52 <- readRDS(file = "TMS_merge_T_cells_clusters_52.rds")

#Highly variable expressed genes
Cluster_52 <- FindVariableFeatures(Cluster_52, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_52),15)
plot1 <- VariableFeaturePlot(Cluster_52)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_52), file = "TMS_merge_T_cells_Cluster_52_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_52, file = "TMS_merge_T_cells_Cluster_52_VariableGenes.rds")
Cluster_52 <- readRDS(file = "TMS_merge_T_cells_Cluster_52_VariableGenes.rds")

all.genes <- rownames(Cluster_52)
Cluster_52 <- ScaleData(Cluster_52, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_52, file = "TMS_merge_T_cells_Cluster_52_ScaledData.rds")
Cluster_52 <- readRDS(file = "TMS_merge_T_cells_Cluster_52_ScaledData.rds")

Cluster_52 <- RunPCA(Cluster_52, npcs = 4, ndims.print = 1:4, nfeatures.print = 5)
DimPlot(Cluster_52, reduction = "pca")

saveRDS(Cluster_52, file = "TMS_merge_T_cells_Cluster_52_PCA.rds")
Cluster_52 <- readRDS(file = "TMS_merge_T_cells_Cluster_52_PCA.rds")

Cluster_52 <- JackStraw(Cluster_52, num.replicate = 100, dims = 4)
Cluster_52 <- ScoreJackStraw(Cluster_52, dims = 1:4)
JackStrawPlot(Cluster_52, dims = 1:4)

ElbowPlot(Cluster_52, ndims = 4)

Cluster_52 <- RunTSNE(Cluster_52, dims = 1:4, perplexity = 10)
DimPlot(Cluster_52, reduction = "tsne")

Cluster_52 <- RunUMAP(Cluster_52, dims = 1:4)
DimPlot(Cluster_52, reduction = "umap")

saveRDS(Cluster_52, file = "TMS_merge_T_cells_Cluster_52_TSNE_UMAP.rds")
Cluster_52 <- readRDS(file = "TMS_merge_T_cells_Cluster_52_TSNE_UMAP.rds")

Cluster_52 <- FindNeighbors(Cluster_52, dims = 1:4)

#louvian
Cluster_52 <- FindClusters(Cluster_52, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_52, prefix = "RNA_snn_res.")
clustree(Cluster_52, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_52) <- Cluster_52$RNA_snn_res.0

#Cluster_52 <- BuildClusterTree(Cluster_52, reorder.numeric = TRUE, reorder = TRUE, dims = 1:4)
#PlotClusterTree(object = Cluster_52)

#Cluster_52.markers <- FindAllMarkers(Cluster_52)
#m <- Cluster_52.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_52, features = m$gene)
#write.csv(Cluster_52.markers, "TMS_T_cells_Cluster_52_merkers.csv")

Cluster_52 <- RenameIdents(Cluster_52, "0" = "CD8 T")

saveRDS(Cluster_52, file = "TMS_merge_T_cells_Cluster_52_clustered.rds")
Cluster_52 <- readRDS(file = "TMS_merge_T_cells_Cluster_52_clustered.rds")

