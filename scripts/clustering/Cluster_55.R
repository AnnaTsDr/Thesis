Cluster_55 <- readRDS(file = "TMS_merge_T_cells_clusters_55.rds")

#Highly variable expressed genes
Cluster_55 <- FindVariableFeatures(Cluster_55, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.1, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_55),15)
plot1 <- VariableFeaturePlot(Cluster_55)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_55), file = "TMS_merge_T_cells_Cluster_55_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_55, file = "TMS_merge_T_cells_Cluster_55_VariableGenes.rds")
Cluster_55 <- readRDS(file = "TMS_merge_T_cells_Cluster_55_VariableGenes.rds")

all.genes <- rownames(Cluster_55)
Cluster_55 <- ScaleData(Cluster_55, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_55, file = "TMS_merge_T_cells_Cluster_55_ScaledData.rds")
Cluster_55 <- readRDS(file = "TMS_merge_T_cells_Cluster_55_ScaledData.rds")

Cluster_55 <- RunPCA(Cluster_55, ndims.print = 1:3, nfeatures.print = 5)
DimPlot(Cluster_55, reduction = "pca")

saveRDS(Cluster_55, file = "TMS_merge_T_cells_Cluster_55_PCA.rds")
Cluster_55 <- readRDS(file = "TMS_merge_T_cells_Cluster_55_PCA.rds")

Cluster_55 <- JackStraw(Cluster_55, num.replicate = 100, dims = 3)
Cluster_55 <- ScoreJackStraw(Cluster_55, dims = 1:3)
JackStrawPlot(Cluster_55, dims = 1:3)

ElbowPlot(Cluster_55, ndims = 3)

Cluster_55 <- RunTSNE(Cluster_55, dims = 1:3, perplexity = 10)
DimPlot(Cluster_55, reduction = "tsne")

Cluster_55 <- RunUMAP(Cluster_55, dims = 1:3)
DimPlot(Cluster_55, reduction = "umap")

saveRDS(Cluster_55, file = "TMS_merge_T_cells_Cluster_55_TSNE_UMAP.rds")
Cluster_55 <- readRDS(file = "TMS_merge_T_cells_Cluster_55_TSNE_UMAP.rds")

Cluster_55 <- FindNeighbors(Cluster_55, dims = 1:3)

#louvian
Cluster_55 <- FindClusters(Cluster_55, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_55, prefix = "RNA_snn_res.")
clustree(Cluster_55, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_55) <- Cluster_55$RNA_snn_res.0

#Cluster_55 <- BuildClusterTree(Cluster_55, reorder.numeric = TRUE, reorder = TRUE, dims = 1:15)
#PlotClusterTree(object = Cluster_55)

#Cluster_55.markers <- FindAllMarkers(Cluster_55)
#m <- Cluster_55.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_55, features = m$gene)
#write.csv(Cluster_55.markers, "TMS_T_cells_Cluster_55_merkers.csv")

Cluster_55 <- RenameIdents(Cluster_55, "0" = "CD4 T")

saveRDS(Cluster_55, file = "TMS_merge_T_cells_Cluster_55_clustered.rds")
Cluster_55 <- readRDS(file = "TMS_merge_T_cells_Cluster_55_clustered.rds")

