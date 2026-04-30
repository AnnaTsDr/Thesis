Cluster_58 <- readRDS(file = "TMS_merge_T_cells_clusters_58.rds")

#Highly variable expressed genes
Cluster_58 <- FindVariableFeatures(Cluster_58, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_58),15)
plot1 <- VariableFeaturePlot(Cluster_58)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_58), file = "TMS_merge_T_cells_Cluster_58_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_58, file = "TMS_merge_T_cells_Cluster_58_VariableGenes.rds")
Cluster_58 <- readRDS(file = "TMS_merge_T_cells_Cluster_58_VariableGenes.rds")

all.genes <- rownames(Cluster_58)
Cluster_58 <- ScaleData(Cluster_58, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_58, file = "TMS_merge_T_cells_Cluster_58_ScaledData.rds")
Cluster_58 <- readRDS(file = "TMS_merge_T_cells_Cluster_58_ScaledData.rds")

Cluster_58 <- RunPCA(Cluster_58, ndims.print = 1:3, nfeatures.print = 5)
DimPlot(Cluster_58, reduction = "pca")

saveRDS(Cluster_58, file = "TMS_merge_T_cells_Cluster_58_PCA.rds")
Cluster_58 <- readRDS(file = "TMS_merge_T_cells_Cluster_58_PCA.rds")

Cluster_58 <- JackStraw(Cluster_58, num.replicate = 100, dims = 3)
Cluster_58 <- ScoreJackStraw(Cluster_58, dims = 1:3)
JackStrawPlot(Cluster_58, dims = 1:3)

ElbowPlot(Cluster_58, ndims = 3)

Cluster_58 <- RunTSNE(Cluster_58, dims = 1:3, perplexity = 10)
DimPlot(Cluster_58, reduction = "tsne")

Cluster_58 <- RunUMAP(Cluster_58, dims = 1:3)
DimPlot(Cluster_58, reduction = "umap")

saveRDS(Cluster_58, file = "TMS_merge_T_cells_Cluster_58_TSNE_UMAP.rds")
Cluster_58 <- readRDS(file = "TMS_merge_T_cells_Cluster_58_TSNE_UMAP.rds")

Cluster_58 <- FindNeighbors(Cluster_58, dims = 1:3)

#louvian
Cluster_58 <- FindClusters(Cluster_58, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_58, prefix = "RNA_snn_res.")
clustree(Cluster_58, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_58) <- Cluster_58$RNA_snn_res.0

#Cluster_58 <- BuildClusterTree(Cluster_58, reorder.numeric = TRUE, reorder = TRUE, dims = 1:2)
#PlotClusterTree(object = Cluster_58)

#Cluster_58.markers <- FindAllMarkers(Cluster_58)
#m <- Cluster_58.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_58, features = m$gene)
#write.csv(Cluster_58.markers, "TMS_T_cells_Cluster_58_merkers.csv")

Cluster_58 <- RenameIdents(Cluster_58, "0" = "CD4 T")

saveRDS(Cluster_58, file = "TMS_merge_T_cells_Cluster_58_clustered.rds")
Cluster_58 <- readRDS(file = "TMS_merge_T_cells_Cluster_58_clustered.rds")

