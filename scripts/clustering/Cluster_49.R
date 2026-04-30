Cluster_49 <- readRDS(file = "TMS_merge_T_cells_clusters_49.rds")

#Highly variable expressed genes
Cluster_49 <- FindVariableFeatures(Cluster_49, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.95, Inf), mean.cutoff = c(0.0125, 7))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_49),15)
plot1 <- VariableFeaturePlot(Cluster_49)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_49), file = "TMS_merge_T_cells_Cluster_49_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_49, file = "TMS_merge_T_cells_Cluster_49_VariableGenes.rds")
Cluster_49 <- readRDS(file = "TMS_merge_T_cells_Cluster_49_VariableGenes.rds")

all.genes <- rownames(Cluster_49)
Cluster_49 <- ScaleData(Cluster_49, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_49, file = "TMS_merge_T_cells_Cluster_49_ScaledData.rds")
Cluster_49 <- readRDS(file = "TMS_merge_T_cells_Cluster_49_ScaledData.rds")

Cluster_49 <- RunPCA(Cluster_49, ndims.print = 1:7, nfeatures.print = 5)
DimPlot(Cluster_49, reduction = "pca")

saveRDS(Cluster_49, file = "TMS_merge_T_cells_Cluster_49_PCA.rds")
Cluster_49 <- readRDS(file = "TMS_merge_T_cells_Cluster_49_PCA.rds")

Cluster_49 <- JackStraw(Cluster_49, num.replicate = 100, dims = 7)
Cluster_49 <- ScoreJackStraw(Cluster_49, dims = 1:7)
JackStrawPlot(Cluster_49, dims = 1:7)

ElbowPlot(Cluster_49, ndims = 7)

Cluster_49 <- RunTSNE(Cluster_49, dims = 1:7, perplexity = 10)
DimPlot(Cluster_49, reduction = "tsne")

Cluster_49 <- RunUMAP(Cluster_49, dims = 1:7)
DimPlot(Cluster_49, reduction = "umap")

saveRDS(Cluster_49, file = "TMS_merge_T_cells_Cluster_49_TSNE_UMAP.rds")
Cluster_49 <- readRDS(file = "TMS_merge_T_cells_Cluster_49_TSNE_UMAP.rds")

Cluster_49 <- FindNeighbors(Cluster_49, dims = 1:7)

#louvian
Cluster_49 <- FindClusters(Cluster_49, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_49, prefix = "RNA_snn_res.")
clustree(Cluster_49, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_49) <- Cluster_49$RNA_snn_res.0

#Cluster_49 <- BuildClusterTree(Cluster_49, reorder.numeric = TRUE, reorder = TRUE, dims = 1:7)
#PlotClusterTree(object = Cluster_49)

#Cluster_49.markers <- FindAllMarkers(Cluster_49)
#m <- Cluster_49.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_49, features = m$gene)
#write.csv(Cluster_49.markers, "TMS_T_cells_Cluster_49_merkers.csv")

Cluster_49 <- RenameIdents(Cluster_49, "0" = "DN T")

saveRDS(Cluster_49, file = "TMS_merge_T_cells_Cluster_49_clustered.rds")
Cluster_49 <- readRDS(file = "TMS_merge_T_cells_Cluster_49_clustered.rds")

