Cluster_66 <- readRDS(file = "TMS_merge_T_cells_clusters_66.rds")

#Highly variable expressed genes
Cluster_66 <- FindVariableFeatures(Cluster_66, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_66),15)
plot1 <- VariableFeaturePlot(Cluster_66)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_66), file = "TMS_merge_T_cells_Cluster_66_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_66, file = "TMS_merge_T_cells_Cluster_66_VariableGenes.rds")
Cluster_66 <- readRDS(file = "TMS_merge_T_cells_Cluster_66_VariableGenes.rds")

all.genes <- rownames(Cluster_66)
Cluster_66 <- ScaleData(Cluster_66, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_66, file = "TMS_merge_T_cells_Cluster_66_ScaledData.rds")
Cluster_66 <- readRDS(file = "TMS_merge_T_cells_Cluster_66_ScaledData.rds")

Cluster_66 <- RunPCA(Cluster_66, ndims.print = 1:9, nfeatures.print = 5)
DimPlot(Cluster_66, reduction = "pca")

saveRDS(Cluster_66, file = "TMS_merge_T_cells_Cluster_66_PCA.rds")
Cluster_66 <- readRDS(file = "TMS_merge_T_cells_Cluster_66_PCA.rds")

Cluster_66 <- JackStraw(Cluster_66, num.replicate = 100, dims = 9)
Cluster_66 <- ScoreJackStraw(Cluster_66, dims = 1:9)
JackStrawPlot(Cluster_66, dims = 1:9)

ElbowPlot(Cluster_66, ndims = 9)

Cluster_66 <- RunTSNE(Cluster_66, dims = 1:9, perplexity = 10)
DimPlot(Cluster_66, reduction = "tsne")

Cluster_66 <- RunUMAP(Cluster_66, dims = 1:9)
DimPlot(Cluster_66, reduction = "umap")

saveRDS(Cluster_66, file = "TMS_merge_T_cells_Cluster_66_TSNE_UMAP.rds")
Cluster_66 <- readRDS(file = "TMS_merge_T_cells_Cluster_66_TSNE_UMAP.rds")

Cluster_66 <- FindNeighbors(Cluster_66, dims = 1:9)

#louvian
Cluster_66 <- FindClusters(Cluster_66, resolution = seq(0,1.2, by = 0.1))

clustree(Cluster_66, prefix = "RNA_snn_res.")
clustree(Cluster_66, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_66) <- Cluster_66$RNA_snn_res.0

#Cluster_66 <- BuildClusterTree(Cluster_66, reorder.numeric = TRUE, reorder = TRUE, dims = 1:15)
#PlotClusterTree(object = Cluster_66)

#Cluster_66.markers <- FindAllMarkers(Cluster_66)
#m <- Cluster_66.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#DoHeatmap(Cluster_66, features = m$gene)
#write.csv(Cluster_66.markers, "TMS_T_cells_Cluster_66_merkers.csv")

Cluster_66 <- RenameIdents(Cluster_66, "0" = "DP T")

saveRDS(Cluster_66, file = "TMS_merge_T_cells_Cluster_66_clustered.rds")
Cluster_66 <- readRDS(file = "TMS_merge_T_cells_Cluster_66_clustered.rds")

