Cluster_53 <- readRDS(file = "TMS_merge_T_cells_clusters_53.rds")

#Highly variable expressed genes
Cluster_53 <- FindVariableFeatures(Cluster_53, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(1.2, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(Cluster_53),15)
plot1 <- VariableFeaturePlot(Cluster_53)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Cluster_53), file = "TMS_merge_T_cells_Cluster_53_Highly_Variable_Expessed_Genes.csv")

saveRDS(Cluster_53, file = "TMS_merge_T_cells_Cluster_53_VariableGenes.rds")
Cluster_53 <- readRDS(file = "TMS_merge_T_cells_Cluster_53_VariableGenes.rds")

all.genes <- rownames(Cluster_53)
Cluster_53 <- ScaleData(Cluster_53, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Cluster_53, file = "TMS_merge_T_cells_Cluster_53_ScaledData.rds")
Cluster_53 <- readRDS(file = "TMS_merge_T_cells_Cluster_53_ScaledData.rds")

Cluster_53 <- RunPCA(Cluster_53, ndims.print = 1:17, nfeatures.print = 5)
DimPlot(Cluster_53, reduction = "pca")

saveRDS(Cluster_53, file = "TMS_merge_T_cells_Cluster_53_PCA.rds")
Cluster_53 <- readRDS(file = "TMS_merge_T_cells_Cluster_53_PCA.rds")

Cluster_53 <- JackStraw(Cluster_53, num.replicate = 100, dims = 17)
Cluster_53 <- ScoreJackStraw(Cluster_53, dims = 1:17)
JackStrawPlot(Cluster_53, dims = 1:17)

ElbowPlot(Cluster_53, ndims = 17)

Cluster_53 <- RunTSNE(Cluster_53, dims = 1:17, perplexity = 10)
DimPlot(Cluster_53, reduction = "tsne")

Cluster_53 <- RunUMAP(Cluster_53, dims = 1:17)
DimPlot(Cluster_53, reduction = "umap")

saveRDS(Cluster_53, file = "TMS_merge_T_cells_Cluster_53_TSNE_UMAP.rds")
Cluster_53 <- readRDS(file = "TMS_merge_T_cells_Cluster_53_TSNE_UMAP.rds")

Cluster_53 <- FindNeighbors(Cluster_53, dims = 1:17)

#louvian
Cluster_53 <- FindClusters(Cluster_53, resolution = seq(0,1.22, by = 0.1))

clustree(Cluster_53, prefix = "RNA_snn_res.")
clustree(Cluster_53, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Cluster_53) <- Cluster_53$RNA_snn_res.1.2

Cluster_53 <- BuildClusterTree(Cluster_53, reorder.numeric = TRUE, reorder = TRUE, dims = 1:15)
PlotClusterTree(object = Cluster_53)

Cluster_53.markers <- FindAllMarkers(Cluster_53)
m <- Cluster_53.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Cluster_53, features = m$gene)
write.csv(Cluster_53.markers, "TMS_T_cells_Cluster_53_merkers.csv")

Cluster_53 <- RenameIdents(Cluster_53, "1" = "Not T", "2" = "DN T", "3" = "CD8 T", "4" = "DN T", "5" = "DN T",
                           "6" = "DN T", "7" = "Not T")

saveRDS(Cluster_53, file = "TMS_merge_T_cells_Cluster_53_clustered.rds")
Cluster_53 <- readRDS(file = "TMS_merge_T_cells_Cluster_53_clustered.rds")

