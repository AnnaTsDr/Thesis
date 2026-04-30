library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(monocle3)
library(clustree)
library(SeuratWrappers)
library(base)
library(ranger)
library(scater)
library(clustifyrdata)
library(clustifyr)

setwd("D:/Anna/HezisData")

library(MazamaCoreUtils)
library(pbapply)
library(prob)
source(file = "AssessNodesAdaptedfromSeurat.R")

raw_counts <- read.csv(file = "RawData1.csv", sep = ",", row.names = 1, header = TRUE)
object.size(raw_counts)
dim(raw_counts)

agingdata <- CreateSeuratObject(counts = raw_counts, project = "rawdata1", min.cells = 3, min.features = 200)
object.size(agingdata)

saveRDS(agingdata, file = "AgingDataSeuratObject.rds")
agingdata <- readRDS("AgingDataSeuratObject.rds")

Idents(object = agingdata) <- "CD4_T_Cells"
agingdata[["percent.mt"]] <- PercentageFeatureSet(agingdata, pattern = "^mt-")
agingdata[['percent.ribo']] <- PercentageFeatureSet(agingdata, pattern = "^Rp[sl]")
p1 <- VlnPlot(agingdata, features = "nFeature_RNA", pt.size = 0) + labs(title = "nGenes", tag = "A") + NoLegend()
p2 <- VlnPlot(agingdata, features = "nCount_RNA", pt.size = 0) + labs(title = "nUMI", tag = "B") + NoLegend()
p3 <- VlnPlot(agingdata, features = "percent.mt", pt.size = 0) + labs(title = "Mito percent", tag = "C") + NoLegend()
p4 <- VlnPlot(agingdata, features = "percent.ribo", pt.size = 0) + labs(title = "Ribo percent", tag = "C") + NoLegend()
p1 <- AugmentPlot(plot = p1) + labs(y = "nGenes")
p2 <- AugmentPlot(plot = p2) + labs(y = "nUMI")
p3 <- AugmentPlot(plot = p3) + labs(y = "Mito percrnt")
p4 <- AugmentPlot(plot = p4) + labs(y = "Ribo percrnt")
p1 + p2 + p3 + p4
VlnPlot(agingdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) + ggtitle(label = c("nGenes", "nUMI", "Mito percrnt"))

Metadata.aging <- read.csv(file = "Metadata.csv")

Mouse <- as.array(Metadata.aging$Mouse)
Batch <- as.array(Metadata.aging$Batch)
Age_group <- as.array(Metadata.aging$Age_group)
Subset <- as.array(Metadata.aging$Subset)

agingdata <- AddMetaData(agingdata, metadata = Mouse, col.name = "Mouse")
agingdata <- AddMetaData(agingdata, metadata = Batch, col.name = "Batch")
agingdata <- AddMetaData(agingdata, metadata = Age_group, col.name = "Age_group")
agingdata <- AddMetaData(agingdata, metadata = Subset, col.name = "Subset_Idan")

saveRDS(agingdata, file = "AgingDataMetadataAdded.rds")
agingdata <- readRDS("AgingDataMetadataAdded.rds")

plot1 <- FeatureScatter(agingdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(agingdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

agingdata <- NormalizeData(agingdata, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(agingdata, file = "AgingDataNormalized.rds")
agingdata <- readRDS(file = "AgingDataNormalized.rds")

agingdata <- FindVariableFeatures(agingdata, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                                  mean.cutoff = c(0.0125, 2), mean.function = ExpMean, dispersion.function = LogVMR)

top10 <- head(VariableFeatures(agingdata), 10)
plot1 <- VariableFeaturePlot(agingdata)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

saveRDS(agingdata, file = "AgingDataHVG.rds")
write.csv(VariableFeatures(agingdata), file = "AgingDataHVG.csv")
agingdata <- readRDS(file = "AgingDataHVG.rds")

all.genes <- rownames(agingdata)
agingdata <- ScaleData(agingdata, features = all.genes)

saveRDS(agingdata, file = "AgingDataScaled.rds")
agingdata <- readRDS(file = "AgingDataScaled.rds")

agingdata <- RunPCA(agingdata, npcs = 20)

DimPlot(agingdata, reduction = "pca")

agingdata <- JackStraw(agingdata, num.replicate = 100, dims = 20)
agingdata <- ScoreJackStraw(agingdata, dims = 1:20)
JackStrawPlot(agingdata, dims = 1:20)

saveRDS(agingdata, file = "AgingDataJackStraw.rds")
agingdata <- readRDS(file = "AgingDataJackStraw.rds")

ElbowPlot(agingdata, ndims = 20)

agingdata <- RunTSNE(object = agingdata, dims = 1:20, perplexity = 30)
TSNEPlot(object = agingdata, label = TRUE, pt.size = 0.5)
FeaturePlot(agingdata, features = "nCount_RNA")

agingdata <- RunUMAP(agingdata, dims = 1:20)
DimPlot(agingdata, reduction = "umap")

saveRDS(agingdata, file = "AgingData&tSNE&UMAP.rds")
agingdata <- readRDS("AgingData&tSNE&UMAP.rds")

agingdata <- FindNeighbors(agingdata, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
agingdata <- FindClusters(agingdata, resolution = res)

clustree(agingdata)

Idents(agingdata) <- agingdata$RNA_snn_res.0.8

agingdata <- BuildClusterTree(agingdata, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = agingdata)

#FindAllMarkers
agingdata.markers.wilcox <- FindAllMarkers(agingdata, test.use = "wilcox", logfc.threshold = 0.4, only.pos = TRUE)
m <- agingdata.markers.wilcox %>% group_by(cluster) %>% top_n(n = 5, wt = p_val_adj)

#Merge nodes
agingdata <- RenameIdents(agingdata, "12" = "3", "7" = "6", "8" = "6", "9" = "6", "10" = "6")
agingdata <- BuildClusterTree(agingdata, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = agingdata)

agingdata$clusters_reordered <- agingdata@active.ident

#clystify
ref_CD4TcellsAging <- seurat_ref(agingdata, cluster_col = "Subset_Idan")
agingdata <- clustify(input = agingdata, ref_mat = ref_CD4TcellsAging, cluster_col = "RNA_snn_res.0.8")
agingdata_lists <- clustify(agingdata, ref_mat = ref_CD4TcellsAging, cluster_col = "clusters_reordered", 
                         seurat_out = FALSE)
write.csv(agingdata_lists, "agingdata_lists.csv")

plot_cor_heatmap(agingdata_lists)

#modify identetis
current.cluster.ids <- c(3, 6, 7, 1, 2, 5, 4)
new.cluster.ids <- c("Cytotoxic", "TEM", "rTreg", "aTreg", "Exhausted", "Naive_Isg15", "Naive")
Idents(object = agingdata) <- plyr::mapvalues(x = Idents(object = agingdata), from = current.cluster.ids, 
                                              to = new.cluster.ids)

#DEG
agingdata.old.markers <- FindMarkers(agingdata, group.by = "Age_group", ident.1 = "Old", test.use = "roc", logfc.threshold = 0.4, only.pos = TRUE)
agingdata.young.markers <- FindMarkers(agingdata, group.by = "Age_group", ident.1 = "Young", test.use = "roc", logfc.threshold = 0.4, only.pos = TRUE)
agingdata.conserved.markers_Cytotoxic <- FindConservedMarkers(agingdata, ident.1 = "Cytotoxic", grouping.var = "Age_group")
agingdata.conserved.markers_TEM <- FindConservedMarkers(agingdata, ident.1 = "TEM", grouping.var = "Age_group")
agingdata.conserved.markers_rTreg <- FindConservedMarkers(agingdata, ident.1 = "rTreg", grouping.var = "Age_group")
agingdata.conserved.markers_aTreg <- FindConservedMarkers(agingdata, ident.1 = "aTreg", grouping.var = "Age_group")
agingdata.conserved.markers_Exhausted <- FindConservedMarkers(agingdata, ident.1 = "Exhausted", grouping.var = "Age_group")
agingdata.conserved.markers_Naive_Isg15 <- FindConservedMarkers(agingdata, ident.1 = "Naive_Isg15", grouping.var = "Age_group")
agingdata.conserved.markers_Naive <- FindConservedMarkers(agingdata, ident.1 = "Naive", grouping.var = "Age_group")

write.csv(agingdata.old.markers, "agingdata_old_markers.csv")
write.csv(agingdata.young.markers, "agingdata_young_markers.csv")
write.csv(agingdata.conserved.markers_Cytotoxic, "agingdata_conserved_markers_Cytotoxic.csv")
write.csv(agingdata.conserved.markers_TEM, "agingdata_conserved_markers_TEM.csv")
write.csv(agingdata.conserved.markers_rTreg, "agingdata_conserved_markers_rTreg.csv")
write.csv(agingdata.conserved.markers_aTreg, "agingdata_conserved_markers_aTreg.csv")
write.csv(agingdata.conserved.markers_Exhausted, "agingdata_conserved_markers_Exhausted.csv")
write.csv(agingdata.conserved.markers_Naive, "agingdata_conserved_markers_Naive_Isg15.csv")
write.csv(agingdata.conserved.markers_Naive_Isg15, "agingdata_conserved_markers_Naive_Isg15.csv")

#FindAllMarkers
agingdata.markers.wilcox <- FindAllMarkers(agingdata, test.use = "wilcox", logfc.threshold = 0.4, only.pos = TRUE)
m <- agingdata.markers.wilcox %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(agingdata, features = m$gene) + NoLegend()

write.csv(agingdata.markers.wilcox, "agingdata_markers_wilcox.csv")

agingdata$Subset_Anna <- agingdata@active.ident

CTL_CD4_Subset_agingdata <- subset(agingdata, subset = Subset_Anna == "Cytotoxic")
saveRDS(CTL_CD4_Subset_agingdata, file = "CTL_CD4_Subset_agingdata.csv")

agingdata <- RunUMAP(agingdata, dims = 1:20)
DimPlot(agingdata, reduction = "umap")

saveRDS(agingdata, file = "AgingDataClustered.rds")
agingdata <- readRDS(file = "AgingDataClustered.rds")

#monocle3
CDS <- as.cell_data_set(agingdata)
CDS <- cluster_cells(cds = CDS, reduction_method = "UMAP")
plot_cells(cds = CDS, label_groups_by_cluster = FALSE, cell_size = 1, color_cells_by = "ident")
CDS <- learn_graph(CDS, use_partition = FALSE)
plot_cells(cds = CDS, cell_size = 1, color_cells_by = "ident")
cell_ids <- colnames(CDS)[CDS$ident ==  "Na?ve"]
closest_vertex <- CDS@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(CDS), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(CDS)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]
CDS <- order_cells(cds = CDS, root_pr_nodes = root_pr_nodes)
plot_cells(CDS, color_cells_by = "pseudotime", graph_label_size = 4)
plot_cells(CDS, color_cells_by = "ident", group_label_size = 5, graph_label_size = 4)
de_res <- graph_test(CDS, neighbor_graph = "principal_graph", cores = 3)

saveRDS(CDS, "aging_CDS.rds")
CDS <- readRDS("aging_CDS.rds")

#CTL_Subset
VlnPlot(CTL_Subset, features = c("Cdkn2a", "Il27ra"), group.by = "Age_group")
CTL_Subset <- FindVariableFeatures(CTL_Subset, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 2), mean.function = ExpMean, dispersion.function = LogVMR)


DotPlot(agingdata, features = c("Il27ra", "Il27"), group.by = "Age_group", cols = c("blue", "yellow"))
RidgePlot(agingdata, features = "Il27ra", group.by = "Age_group")
VlnPlot(agingdata, features = "Il27ra", group.by = "Age_group")
DotPlot(Old_Mice, features = "Il27ra", group.by = "Subset", cols = c("blue", "yellow"))
DotPlot(Cytotoxic, features = "Il27ra", group.by = "Age_group", cols = c("blue", "yellow"))
VlnPlot(agingdata, features = c("Il27", "IL27ra"), group.by = "Subset",pt.size = 0)
VlnPlot(agingdata, features = c("Il27", "Il27ra"), group.by = "Age_group",pt.size = 0, ncol = 1)
VlnPlot(agingdata, features = c("Tnf"), group.by = "Age_group",pt.size = 0, ncol = 1)
VlnPlot(agingdata, features = c("Gzmk"), group.by = "Subset",pt.size = 0, ncol = 1)
VlnPlot(agingdata, features = c("Tnf"), group.by = "Age_group",pt.size = 0, ncol = 1)
VlnPlot(agingdata, features = c("Tnf"), group.by = "Subset",pt.size = 0, ncol = 1)
DotPlot(agingdata, features = "Tnf", group.by = "Subset")
DotPlot(Old_Mice, features = "Il27ra", group.by = "Subset")
DotPlot(Old_Mice, features = "Il27ra", group.by = "Subset", cols = c("blue", "yellow"))
VlnPlot(Cytotoxic, features = "Il27ra", group.by = "Age_group")
VlnPlot(Cytotoxic, features = "Il27ra")
DotPlot(Cytotoxic, features = "Il27ra", cols = c("blue", "yellow"))
VlnPlot(Cytotoxic, features = "Il27ra", group.by = "Age_group")
DotPlot(Cytotoxic, features = "Il27ra", cols = c("blue", "yellow"))
DotPlot(Old_Mice, features = "Il27ra", cols = c("blue", "yellow"))
RidgePlot(Cytotoxic, features = "Il27ra")
RidgePlot(agingdata, features = "Il27ra")
RidgePlot(Old_Mice, features = "Il27ra")
RidgePlot(Old_Mice, features = "Il27ra", group.by = "Subset")
RidgePlot(agingdata, features = "Il27ra", group.by = "Age_group")
RidgePlot(agingdata, features = "Il27ra", group.by = "Subset")
RidgePlot(Cytotoxic, features = "Il27ra", group.by = "Age_group")
RidgePlot(Young_Mice, features = "Il27ra", group.by = "Subset")
RidgePlot(Young_Mice, features = "Il27ra")
DotPlot(Young_Mice, features = "Il27ra", cols = c("blue", "yellow"))
DotPlot(Young_Mice, features = "Il27ra", cols = c("blue", "yellow"), group.by = "Subset")
DotPlot(Cytotoxic_Young, features = "Il27ra", cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = "Il27ra", cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Tnf", "Il1", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Tnf", "Il1a", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Il1b", "Il1a", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Tp53", "Il1a", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Ifn", "Il1a", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Ifng", "Il1a", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Ifng", "Il1nr", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Ifng", "Il1rn", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Ifng", "Il6", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Ifng", "Tnfa", "Cdkn2a"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Ifng", "Cdkn2b"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Ifng", "Cdkn1a"), cols = c("blue", "yellow"))
DotPlot(agingdata, features = c("Il27ra", "Ifng", "Cdkn1a", "Tnf", "Cdkn2a", "Cdkn2b"), cols = c("blue", "yellow"), group.by = "Subset")
DotPlot(agingdata, features = c("Il27ra", "Ifng", "Cdkn1a", "Tnf", "Cdkn2a", "Cdkn2b"), cols = c("blue", "yellow"), group.by = "Age_group")
DotPlot(Old_Mice, features = c("Il27ra", "Ifng", "Cdkn1a", "Tnf", "Cdkn2a", "Cdkn2b"), cols = c("blue", "yellow"), group.by = "Subset")
DotPlot(Cytotoxic, features = c("Il27ra", "Ifng", "Cdkn1a", "Tnf", "Cdkn2a", "Cdkn2b"), cols = c("blue", "yellow"), group.by = "Age_group")
DotPlot(Cytotoxic_Young, features = c("Il27ra", "Ifng", "Cdkn1a", "Tnf", "Cdkn2a", "Cdkn2b"), cols = c("blue", "yellow"))
DotPlot(Cytotoxic_Old, features = c("Il27ra", "Ifng", "Cdkn1a", "Tnf", "Cdkn2a", "Cdkn2b"), cols = c("blue", "yellow"))
DotPlot(agingdata, features = c("Il27ra", "Ifng", "Cdkn1a", "Tnf", "Cdkn2a", "Cdkn2b"), cols = c("blue", "yellow"))
DoHeatmap(agingdata, features = VariableFeatures(agingdata), group.by = "Age_group", group.bar = TRUE)
DoHeatmap(agingdata, features = VariableFeatures(agingdata), group.by = "Subset")

Old.markers <- FindMarkers(agingdata, ident.1 = c("o1", "o2", "o3", "o4"), ident.2 = c("o1", "o2", "o3", "o4"))



#top10 <- Old.markers %>% group_by(Age_group) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(agingdata, features = rownames(head(Old.markers, 10)), group.by = "Age_group") + NoLegend()

#scater
agingdata.sce <- as.SingleCellExperiment(agingdata)
p1 <- plotExpression(agingdata.sce, features = "Gzmk", x = "ident") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- plotPCA(agingdata.sce, colour_by = "ident")
CombinePlots(plots = list(p1, p2))

#here should be an assessnode dont know whats not working needs to fix
node.scores <- AssessNodes(object = agingdata)#, genes.training = VariableFeatures(agingdata)
node.scores[order(node.scores$oobe,decreasing = TRUE),] -> node.scores
nodes.merge <- node.scores[which(node.scores$oobe > 0.05),]
nodes.to.merge <- sort(nodes.merge$node)
agingdata.merged <- agingdata

#Get seurat v2 clusters and set clusters
seurat_v2_clusters <- read.csv(file = "AgingDataMetaData_res.csv", sep = ";")
seurat_v2_clusters <- as.array(seurat_v2_clusters$seurat_v2_clusters)
agingdata <- AddMetaData(agingdata, metadata = seurat_v2_clusters, col.name = "seurat_v2_clusters")
agingdata$seurat_clusters <- agingdata$seurat_v2_clusters
#agingdata <- BuildClusterTree(agingdata, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
agingdata$tree.ident <- agingdata$seurat_v2_clusters
Idents(agingdata) <- agingdata$seurat_v2_clusters
PlotClusterTree(object = agingdata)

Idents(agingdata) <- agingdata$Subset_Idan
agingdata.markers_Cytotoxic <- FindMarkers(agingdata, ident.1 = "Cytotoxic", 
                                           logfc.threshold = -Inf, grouping.var = "Subset_Idan")
# add a column of NAs
agingdata.markers_Cytotoxic$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
agingdata.markers_Cytotoxic$diffexpressed[agingdata.markers_Cytotoxic$avg_logFC > 0.4 & agingdata.markers_Cytotoxic$p_val_adj < 0.05] <- "Cytotoxic"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
agingdata.markers_Cytotoxic$diffexpressed[agingdata.markers_Cytotoxic$avg_logFC < -0.4 & agingdata.markers_Cytotoxic$p_val_adj < 0.05] <- "CD4 T cells not Cytotoxic"
agingdata.markers_Cytotoxic$delabel <- NA
agingdata.markers_Cytotoxic$delabel[agingdata.markers_Cytotoxic$diffexpressed != "NO"] <- rownames(agingdata.markers_Cytotoxic)[agingdata.markers_Cytotoxic$diffexpressed != "NO"]
ggplot(agingdata.markers_Cytotoxic, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(agingdata.markers_Cytotoxic, "agingdata_markers_Cytotoxic.csv")

# add a column of NAs
CD4_T_Cytotoxic_vs_other $diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_T_Cytotoxic_vs_other $diffexpressed[CD4_T_Cytotoxic_vs_other $avg_logFC > 0.4 & CD4_T_Cytotoxic_vs_other $p_val_adj < 0.05] <- "Cytotoxic"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_T_Cytotoxic_vs_other $diffexpressed[CD4_T_Cytotoxic_vs_other $avg_logFC < -0.4 & CD4_T_Cytotoxic_vs_other $p_val_adj < 0.05] <- "CD4 T cells not Cytotoxic"
CD4_T_Cytotoxic_vs_other $delabel <- NA
CD4_T_Cytotoxic_vs_other $delabel[CD4_T_Cytotoxic_vs_other $diffexpressed != "NO"] <- rownames(CD4_T_Cytotoxic_vs_other )[CD4_T_Cytotoxic_vs_other $diffexpressed != "NO"]
ggplot(CD4_T_Cytotoxic_vs_other , aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(CD4_T_Cytotoxic_vs_other , "CD4_TMS_markers_Cytotoxic.csv")

agingdata.markers_Cytotoxic_vs_TEM <- FindMarkers(agingdata, ident.1 = "Cytotoxic", ident.2 = "TEM",
                                           logfc.threshold = -Inf, grouping.var = "Subset_Idan")
# add a column of NAs
agingdata.markers_Cytotoxic_vs_TEM$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
agingdata.markers_Cytotoxic_vs_TEM$diffexpressed[agingdata.markers_Cytotoxic_vs_TEM$avg_logFC > 0.4 & 
                                                   agingdata.markers_Cytotoxic_vs_TEM$p_val_adj < 0.05] <- "Cytotoxic"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
agingdata.markers_Cytotoxic_vs_TEM$diffexpressed[agingdata.markers_Cytotoxic_vs_TEM$avg_logFC < -0.4 & 
                                       agingdata.markers_Cytotoxic_vs_TEM$p_val_adj < 0.05] <- "TEM"
agingdata.markers_Cytotoxic_vs_TEM$delabel <- NA
agingdata.markers_Cytotoxic_vs_TEM$delabel[agingdata.markers_Cytotoxic_vs_TEM$diffexpressed != "NO"] <- 
  rownames(agingdata.markers_Cytotoxic_vs_TEM)[agingdata.markers_Cytotoxic_vs_TEM$diffexpressed != "NO"]
ggplot(agingdata.markers_Cytotoxic_vs_TEM, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+
  geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") + 
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")

agingdata.markers_Exhausted_vs_TEM <- FindMarkers(agingdata, ident.1 = "Exhausted", ident.2 = "TEM",
                                                  logfc.threshold = -Inf, grouping.var = "Subset_Idan")
# add a column of NAs
agingdata.markers_Exhausted_vs_TEM$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
agingdata.markers_Exhausted_vs_TEM$diffexpressed[agingdata.markers_Exhausted_vs_TEM$avg_logFC > 0.4 & 
                                                   agingdata.markers_Exhausted_vs_TEM$p_val_adj < 0.05] <- "Exhausted"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
agingdata.markers_Exhausted_vs_TEM$diffexpressed[agingdata.markers_Exhausted_vs_TEM$avg_logFC < -0.4 & 
                                       agingdata.markers_Exhausted_vs_TEM$p_val_adj < 0.05] <- "TEM"
agingdata.markers_Exhausted_vs_TEM$delabel <- NA
agingdata.markers_Exhausted_vs_TEM$delabel[agingdata.markers_Exhausted_vs_TEM$diffexpressed != "NO"] <- 
  rownames(agingdata.markers_Exhausted_vs_TEM)[agingdata.markers_Exhausted_vs_TEM$diffexpressed != "NO"]
ggplot(agingdata.markers_Exhausted_vs_TEM, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+
  geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") + geom_hline(yintercept = -log10(0.05), col="red")+ 
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")

agingdata.markers_aTregs_vs_rTregs <- FindMarkers(agingdata, ident.1 = "aTregs", ident.2 = "rTregs",
                                                  logfc.threshold = -Inf, grouping.var = "Subset_Idan")
# add a column of NAs
agingdata.markers_aTregs_vs_rTregs$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
agingdata.markers_aTregs_vs_rTregs$diffexpressed[agingdata.markers_aTregs_vs_rTregs$avg_logFC > 0.4 & 
                                                   agingdata.markers_aTregs_vs_rTregs$p_val_adj < 0.05] <- "aTregs"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
agingdata.markers_aTregs_vs_rTregs$diffexpressed[agingdata.markers_aTregs_vs_rTregs$avg_logFC < -0.4 & 
                                                   agingdata.markers_aTregs_vs_rTregs$p_val_adj < 0.05] <- "rTregs"
agingdata.markers_aTregs_vs_rTregs$delabel <- NA
agingdata.markers_aTregs_vs_rTregs$delabel[agingdata.markers_aTregs_vs_rTregs$diffexpressed != "NO"] <- 
  rownames(agingdata.markers_aTregs_vs_rTregs)[agingdata.markers_aTregs_vs_rTregs$diffexpressed != "NO"]
ggplot(agingdata.markers_aTregs_vs_rTregs, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+
  geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") + geom_hline(yintercept = -log10(0.05), col="red")+ 
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")
