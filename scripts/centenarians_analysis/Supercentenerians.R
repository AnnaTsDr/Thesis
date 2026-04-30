library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(ggplot2)

setwd("D:/Anna/Centenarians")

genes <- read.csv(file = "genes.csv", row.names = 1)
Supercentenarians.raw <- read.table(file = "01.UMI.txt.gz")
row.names(Supercentenarians.raw) <- make.unique(genes[,2])
#Supercentenarians.raw <- read.csv(file = "UMI+gene_names.csv", sep = ",", row.names = 1)
#Supercentenarians.raw <- read.csv(file = "UMI.csv", sep = ";", row.names = 1)
#writeMM(as.sparse(Supercentenarians.raw), file = "matrix.mtx")
#Supercentenarians.raw <- Read10X(data.dir = "D:/Anna/Centenarians/")
#rownames(Supercentenarians.raw) <- genes[,2]

saveRDS(Supercentenarians.raw, "Supercentenarians_raw_with_gene_symbols.rds")
Supercentenarians.raw <- readRDS("Supercentenarians_raw_with_gene_symbols.rds")

Supercentenarians <- CreateSeuratObject(counts = Supercentenarians.raw, project = "Supercentenarians")

saveRDS(Supercentenarians, "Supercentenarians_SeuratObject.rds")
Supercentenarians <- readRDS("Supercentenarians_SeuratObject.rds")

Metadata.Supercentenerians <- read.table(file = "Centenarians/03.Cell.Barcodes.txt", row.names = 1)

Supercentenarians <- AddMetaData(Supercentenarians, metadata = Metadata.Supercentenerians)

Supercentenarians[["percent.MT"]] <- PercentageFeatureSet(Supercentenarians, pattern = "^MT-")
Supercentenarians[['percent.RIBO']] <- PercentageFeatureSet(Supercentenarians, pattern = "^RP[SL]")
VlnPlot(Supercentenarians, features = c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.RIBO"), ncol = 2)# + ggtitle(label = c("nGenes", "nUMI", "Mito percrnt"))

plot1 <- FeatureScatter(Supercentenarians, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(Supercentenarians, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

max(Supercentenarians@meta.data$nCount_RNA)

Supercentenarians <- subset(Supercentenarians, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.MT < 10)

Supercentenarians <- NormalizeData(Supercentenarians, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(Supercentenarians, "Supercentenarians_NormalizeData.rds")
Supercentenarians <- readRDS("Supercentenarians_NormalizeData.rds")

Supercentenarians <- FindVariableFeatures(Supercentenarians, selection.method = "mean.var.plot", 
                                          dispersion.cutoff = c(0, Inf), mean.cutoff = c(0.0125, 4))

top10 <- head(VariableFeatures(Supercentenarians), 10)
plot1 <- VariableFeaturePlot(Supercentenarians)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(Supercentenarians)
Supercentenarians <- ScaleData(Supercentenarians, features = all.genes)

Supercentenarians <- RunPCA(Supercentenarians, npcs = 24)

DimPlot(Supercentenarians, reduction = "pca")

Supercentenarians <- JackStraw(Supercentenarians, num.replicate = 100, dims = 24)
Supercentenarians <- ScoreJackStraw(Supercentenarians, dims = 1:24)
JackStrawPlot(Supercentenarians, dims = 1:24)

ElbowPlot(Supercentenarians, ndims = 24)

Supercentenarians <- RunTSNE(object = Supercentenarians, dims = 1:24, perplexity = 30)
TSNEPlot(object = Supercentenarians, label = TRUE, pt.size = 0.5)

Supercentenarians <- RunUMAP(object = Supercentenarians, dims = 1:24)
DimPlot(Supercentenarians, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(Supercentenarians, file = "Supercentenarians&tSNE.rds")
Supercentenarians <- readRDS("Supercentenarians&tSNE.rds")

Supercentenarians <- FindNeighbors(Supercentenarians, dims = 1:24)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
Supercentenarians <- FindClusters(Supercentenarians, resolution = res)

clustree(Supercentenarians, prefix = "RNA_snn_res.")
clustree(Supercentenarians, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD14",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD3D",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD3E",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD3G",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "TRAC",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "MS4A1",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD19",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "KLRF1",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "FCGR3A",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "HBA1",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD8A",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD8B",
         exprs = 'scale.data') + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(Supercentenarians) <- Supercentenarians$RNA_snn_res.1.2

Supercentenarians <- BuildClusterTree(Supercentenarians, reorder.numeric = TRUE, reorder = TRUE, dims = 1:24)
PlotClusterTree(object = Supercentenarians)

Supercentenarians.markers_CT_vs_SC <- list()

for (n in 1:29) {
  Supercentenarians.markers_CT_vs_SC[[n]] <- FindConservedMarkers(Supercentenarians, ident.1 = n, grouping.var = "V3")
  write.csv(Supercentenarians.markers_CT_vs_SC[[n]], sprintf("Supercentenarians_markers_CT_vs_SC_%s.csv", n))
}

saveRDS(Supercentenarians, "Supercentenarians_clustered_without_names.rds")
Supercentenarians <- readRDS("Supercentenarians_clustered_without_names.rds")

Supercentenarians.all_markers <- FindAllMarkers(Supercentenarians)

write.csv(Supercentenarians.all_markers, "Supercentenarians_all_markers.csv")

Supercentenarians <- RenameIdents(Supercentenarians, "29" = "27", "28" = "27", "18" = "10", "17" = "10", "16" = "10",
                                  "15" = "10", "14" = "10", "12" = "10", "11" = "10", "21" = "20", 
                                  "26" = "20", "25" = "20", "24" = "20", "23" = "20", "22" = "20", "6" = "3", 
                                  "5" = "3", "4" = "3", "13" = "1")#19 - mki67, 9 - erythro, 2 - pdc, 1 - Mk, 7 - cd16, 8 - mdc

Supercentenarians <- BuildClusterTree(Supercentenarians, reorder.numeric = TRUE, reorder = TRUE, dims = 1:24)
PlotClusterTree(object = Supercentenarians)

Supercentenarians.markers <- list()

for (n in 1:10) {
  Supercentenarians.markers[[n]] <- FindMarkers(Supercentenarians, ident.1 = n, only.pos = TRUE)
  write.csv(Supercentenarians.markers[[n]], sprintf("Supercentenarians_markers_%s_clusters_merged.csv", n))
}

new.cluster.ids <- c("pDC", "mDC", "CD14+ Mono", "CD16+ Mono", "B", "Erythro", "Mk", "MKI67+", "T", "NK")
names(new.cluster.ids) <- levels(Supercentenarians)
Supercentenarians <- RenameIdents(Supercentenarians, new.cluster.ids)

saveRDS(Supercentenarians, file = "Supercentenarians_clutered.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_clutered.rds")

Supercentenarians$cluster_names_pbms <- Supercentenarians@active.ident

Supercentenarians_verification <- clustify(Supercentenarians, cbmc_ref, cluster_col = "cluster_names_pbms", 
                                           obj_out = FALSE)
plot_cor_heatmap(Supercentenarians_verification)

#re-clustering of t cells
Supercentenarians <- subset(Supercentenarians, subset = cluster_names_pbms == "T")

saveRDS(Supercentenarians, file = "Supercentenarians_T_cells.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_T_cells.rds")

Supercentenarians <- subset(Supercentenarians, subset = CD8A == 0)
Supercentenarians <- subset(Supercentenarians, subset = TRDC == 0)
Supercentenarians <- subset(Supercentenarians, subset = CD8B == 0)

Supercentenarians <- FindVariableFeatures(Supercentenarians, selection.method = "mean.var.plot", 
                                          dispersion.cutoff = c(0.4, Inf), mean.cutoff = c(0.0125, 4))

top10 <- head(VariableFeatures(Supercentenarians), 10)
plot1 <- VariableFeaturePlot(Supercentenarians)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(Supercentenarians)
Supercentenarians <- ScaleData(Supercentenarians, features = all.genes)

Supercentenarians <- RunPCA(Supercentenarians, npcs = 20)

DimPlot(Supercentenarians, reduction = "pca")

Supercentenarians <- JackStraw(Supercentenarians, num.replicate = 100, dims = 20)
Supercentenarians <- ScoreJackStraw(Supercentenarians, dims = 1:20)
JackStrawPlot(Supercentenarians, dims = 1:20)

ElbowPlot(Supercentenarians, ndims = 20)

Supercentenarians <- RunTSNE(object = Supercentenarians, dims = 1:20, perplexity = 30)
TSNEPlot(object = Supercentenarians, label = TRUE, pt.size = 0.5)

Supercentenarians <- RunUMAP(object = Supercentenarians, dims = 1:20)
DimPlot(Supercentenarians, label = TRUE, pt.size = 0.5, reduction = "umap")

Supercentenarians <- FindNeighbors(Supercentenarians, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
Supercentenarians <- FindClusters(Supercentenarians, resolution = res)

clustree(Supercentenarians, prefix = "RNA_snn_res.")
clustree(Supercentenarians, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

saveRDS(Supercentenarians, file = "Supercentenarians_T_cells_clustered.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_T_cells_clustered.rds")

Idents(Supercentenarians) <- Supercentenarians$RNA_snn_res.0.8

Supercentenarians <- BuildClusterTree(Supercentenarians, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = Supercentenarians)

#Supercentenarians.clusters <- list()

#for (n in 1:12) {
#  Supercentenarians.clusters[[n]] <- subset(Supercentenarians, idents = n)
#  saveRDS(Supercentenarians.clusters[[n]], sprintf("Supercentenarians_T_cells_clusters_%s.rds", n))
#}

Supercentenarians.markers <- FindAllMarkers(Supercentenarians)
write.csv(Supercentenarians.markers, "Supercentenarians_markers_T_without_CD8A.csv")

Supercentenarians.markers_CT_SC <- list()

for (n in 1:13) {
  Supercentenarians.markers_CT_SC[[n]] <- FindConservedMarkers(Supercentenarians, ident.1 = n, only.pos = TRUE, 
                                                         grouping.var = "V3", min.cells.group = 0, min.cells.feature = 0)
  write.csv(Supercentenarians.markers_CT_SC[[n]], sprintf("Supercentenarians_T_cells_markers_%s_CT_SC_without_CD8A.csv", n))
}

top10 <- Supercentenarians.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Supercentenarians, features = top10$gene)

saveRDS(Supercentenarians, file = "Supercentenarians_T_cells_clusters.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_T_cells_clusters.rds")

Supercentenarians <- RenameIdents(Supercentenarians, "13" = "10", "12" = "10", "5" = "1", "4" = "1", "3" = "1", 
                                  "2" = "1", "6" = "1", "9" = "8", "11" = "10")

Supercentenarians <- RenameIdents(Supercentenarians, "10" = "Naive", "1" = "Cytotoxic", "8" = "TEM", 
                                  "7" = "Tregs")

Supercentenarians <- BuildClusterTree(Supercentenarians, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = Supercentenarians)

saveRDS(Supercentenarians, file = "Supercentenarians_T_cells_clusters_names.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_T_cells_clusters_names.rds")

Supercentenarians$clusters_names_T_cells <- Supercentenarians@active.ident

Supercentenarians.all_markers <- FindAllMarkers(Supercentenarians)
m <- Supercentenarians.all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Supercentenarians, features = m$gene)
write.csv(Supercentenarians.markers, "Supercentenarians_markers_T_without_CD8A_names.csv")

saveRDS(Supercentenarians, file = "Supercentenarians_CD4.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_CD4.rds")

Supercentenarians_Cytotoxic_CD4 <- subset(Supercentenarians, subset = clusters_names_T_cells == "Cytotoxic")

saveRDS(Supercentenarians_Cytotoxic_CD4, "Supercentenarians_Cytotoxic_CD4.rds")

#monocle3
CDS <- as.cell_data_set(Supercentenarians)
CDS <- cluster_cells(cds = CDS, reduction_method = "UMAP")
plot_cells(cds = CDS, label_groups_by_cluster = FALSE, cell_size = 1, color_cells_by = "clusters_names_T_cells")
CDS <- learn_graph(CDS, use_partition = FALSE)
cell_ids <- colnames(CDS)[CDS$clusters_names_T_cells ==  "Naive"]
closest_vertex <- CDS@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(CDS), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(CDS)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]
CDS <- order_cells(cds = CDS, root_pr_nodes = root_pr_nodes)
plot_cells(CDS, color_cells_by = "pseudotime", graph_label_size = 4)
plot_cells(CDS, color_cells_by = "clusters_names_T_cells", group_label_size = 5, graph_label_size = 4)
de_res <- graph_test(CDS, neighbor_graph = "principal_graph", cores = 3)

saveRDS(CDS, "Supercentenarians_CDS.rds")

