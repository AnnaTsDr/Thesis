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

Supercentenarians <- CreateSeuratObject(counts = Supercentenarians.raw, project = "Supercentenarians", min.cells = 3, min.features = 200)

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

Supercentenarians <- FindVariableFeatures(Supercentenarians, selection.method = "mean.var.plot", 
                                          dispersion.cutoff = c(0.3, Inf), mean.cutoff = c(0.0125, 4))

top10 <- head(VariableFeatures(Supercentenarians), 10)
plot1 <- VariableFeaturePlot(Supercentenarians)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

all.genes <- rownames(Supercentenarians)
Supercentenarians <- ScaleData(Supercentenarians, features = all.genes)

Supercentenarians <- RunPCA(Supercentenarians)

DimPlot(Supercentenarians, reduction = "pca")

Supercentenarians <- JackStraw(Supercentenarians, num.replicate = 100, dims = 30)
Supercentenarians <- ScoreJackStraw(Supercentenarians, dims = 1:30)
JackStrawPlot(Supercentenarians, dims = 1:30)

ElbowPlot(Supercentenarians, ndims = 30)

Supercentenarians <- RunTSNE(object = Supercentenarians, dims = 1:30, perplexity = 30)
TSNEPlot(object = Supercentenarians, label = TRUE, pt.size = 0.5)

saveRDS(Supercentenarians, file = "Supercentenarians&tSNE.rds")
Supercentenarians <- readRDS("Supercentenarians&tSNE.rds")

Supercentenarians <- FindNeighbors(Supercentenarians, dims = 1:30)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
Supercentenarians <- FindClusters(Supercentenarians, resolution = res)

clustree(Supercentenarians, prefix = "RNA_snn_res.")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD14")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD3D")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD3E")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD3G")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "TRAC")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "MS4A1")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD19")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "KLRF1")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "FCGR3A")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "HBA1")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD8A")
clustree(Supercentenarians, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "CD8B")

Idents(Supercentenarians) <- Supercentenarians$RNA_snn_res.1.1

Supercentenarians <- BuildClusterTree(Supercentenarians, reorder.numeric = TRUE, reorder = TRUE, dims = 1:30)
PlotClusterTree(object = Supercentenarians)

Supercentenarians.markers <- list()

for (n in 1:27) {
  Supercentenarians.markers[[n]] <- FindMarkers(Supercentenarians, ident.1 = n, only.pos = TRUE)
  write.csv(Supercentenarians.markers[[n]], sprintf("Supercentenarians_markers_%s.csv", n))
}

Supercentenarians.markers_CT_vs_SC <- list()

for (n in 1:27) {
  Supercentenarians.markers_CT_vs_SC[[n]] <- FindConservedMarkers(Supercentenarians, ident.1 = n, grouping.var = "V3")
  write.csv(Supercentenarians.markers_CT_vs_SC[[n]], sprintf("Supercentenarians_markers_CT_vs_SC_%s.csv", n))
}

saveRDS(Supercentenarians, "Supercentenarians_clustered_without_names.rds")
Supercentenarians <- readRDS("Supercentenarians_clustered_without_names.rds")

Supercentenarians.all_markers <- FindAllMarkers(Supercentenarians)

write.csv(Supercentenarians.all_markers, "Supercentenarians_all_markers.csv")

Supercentenarians <- RenameIdents(Supercentenarians, "27" = "25", "26" = "25", "24" = "16", "23" = "16", "22" = "16",
                                  "21" = "16", "20" = "16", "19" = "16", "18" = "16", "17" = "16", "14" = "9",
                                  "13" = "9", "12" = "9", "11" = "9", "10" = "9", "4" = "2", "3" = "2")

Supercentenarians <- BuildClusterTree(Supercentenarians, reorder.numeric = TRUE, reorder = TRUE, dims = 1:30)
PlotClusterTree(object = Supercentenarians)

Supercentenarians.markers <- list()

for (n in 1:10) {
  Supercentenarians.markers[[n]] <- FindMarkers(Supercentenarians, ident.1 = n, only.pos = TRUE)
  write.csv(Supercentenarians.markers[[n]], sprintf("Supercentenarians_markers_%s_clusters_merged.csv", n))
}

new.cluster.ids <- c( "DC", "moDC", "CD14+ Mono", "FCGR3A+ Mono", "B", "MKI67+", "Erithoid", "Platelet","T", "NK")
names(new.cluster.ids) <- levels(Supercentenarians)
Supercentenarians <- RenameIdents(Supercentenarians, new.cluster.ids)

Supercentenarians_verification <- clustify(Supercentenarians, cbmc_ref, cluster_col = "cluster_names_pbms", 
                                           obj_out = FALSE)
plot_cor_heatmap(Supercentenarians_verification)

saveRDS(Supercentenarians, file = "Supercentenarians_clutered.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_clutered.rds")

Supercentenarians$cluster_names_pbms <- Supercentenarians@active.ident

#reclustering of t cells
Supercentenarians <- subset(Supercentenarians, subset = cluster_names_pbms == "T")

saveRDS(Supercentenarians, file = "Supercentenarians_T_cells.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_T_cells.rds")

Supercentenarians <- FindVariableFeatures(Supercentenarians, selection.method = "mean.var.plot", 
                                          dispersion.cutoff = c(0.45, Inf), mean.cutoff = c(0.0125, 4))

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

Supercentenarians <- FindNeighbors(Supercentenarians, dims = 1:20)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
Supercentenarians <- FindClusters(Supercentenarians, resolution = res)

clustree(Supercentenarians, prefix = "RNA_snn_res.")

Idents(Supercentenarians) <- Supercentenarians$RNA_snn_res.1

Supercentenarians <- BuildClusterTree(Supercentenarians, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = Supercentenarians)

Supercentenarians.markers <- list()

for (n in 1:19) {
  Supercentenarians.markers[[n]] <- FindMarkers(Supercentenarians, ident.1 = n, only.pos = TRUE)
  write.csv(Supercentenarians.markers[[n]], sprintf("Supercentenarians_T_cells_markers_%s.csv", n))
}

saveRDS(Supercentenarians, file = "Supercentenarians_T_cells_clusters.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_T_cells_clusters.rds")



Supercentenarians <- RenameIdents(Supercentenarians, "18" = "6", "17" = "6", "16" = "6", "13" ="6", "11" = "6", 
                                  "19" = "6", "15" = "1", "14" = "1", "12" = "1", "8" = "1", "7" = "1", "3" = "1", 
                                  "4" = "1", "5" = "1", "2" = "1", "9" = "6")

new.cluster.ids <- c( "CD8 T", "CD4 T", "NK")
names(new.cluster.ids) <- levels(Supercentenarians)
Supercentenarians <- RenameIdents(Supercentenarians, new.cluster.ids)

saveRDS(Supercentenarians, file = "Supercentenarians_T_cells_clusters_names.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_T_cells_clusters_names.rds")

Supercentenarians$clusters_names_T_cells <- Supercentenarians@active.ident

Supercentenarians <- subset(Supercentenarians, subset = clusters_names_T_cells == "CD4 T")

saveRDS(Supercentenarians, file = "Supercentenarians_CD4_for_integrated.rds")
Supercentenarians <- readRDS(file = "Supercentenarians_CD4_for_integrated.rds")

Supercentenarians <- BuildClusterTree(Supercentenarians, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = Supercentenarians)

Supercentenarians.markers <- FindAllMarkers(Supercentenarians)
Supercentenarians <- BuildClusterTree(Supercentenarians, dims = 1:20)
PlotClusterTree(object = Supercentenarians)

Supercentenarians.wilcox <- FindAllMarkers(Supercentenarians, test.use = "wilcox", only.pos = TRUE, return.thresh = 0.001, logfc.threshold = 0.4)
Supercentenarians.bimod <- FindAllMarkers(Supercentenarians, test.use = "bimod", only.pos = TRUE, return.thresh = 0.001, logfc.threshold = 0.4)
Supercentenarians.roc <- FindAllMarkers(Supercentenarians, test.use = "roc", only.pos = TRUE, return.thresh = 0.001, logfc.threshold = 0.4)
Supercentenarians.t <- FindAllMarkers(Supercentenarians, test.use = "t", only.pos = TRUE, return.thresh = 0.001, logfc.threshold = 0.4)
Supercentenarians.negbinom <- FindAllMarkers(Supercentenarians, test.use = "negbinom", only.pos = TRUE, return.thresh = 0.001, logfc.threshold = 0.4)
Supercentenarians.poisson <- FindAllMarkers(Supercentenarians, test.use = "poisson", only.pos = TRUE, return.thresh = 0.001, logfc.threshold = 0.4)
Supercentenarians.LR <- FindAllMarkers(Supercentenarians, test.use = "LR", only.pos = TRUE, return.thresh = 0.001, logfc.threshold = 0.4)
Supercentenarians.MAST <- FindAllMarkers(Supercentenarians, test.use = "MAST", only.pos = TRUE, return.thresh = 0.001, logfc.threshold = 0.4)
Supercentenarians.DESeq2 <- FindAllMarkers(Supercentenarians, test.use = "DESeq2", only.pos = TRUE, return.thresh = 0.001, logfc.threshold = 0.4)

saveRDS(Supercentenarians.wilcox, file = "Supercentenarians_wilcox.rds")
saveRDS(Supercentenarians.bimod, file = "Supercentenarians_bimod.rds")
saveRDS(Supercentenarians.roc, file = "Supercentenarians_roc.rds")
saveRDS(Supercentenarians.t, file = "Supercentenarians_t.rds")
saveRDS(Supercentenarians.negbinom, file = "Supercentenarians_negbinom.rds")
saveRDS(Supercentenarians.poisson, file = "Supercentenarians_poisson.rds")
saveRDS(Supercentenarians.LR, file = "Supercentenarians_LR.rds")
saveRDS(Supercentenarians.MAST, file = "Supercentenarians_MAST.rds")
saveRDS(Supercentenarians.DESeq2, file = "Supercentenarians_DESeq2.rds")

saveRDS(clu_ann, file = "clu_ann.rds")
saveRDS(clu_marker, file = "clu_marker.rds")

BiocManager::install("DESeq2")
#here should be an assessnode dont know whats not working needs to fix
node.scores <- AssessNodes(object = Supercentenarians)#, genes.training = VariableFeatures(Supercentenarians)
node.scores[order(node.scores$oobe,decreasing = TRUE),] -> node.scores
nodes.merge <- node.scores[which(node.scores$oobe > 0.05),]
nodes.to.merge <- sort(nodes.merge$node)
Supercentenarians.merged <- Supercentenarians

#Merge nodes
Supercentenarians <- RenameIdents(Supercentenarians, "3" = "2", "7" = "6")
Supercentenarians <- BuildClusterTree(Supercentenarians, reorder.numeric = TRUE, reorder = TRUE, dims = 1:19)
PlotClusterTree(object = Supercentenarians)