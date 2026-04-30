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

setwd("D:/Anna/GSE152522")

library(MazamaCoreUtils)
library(pbapply)
library(prob)
source(file = "D:/Anna/HezisData/AssessNodesAdaptedfromSeurat.R")

raw_counts <- read.csv(file = "cd4t0n6_umi.csv", sep = "\t", row.names = 1)
object.size(raw_counts)
dim(raw_counts)

cd4t0n6 <- CreateSeuratObject(counts = raw_counts, project = "cd4t0n6", min.cells = 3, min.features = 200)
object.size(cd4t0n6)

saveRDS(cd4t0n6, file = "AgingDataSeuratObject.rds")
cd4t0n6 <- readRDS("AgingDataSeuratObject.rds")

Idents(object = cd4t0n6) <- "CD4_T_Cells"
cd4t0n6[["percent.mt"]] <- PercentageFeatureSet(cd4t0n6, pattern = "^MT-")
cd4t0n6[['percent.ribo']] <- PercentageFeatureSet(cd4t0n6, pattern = "^RP[SL]")
p1 <- VlnPlot(cd4t0n6, features = "nFeature_RNA", pt.size = 0) + labs(title = "nGenes", tag = "A") + NoLegend()
p2 <- VlnPlot(cd4t0n6, features = "nCount_RNA", pt.size = 0) + labs(title = "nUMI", tag = "B") + NoLegend()
p3 <- VlnPlot(cd4t0n6, features = "percent.mt", pt.size = 0) + labs(title = "Mito percrnt", tag = "C") + NoLegend()
p1 <- AugmentPlot(plot = p1) + labs(y = "nGenes")
p2 <- AugmentPlot(plot = p2) + labs(y = "nUMI")
p3 <- AugmentPlot(plot = p3) + labs(y = "Mito percrnt")
p1 + p2 + p3
VlnPlot(cd4t0n6, features = "percent.ribo")
VlnPlot(cd4t0n6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) + ggtitle(label = c("nGenes", "nUMI", "Mito percrnt"))

#Metadata.aging <- read.csv(file = "Metadata.csv")

#Mouse <- as.array(Metadata.aging$Mouse)
#Batch <- as.array(Metadata.aging$Batch)
#Age_group <- as.array(Metadata.aging$Age_group)
#Subset <- as.array(Metadata.aging$Subset)

#agingdata <- AddMetaData(agingdata, metadata = Mouse, col.name = "Mouse")
#agingdata <- AddMetaData(agingdata, metadata = Batch, col.name = "Batch")
#agingdata <- AddMetaData(agingdata, metadata = Age_group, col.name = "Age_group")
#agingdata <- AddMetaData(agingdata, metadata = Subset, col.name = "Subset")

saveRDS(cd4t0n6, file = "AgingDataMetadataAdded.rds")
cd4t0n6 <- readRDS("AgingDataMetadataAdded.rds")

plot1 <- FeatureScatter(cd4t0n6, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cd4t0n6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

cd4t0n6 <- subset(cd4t0n6, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

cd4t0n6 <- NormalizeData(cd4t0n6, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(cd4t0n6, file = "AgingDataNormalized.rds")
cd4t0n6 <- readRDS(file = "AgingDataNormalized.rds")

cd4t0n6 <- FindVariableFeatures(agingdata, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 2), mean.function = ExpMean, dispersion.function = LogVMR)

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

agingdata <- RunPCA(agingdata, features = VariableFeatures(object = agingdata))

DimPlot(agingdata, reduction = "pca")

agingdata <- JackStraw(agingdata, num.replicate = 240, dims = 30)
agingdata <- ScoreJackStraw(agingdata, dims = 1:30)
JackStrawPlot(agingdata, dims = 1:30)

saveRDS(agingdata, file = "AgingDataJackStraw.rds")
agingdata <- readRDS(file = "AgingDataJackStraw.rds")

ElbowPlot(agingdata, ndims = 30)# + geom_point(data = agingdata@reductions$pca@stdev, mapping = aes(colour = 0.001))+geom_hline(aes(yintercept = 0.001))

agingdata <- RunTSNE(object = agingdata, dims = 1:20, perplexity = 30)
TSNEPlot(object = agingdata, label = TRUE, pt.size = 0.5)

saveRDS(agingdata, file = "AgingData&tSNE.rds")
agingdata <- readRDS("AgingData&tSNE.rds")