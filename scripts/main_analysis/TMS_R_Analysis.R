#Librarys anounsment
library(reticulate)
library(scater)
library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)
library(dbplyr)
library(data.table)
library(EnhancedVolcano)
library(hdf5r)
library(SeuratDisk)
library(base)
library(renv)
library(ggplot2)
library(ggrepel)
library(SeuratWrappers)
library(clustifyr)
library(clustifyrdata)

#Setting directory
setwd("G:/My Drive/SenescenceProject/MyTMSAnalysis")
setwd("D:/Anna")

#Load raw data
Senis.big <- ReadH5AD(file = "tabula-muris-senis-droplet-official-raw-obj.h5ad")

#Save rds file
saveRDS(Senis.big, file = "RawData.rds")
Senis.big <- readRDS(file = "RawData.rds")

#Loading metadata
mdata <- read.csv(file = "GSM4505404_tabula-muris-senis-droplet-official-raw-obj-metadata.csv")
rownames(mdata) <- mdata[,1]

#Adding metadata
Senis.big <- AddMetaData(Senis.big, metadata = mdata)

saveRDS(Senis.big, file = "MetadataAddedData.rds")
Senis.big <- readRDS(file = "MetadataAddedData.rds")

#Creating Seurat Object
Senis.big <- CreateSeuratObject(counts = Senis.big@assays$RNA@counts, project = "TMS_R_Analysis", 
                                meta.data = Senis.big@meta.data)

#Save rds file
saveRDS(Senis.big, file = "SeuratObject.rds")
Senis.big <- readRDS(file = "SeuratObject.rds")

Senis.big[["Age_group"]] <- plyr::mapvalues(x = Senis.big$age, from = c("1m", "3m", "18m", "21m", "24m", "30m"),
                                            to = c("Young", "Young", "Old", "Old", "Old", "Supercentenarian"))

Senis.big[["age_numbers"]] <- plyr::mapvalues(x = Senis.big$age, from = c("1m", "3m", "18m", "21m", "24m", "30m"),
                                        to = c("1", "3", "18", "21", "24", "30"))

Senis.big[['nFeatures_RNA']] <- NULL
Senis.big@meta.data[["cell.ontology.class"]] <- NULL
Senis.big@meta.data[["cell.ontology.id"]] <- NULL
Senis.big@meta.data[["free.annotation"]] <- NULL
Senis.big@meta.data[["tissue.free.annotation"]] <- NULL
Senis.big@meta.data[["cell_ontology_id"]] <- NULL

Senis.big[['percent_ribo']] <- PercentageFeatureSet(Senis.big, pattern = "^Rp[sl]")
VlnPlot(Senis.big, features = c("nFeature_RNA", "nCount_RNA", "percent_ribo"), ncol = 3, pt.size = 0.3)
p1 <- VlnPlot(Senis.big, features = "nFeature_RNA", pt.size = 0) + labs(title = "Number of genes detected in cell", tag = "A", x = NULL, y = NULL) + NoLegend()
p2 <- VlnPlot(Senis.big, features = "nCount_RNA", pt.size = 0) + labs(title = "Total number of molecules detected in cell", tag = "B") + NoLegend()
p3 <- VlnPlot(Senis.big, features = "percent_ribo", pt.size = 0) + labs(title = "Percent of ribosomal genes", tag = "D") + NoLegend()
p1 <- AugmentPlot(plot = p1) + labs(y = "nGenes")
p2 <- AugmentPlot(plot = p2) + labs(y = "nUMI")
p3 <- AugmentPlot(plot = p3) + labs(y = "Percent of ribosomal genes")
p1 + p2 + p3 
plot1 <- FeatureScatter(Senis.big, feature1 = "nCount_RNA", feature2 = "percent_ribo")
plot2 <- FeatureScatter(Senis.big, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + NoLegend() + plot2 + NoLegend()

#Normalization
Senis.big <- NormalizeData(Senis.big, normalization.method = "LogNormalize", scale.factor = 10000)
 
saveRDS(Senis.big, file = "Normalized.rds")
Senis.big <- readRDS(file = "Normalized.rds")

#Set senescence subset
p16.stats <- as.data.frame(summary(FetchData(object = Senis.big, vars = "Cdkn2a")))
p21.stats <- as.data.frame(summary(FetchData(object = Senis.big, vars = "Cdkn1a")))
Senis.big$Senescence.group <- "Senescence.neg"
Senis.big$Senescence.group[WhichCells(Senis.big, expression = Cdkn2a > 0 & Cdkn1a > 0)] <- "Senescence.pos"
Idents(Senis.big) <- "Senescence.group"
Senescence.list <- SplitObject(Senis.big, split.by = "Senescence.group")
saveRDS(Senescence.list, file = "Senescence_list.rds")
saveRDS(Senescence.list$Senescence.pos, file = "Senescence_subset.rds")

write.csv(p21.stats, file = "p21_stats.csv")
write.csv(p16.stats, file = "p16_stats.csv")

Senis.big$p16.group <- "p16.neg"
Senis.big$p16.group[WhichCells(Senis.big, expression = Cdkn2a > 0)] <- "p16.pos"
Idents(Senis.big) <- "p16.group"
p16.list <- SplitObject(Senis.big, split.by = "p16.group")

saveRDS(p16.list, file = "p16_list.rds")

saveRDS(p16.list$p16.pos, file = "p16_subset.rds")

Senis.big$p21.group <- "p21.neg"
Senis.big$p21.group[WhichCells(Senis.big, expression = Cdkn1a > 0)] <- "p21.pos"
Idents(Senis.big) <- "p21.group"
p21.list <- SplitObject(Senis.big, split.by = "p21.group")

saveRDS(p21.list, file = "p21_list.rds")

saveRDS(p21.list$p21.pos, file = "p21_subset.rds")

#Splitting the data by tissues,,,,continue after
Tissue.list <- SplitObject(Senis.big, split.by = "tissue")
saveRDS(Tissue.list, "Tissue_list.rds")

#Highly variable expressed genes
Senis.big <- FindVariableFeatures(Senis.big, selection.method = "mean.var.plot", dispersion.cutoff = c(0.5, Inf), 
                                  mean.cutoff = c(0.0125, 3))
top10 <- head(VariableFeatures(Senis.big),15)
plot1 <- VariableFeaturePlot(Senis.big)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Senis.big), file = "Highly_Variable_Expessed_Genes.csv")

saveRDS(Senis.big, file = "VariableGenes.rds")
Senis.big <- readRDS(file = "VariableGenes.rds")

Senis.big <- ScaleData(Senis.big, features = Senescence_genes_list$REACTOME_CELLULAR_SENESCENCE)

saveRDS(Senis.big, file = "ScaledData.rds")
Senis.big <- readRDS(file = "ScaledData.rds")

Senis.big <- RunPCA(Senis.big, features = Senescence_genes_list$REACTOME_CELLULAR_SENESCENCE
                    , ndims.print = 1:5, nfeatures.print = 5)
DimPlot(Senis.big, reduction = "pca")

VizDimLoadings(Senis.big, dims = c(1:3, 70:75), reduction = "pca")
DimHeatmap(Senis.big, dims = c(1:3, 70:75), cells = 50000, balanced = TRUE)

saveRDS(Senis.big, file = "PCA.rds")
Senis.big <- readRDS(file = "PCA.rds")

#
Senis.big <- JackStraw(Senis.big, num.replicate = 50, dims = 20)
Senis.big <- ScoreJackStraw(Senis.big, dims = 1:20)
JackStrawPlot(Senis.big, dims = 1:15)

saveRDS(Senis.big, file = "JackStraw.rds")
Senis.big <- readRDS(file = "JackStraw.rds")

ElbowPlot(Senis.big, ndims = 15)

Senis.big <- FindNeighbors(Senis.big, dims = 1:15, nn.eps = 0.5)

#louvian
Senis.big <- FindClusters(Senis.big, resolution = seq(0,1.2,by = 0.1))

Senis.big <- RunTSNE(Senis.big, dims = 1:15, periplexity = 50)
DimPlot(Senis.big, reduction = "tsne")

Senis.big <- RunUMAP(Senis.big, min.dist = 0.3, dims = 1:15)
DimPlot(Senis.big, reduction = "umap")

Senis.big <- BuildClusterTree(Senis.big, reorder.numeric = TRUE, reorder = TRUE, dims = 1:15)
PlotClusterTree(object = Senis.big)

saveRDS(Senis.big, file = "Clustred_TMS.rds")

Senis.big <- clustify(input = Senis.big, cluster_col = "RNA_snn_res.3", ref_mat = ref_tabula_muris_drop, seurat_output = TRUE)
Idents(Senis.big) <- Senis.big$type
Senis.big.markers <- FindAllMarkers(Senis.big, logfc.threshold = 0.4, only.pos = TRUE)
Senis.big.markers %>% group_by(cluster) %>% top_n(n = 5, wt = p_val_adj)
Senis.big.list <- clustify_lists(input = Senis.big, marker = Senis.big.markers, cluster_col = "classified", metric = "jaccard")

Senescence.markers.pos <- FindMarkers(Senis.big, group.by = "Senescence.group", ident.1 = "Senescence.pos", ident.2 = "Senescence.neg", logfc.threshold = 0.4, only.pos = TRUE)
Senescence.markers.neg <- FindMarkers(Senis.big, group.by = "Senescence.group", ident.1 = "Senescence.neg", ident.2 = "Senescence.pos", logfc.threshold = 0.4, only.pos = TRUE)
Senescence.markers.cons <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "Senescence.pos")
age.markers.1m <- FindConservedMarkers(Senis.big, grouping.var = "Senescence.group", ident.1 = "1m")
age.markers.3m <- FindConservedMarkers(Senis.big, grouping.var = "Senescence.group", ident.1 = "3m")
age.markers.18m <- FindConservedMarkers(Senis.big, grouping.var = "Senescence.group", ident.1 = "18m")
age.markers.21m <- FindConservedMarkers(Senis.big, grouping.var = "Senescence.group", ident.1 = "21m")
age.markers.24m <- FindConservedMarkers(Senis.big, grouping.var = "Senescence.group", ident.1 = "24m")
age.markers.30m <- FindConservedMarkers(Senis.big, grouping.var = "Senescence.group", ident.1 = "30m")

Senescence.markers.pos[order(Senescence.markers.pos$p_val_adj),]
top20.pos <- head(rownames(Senescence.markers.pos), 20)
top20.neg <- head(rownames(Senescence.markers.pos), 20)
top20.sen <- top20.pos + top20.neg
DoHeatmap(Senis.big, features = top20.pos)

Tissue.list <- SplitObject(Senis.big, split.by = "tissue")
Age.list <- SplitObject(Senis.big, split.by = "age")
Sex.list <- SplitObject(Senis.big, split.by = "sex")
Mouse.list <- SplitObject(Senis.big, split.by = "mouse.id")

saveRDS(Tissue.list$Tongue, file = "Tongue.rds")
saveRDS(Tissue.list$Heart_and_Aorta, file = "Heart_and_Aorta.rds")
saveRDS(Tissue.list$Marrow, file = "Marrow.rds")
saveRDS(Tissue.list$Mammary_Gland, file = "Mammary_Gland.rds")
saveRDS(Tissue.list$Fat, file = "Fat.rds")
saveRDS(Tissue.list$Kidney, file = "Kidney.rds")
saveRDS(Tissue.list$Liver, file = "Liver.rds")
saveRDS(Tissue.list$Lung, file = "Lung.rds")
saveRDS(Tissue.list$Limb_Muscle, file = "Limb_Muscle.rds")
saveRDS(Tissue.list$Pancreas, file = "Pancreas.rds")
saveRDS(Tissue.list$Spleen, file = "Spleen.rds")
saveRDS(Tissue.list$Thymus, file = "Thymus.rds")
saveRDS(Tissue.list$Bladder, file = "Bladder.rds")
saveRDS(Tissue.list$Skin, file = "Skin.rds")
saveRDS(Tissue.list$Large_Intestine, file = "Large_Intestine.rds")
saveRDS(Tissue.list$Trachea, file = "Trachea.rds")
saveRDS(Age.list$`1m`, file = "1m.rds")
saveRDS(Age.list$`3m`, file = "3m.rds")
saveRDS(Age.list$`18m`, file = "18m.rds")
saveRDS(Age.list$`21m`, file = "21m.rds")
saveRDS(Age.list$`24m`, file = "24m.rds")
saveRDS(Age.list$`30m`, file = "30m.rds")
saveRDS(Sex.list$male, file = "Male.rds")
saveRDS(Sex.list$female, file = "Female.rds")
saveRDS(Mouse.list$`24-M-60`, file = "Mouse_24_m_60.rds")
saveRDS(Mouse.list$`18-F-50`, file = "Mouse_18_f_50.rds")
saveRDS(Mouse.list$`18-F-51`, file = "Mouse_18_f_51.rds")
saveRDS(Mouse.list$`18-M-52`, file = "Mouse_18_m_52.rds")
saveRDS(Mouse.list$`18-M-53`, file = "Mouse_18_m_53.rds")
saveRDS(Mouse.list$`21-F-54`, file = "Mouse_21_f_54.rds")
saveRDS(Mouse.list$`21-F-55`, file = "Mouse_21_f_55.rds")
saveRDS(Mouse.list$`24-M-58`, file = "Mouse_24_m_58.rds")
saveRDS(Mouse.list$`24-M-59`, file = "Mouse_24_m_59.rds")
saveRDS(Mouse.list$`24-M-61`, file = "Mouse_24_m_61.rds")
saveRDS(Mouse.list$`1-M-63`, file = "Mouse_1_m_63.rds")
saveRDS(Mouse.list$`30-M-2`, file = "Mouse_30_m_2.rds")
saveRDS(Mouse.list$`30-M-3`, file = "Mouse_30_m_3.rds")
saveRDS(Mouse.list$`30-M-4`, file = "Mouse_30_m_4.rds")
saveRDS(Mouse.list$`30-M-5`, file = "Mouse_30_m_5.rds")
saveRDS(Mouse.list$`1-M-62`, file = "Mouse_1_m_62.rds")
saveRDS(Mouse.list$`3-M-8`, file = "Mouse_3_m_8.rds")
saveRDS(Mouse.list$`3-M-9`, file = "Mouse_3_m_9.rds")
saveRDS(Mouse.list$`3-M-8/9`, file = "Mouse_3_m_8_9.rds")
saveRDS(Mouse.list$`3-F-56`, file = "Mouse_3_f_56.rds")
saveRDS(Mouse.list$`3-F-57`, file = "Mouse_3_f_57.rds")
saveRDS(Mouse.list$`3-M-5/6`, file = "Mouse_3_m_5_6.rds")
saveRDS(Mouse.list$`3-M-7/8`, file = "Mouse_3_m_7_8.rds")

CD4_T_cells <- subset(Senis.big, subset = Cd4 > 0 & Cd3 > 0 & Cd45 > 0)
CTL_CD4_T_cells <- subset(CD4_T_cells, subset = Gzmk > 0)

DotPlot(Senis.big, features = c("Il27", "Il27ra", "Cdkn2a", "Cdkn1a", "Cdkn2b", "Cd4", "Gzmk"), group.by = "age")
DotPlot(Senis.big, features = c("Il27", "Il27ra", "Cdkn2a", "Cdkn1a", "Cdkn2b", "Cd4", "Gzmk"), group.by = "tissue")
DotPlot(Senis.big, features = c("Il27", "Il27ra", "Cdkn2a", "Cdkn1a", "Cdkn2b", "Cd4", "Gzmk"), group.by = "age", split.by = "sex")
DotPlot(Senis.big, features = c("Il27", "Il27ra", "Cdkn2a", "Cdkn1a", "Cdkn2b", "Cd4", "Gzmk"), group.by = "tissue", split.by = "sex")
VlnPlot(Senis.big, features = "Gzmk", group.by = "age")
VlnPlot(Senis.big, features = "Cdkn2a", group.by = "age", slot = "counts")
DotPlot(Senis.big, features = c("Il27", "Il27ra", "Cdkn2a", "Cdkn1a", "Cdkn2b", "Cd4", "Gzmk"), group.by = "age")
DotPlot(Senis.big, features = c("Il27", "Il27ra", "Cdkn2a", "Cdkn1a", "Cdkn2b", "Cd4", "Gzmk"), group.by = "tissue")
DotPlot(Senis.big, features = c("Il27", "Il27ra", "Cdkn2a", "Cdkn1a", "Cdkn2b", "Cd4", "Gzmk"), group.by = "sex")

m1.markers <- FindMarkers(Senis.big, group.by = "age", ident.1 = "1m", test.use = "roc", features = VariableFeatures(Senis.big))
m3.markers <- FindMarkers(Senis.big, group.by = "age", ident.1 = "3m", test.use = "roc", features = VariableFeatures(Senis.big))
m18.markers <- FindMarkers(Senis.big, group.by = "age", ident.1 = "18m", test.use = "roc", features = VariableFeatures(Senis.big))
m21.markers <- FindMarkers(Senis.big, group.by = "age", ident.1 = "21m", test.use = "roc", features = VariableFeatures(Senis.big))
m24.markers <- FindMarkers(Senis.big, group.by = "age", ident.1 = "24m", test.use = "roc", features = VariableFeatures(Senis.big))
m30.markers <- FindMarkers(Senis.big, group.by = "age", ident.1 = "30m", test.use = "roc", features = VariableFeatures(Senis.big))
conserved.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.1m_3m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.1m_18m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.1m_21m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.1m_24m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.1m_30m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.3m_18m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.3m_21m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.3m_24m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.3m_30m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.18m_21M.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.18m_24m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.18m_30m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.21m_24m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.21m_30m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")
conserved.24m_30m.markers <- FindConservedMarkers(Senis.big, grouping.var = "age", ident.1 = "1m")


saveRDS(m1.markers, file = "m1_markers.rds")
saveRDS(m3.markers, file = "m3_markers.rds")
saveRDS(m18.markers, file = "m18_markers.rds")
saveRDS(m21.markers, file = "m21_markers.rds")
saveRDS(m24.markers, file = "m24_markers.rds")
saveRDS(m30.markers, file = "m30_markers.rds")

DoHeatmap(Senis.big, features = VariableFeatures(Senis.big), group.by = "age")

#scater
sce <- as.SingleCellExperiment()
