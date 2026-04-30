#Librarys anounsment
library(scater)
library(Seurat)
library(cowplot)
library(dplyr)
library(patchwork)
library(dbplyr)
library(EnhancedVolcano)
library(base)
library(ggplot2)
library(ggrepel)
library(SeuratWrappers)
library(clustifyr)
library(clustifyrdata)
library(clustree)

#Setting directory
setwd("G:/My Drive/SenescenceProject/MyTMSAnalysis")
setwd("D:/Anna")
setwd("F:/Immunaging")

Senis.big <- ReadH5AD("tabula-muris-senis-droplet-official-raw-obj.h5ad")

mdata <- read.csv(file = "GSM4505404_tabula-muris-senis-droplet-official-raw-obj-metadata.csv", row.names = 1)

#Adding metadata
Senis.big <- AddMetaData(Senis.big, metadata = mdata)

saveRDS(Senis.big, file = "MetadataAddedData.rds")
Senis.big <- readRDS(file = "MetadataAddedData.rds")

#Creating Seurat Object
Senis.big <- CreateSeuratObject(counts = Senis.big@assays$RNA@counts, project = "TMS_R_Analysis", 
                                meta.data = Senis.big@meta.data, min.cells = 3, min.features = 250)

Senis.big[["Age_group"]] <- plyr::mapvalues(x = Senis.big$age, from = c("1m", "3m", "18m", "21m", "24m", "30m"),
                                            to = c("Young", "Young", "Old", "Old", "Old", "supercentenarians"))

Senis.big@meta.data[["cell.ontology.class"]] <- NULL
Senis.big@meta.data[["cell.ontology.id"]] <- NULL
Senis.big@meta.data[["free.annotation"]] <- NULL
Senis.big@meta.data[["subtissue"]] <- NULL
Senis.big@meta.data[["tissue.free.annotation"]] <- NULL
Senis.big@meta.data[["cell_ontology_id"]] <- NULL

#Save rds file
saveRDS(Senis.big, file = "TMS+Metadata.rds")
Senis.big <- readRDS(file = "TMS+Metadata.rds")

#QC dont needed, the raw data after qc
#Senis.big[['percent.ribo']] <- PercentageFeatureSet(Senis.big, pattern = "^Rp[sl]")
#VlnPlot(Senis.big, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo"), ncol = 3)
#plot1 <- FeatureScatter(Senis.big, feature1 = "nCount_RNA", feature2 = "percent.ribo")
#plot2 <- FeatureScatter(Senis.big, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2

#Normalization
Senis.big <- NormalizeData(Senis.big, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(Senis.big, file = "TMS+Normalized.rds")
Senis.big <- readRDS(file = "TMS+Normalized.rds")

#Set senescence subset
p16.stats <- as.data.frame(summary(FetchData(object = Senis.big, vars = "Cdkn2a")))
p21.stats <- as.data.frame(summary(FetchData(object = Senis.big, vars = "Cdkn1a")))
Senis.big$Senescence.group <- "Senescence.neg"
Senis.big$Senescence.group[WhichCells(Senis.big, expression = (Cdkn2a > 0 & Cdkn1a > 0))] <- "Senescence.pos"
Idents(Senis.big) <- "Senescence.group"
Senescence.list <- SplitObject(Senis.big, split.by = "Senescence.group")
saveRDS(Senescence.list, file = "Senescence_list.rds")
saveRDS(Senescence.list$Senescence.pos, file = "Senescence_subset.rds")

write.table(p21.stats, file = "p21_stats.txt")
write.table(p16.stats, file = "p16_stats.txt")


saveRDS(Senis.big, "TMS+senescence.rds")
Tissue.list <- SplitObject(Senis.big, split.by = "tissue")
saveRDS(Tissue.list, "Tissue_list.rds")
Senis.big <- readRDS(file = "TMS+senescence.rds")

#Highly variable expessed genes
Senis.big <- FindVariableFeatures(Senis.big, selection.method = "vst")#HVG by TMS
top10 <- head(VariableFeatures(Senis.big),15)
plot1 <- VariableFeaturePlot(Senis.big)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Senis.big), file = "TMS_Highly_Variable_Expessed_Genes.csv")

saveRDS(Senis.big, file = "TMS+VariableGenes.rds")
Senis.big <- readRDS(file = "TMS+VariableGenes.rds")

Senis.big <- ScaleData(Senis.big)

saveRDS(Senis.big, file = "TMS+ScaledData.rds")
Senis.big <- readRDS(file = "TMS+ScaledData.rds")

Senis.big <- RunPCA(Senis.big, npcs = 150, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(Senis.big, reduction = "pca")

saveRDS(Senis.big, file = "TMS+PCA.rds")
Senis.big <- readRDS(file = "TMS+PCA.rds")

ElbowPlot(Senis.big, ndims = 150)

Senis.big <- FindNeighbors(Senis.big, dims = 1:150)

#louvian
Senis.big <- FindClusters(Senis.big, resolution = 0.8)

Senis.big <- RunTSNE(Senis.big, dims = 1:150, perplexity = 50)
DimPlot(Senis.big, reduction = "tsne")

Senis.big <- RunUMAP(Senis.big, dims = 1:150, perplexity = 50)
DimPlot(Senis.big, reduction = "umap")

clustree(Senis.big, prefix = "RNA_snn_res.")

Senis.big <- BuildClusterTree(Senis.big, reorder.numeric = TRUE, reorder = TRUE, dims = 1:150)
PlotClusterTree(object = Senis.big)

sig <- preprocess.signatures("T_cell_markers.csv")
r <- SCINA(Senis.big@assays$RNA@data, sig)

saveRDS(Senis.big.markers, "Senis_big_markers.rds")
saveRDS(Senis.big.age.cons.markers, "Senis_big_age_cons_markers.rds")
saveRDS(Senis.big.tissues.cons.markers, "Senis_big_tissues_cons_markers.rds")
saveRDS(Senis.big, file = "Clustred_TMS.rds")
Senis.big <- readRDS(file = "Clustred_TMS.rds")

Senis.big$t_or_macrophages <- "Cells"
Senis.big$t_or_macrophages <- plyr::mapvalues(x = Senis.big$cell.ontology.class, from = c("leukocyte", "naive T cell", 
                                                                             "hematopoietic precursor cell", "monocyte",
                                                                             "promonocyte", "NK cell","macrophage", 
                                                                             "T cell", "myeloid cell", "lymphocyte",
                                                                             "Kupffer cell", "myeloid leukocyte",
                                                                             "non-classical monocyte", 
                                                                             "alveolar macrophage", 
                                                                             "myeloid dendritic cell", 
                                                                             "CD8-positive, alpha-beta T cell",
                                                                             "classical monocyte", "mature NK T cell",
                                                                             "intermediate monocyte", 
                                                                             "CD4-positive, alpha-beta T cell",
                                                                             "regulatory T cell", "immature NKT cell",
                                                                             "macrophage dendritic cell progenitor",
                                                                             "DN3 thymocyte", "thymocyte", 
                                                                             "immature T cell", 
                                                                             "double negative T cell", 
                                                                             "DN4 thymocyte", 
                                                                             "professional antigen presenting cell",
                                                                             "lung macrophage", 
                                                                             "hematopoietic stem cell", "blood cell"), 
                                          to = c("t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", 
                                                 "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", 
                                                 "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", 
                                                 "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", 
                                                 "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages",
                                                 "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", 
                                                 "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", 
                                                 "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages"))
Senis.big$t_or_macrophages_2 <- plyr::mapvalues(x = Senis.big$free.annotation, from = c("granulocyte-monocyte progenitor", 
                                                                             "CD45", "monocyte", "T cell", 
                                                                             "leukocyte"), 
                                 to = c("t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages", "t_or_macrophages"))

t_or_macrophages_subset <- subset(Senis.big, subset = t_or_macrophages == "t_or_macrophages" | t_or_macrophages_2 == "t_or_macrophages")

#Highly variable expressed genes
t_or_macrophages_subset <- FindVariableFeatures(t_or_macrophages_subset, selection.method = "mean.var.plot", 
                                                dispersion.cutoff = c(0.5, 10), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(t_or_macrophages_subset),15)
plot1 <- VariableFeaturePlot(t_or_macrophages_subset)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(t_or_macrophages_subset), file = "t_or_macrophages_subset_Highly_Variable_Expessed_Genes.csv")

saveRDS(t_or_macrophages_subset, file = "t_or_macrophages_subsetVariableGenes.rds")
t_or_macrophages_subset <- readRDS(file = "t_or_macrophages_subsetVariableGenes.rds")

t_or_macrophages_subset <- ScaleData(t_or_macrophages_subset)

saveRDS(t_or_macrophages_subset, file = "t_or_macrophages_subsetScaledData.rds")
t_or_macrophages_subset <- readRDS(file = "t_or_macrophages_subsetScaledData.rds")

Idents(t_or_macrophages_subset) <- t_or_macrophages_subset$orig.ident
t_or_macrophages_subset <- RunPCA(t_or_macrophages_subset, npcs = 20, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(t_or_macrophages_subset, reduction = "pca")

saveRDS(t_or_macrophages_subset, file = "t_or_macrophages_subsetPCA.rds")
t_or_macrophages_subset <- readRDS(file = "t_or_macrophages_subsetPCA.rds")

ElbowPlot(t_or_macrophages_subset, ndims = 20)

t_or_macrophages_subset <- FindNeighbors(t_or_macrophages_subset, dims = 1:20)

res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
t_or_macrophages_subset <- FindClusters(t_or_macrophages_subset, resolution = res)

t_or_macrophages_subset <- RunTSNE(t_or_macrophages_subset, dims = 1:20)
DimPlot(t_or_macrophages_subset, reduction = "tsne")

t_or_macrophages_subset <- RunUMAP(t_or_macrophages_subset, dims = 1:20)
DimPlot(t_or_macrophages_subset, reduction = "umap")

clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.")
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Ptprc", node_colour_aggr = "median")#cd45 all wbc
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Cd3d", node_colour_aggr = "median")# t cells
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Cd3g", node_colour_aggr = "median")# t cells
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Cd3e", node_colour_aggr = "median")#t cells
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Cd4", node_colour_aggr = "median")#
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Cd34", node_colour_aggr = "median")#
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Cd8a", node_colour_aggr = "median")#cd8 t cells
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Cd8b1", node_colour_aggr = "median")#cd8 t cells
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Ncam1", node_colour_aggr = "median")#cd4 t cells
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Fcgr3", node_colour_aggr = "median")#phages, cd16
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Cd14", node_colour_aggr = "median")#phages
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Cd19", node_colour_aggr = "median")#b cells
clustree(t_or_macrophages_subset, prefix = "RNA_snn_res.", node_colour = "Hba-a1", node_colour_aggr = "median")#rbc

t_or_macrophages_subset_ <- clustify(t_or_macrophages_subset, ref_mouse.rnaseq, cluster_col = "RNA_snn_res.1.1", 
                                    obj_out = TRUE, threshold = 0.35)

t_cells <- subset(t_or_macrophages_subset_, subset = type == "T cells")

#Highly variable expressed genes
t_cells <- FindVariableFeatures(t_cells, selection.method = "mean.var.plot", 
                                                dispersion.cutoff = c(0.5, 10), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(t_cells),15)
plot1 <- VariableFeaturePlot(t_cells)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(t_cells), 
          file = "t_cells_Highly_Variable_Expessed_Genes_recluster.csv")

saveRDS(t_cells, file = "t_cells_VariableGenes_recluster.rds")
t_cells <- readRDS(file = "t_cells_VariableGenes_recluster.rds")

t_cells <- ScaleData(t_cells)

saveRDS(t_cells, file = "t_cells_ScaledData_recluster.rds")
t_cells <- readRDS(file = "t_cells_ScaledData_recluster.rds")

Idents(t_cells) <- t_cells$orig.ident
t_cells <- RunPCA(t_cells, npcs = 20, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(t_cells, reduction = "pca")

saveRDS(t_cells, file = "t_cells_PCA_recluster.rds")
t_cells <- readRDS(file = "t_cells_PCA_recluster.rds")

ElbowPlot(t_cells, ndims = 20)

t_cells <- FindNeighbors(t_cells, dims = 1:20)

res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
t_cells <- FindClusters(t_cells, resolution = res)

t_cells <- RunTSNE(t_cells, dims = 1:20)
DimPlot(t_cells, reduction = "tsne")

t_cells <- RunUMAP(t_cells, dims = 1:20)
DimPlot(t_cells, reduction = "umap")

clustree(t_cells, prefix = "RNA_snn_res.")

clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Ptprc", node_colour_aggr = "median")#cd45 all wbc
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Cd3d", node_colour_aggr = "median")# t cells
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Cd3g", node_colour_aggr = "median")# t cells
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Cd3e", node_colour_aggr = "median")#t cells
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Cd4", node_colour_aggr = "median")#
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Cd34", node_colour_aggr = "median")#
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Cd8a", node_colour_aggr = "median")#cd8 t cells
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Cd8b1", node_colour_aggr = "median")#cd8 t cells
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Ncam1", node_colour_aggr = "median")#cd4 t cells
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Fcgr3", node_colour_aggr = "median")#phages, cd16
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Cd14", node_colour_aggr = "median")#phages
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Cd19", node_colour_aggr = "median")#b cells
clustree(t_cells, prefix = "RNA_snn_res.", node_colour = "Hba-a1", node_colour_aggr = "median")#rbc

t_cells_ <- clustify(t_cells, ref_immgen, cluster_col = "RNA_snn_res.0.8", obj_out = TRUE, threshold = 0.5)
DimPlot(t_cells_, group.by = "type.clustify")

Idents(t_cells_) <- t_cells_$type.clustify

t_cells_ <- RenameIdents(t_cells_, "T cells (T.8MEM.OT1.D45.LISOVA)" = "CD8 T", "T cells (T.8EFF.OT1LISO)" = "CD8 T",
                         "T cells (T.DP.69-)" = "DP T", "T cells (T.4MEM44H62L)" ="CD4 T", "T cells (T.Tregs)" = "CD4 T",
                         "T cells (T.8MEMKLRG1-CD127+.D8.LISOVA)" = "CD8 T", "T cells (T.CD8.1H)" = "CD8 T", 
                         "T cells (T.CD4TESTCJ)" = "CD4 T", "T cells (T.4MEM49D+11A+.D30.LCMV)" = "CD4 T", 
                         "T cells (T.8SP24int)" = "CD8 T", "T cells (T.8NVE.OT1)" = "CD8 T", "T cells (T.DP69+)" = "DP T", 
                         "T cells (T.DN3A)" = "DN T", "T cells (T.DP)" = "DP T", "Tgd (Tgd.VG2+)" = "CD8 T") 

t_cells_$my_clusters <- t_cells_@active.ident
Idents(t_or_macrophages_subset) <- t_or_macrophages_subset$RNA_snn_res.0.1

saveRDS(t_or_macrophages_subset, file = "t_cells_recluster.rds")
t_or_macrophages_subset <- readRDS(file = "t_cells_recluster.rds")

t_or_macrophages_subset <- BuildClusterTree(object = t_or_macrophages_subset, dims = 1:20, reorder = TRUE, 
                                            reorder.numeric = TRUE)
PlotClusterTree(object = t_or_macrophages_subset)

t_or_macrophages_subset.markers <- list()

for (n in 1:17) {
  t_or_macrophages_subset.markers[[n]] <- FindMarkers(t_or_macrophages_subset, ident.1 = n, only.pos = TRUE)
}

t_or_macrophages_subset.age.cons.markers <- list()

for (n in 1:17) {
  t_or_macrophages_subset.age.cons.markers[[n]] <- FindConservedMarkers(t_or_macrophages_subset, ident.1 = n, 
                                                                        grouping.var = "Age_group")
}

saveRDS(t_or_macrophages_subset.markers, file = "t_or_macrophages_subset_markers.rds")
saveRDS(t_or_macrophages_subset.age.cons.markers, file = "t_or_macrophages_subset_age_cons_markers.rds")
write.csv(t_or_macrophages_subset.markers[[1]], "t_or_macrophages_subset_markers_1.csv")
write.csv(t_or_macrophages_subset.markers[[2]], "t_or_macrophages_subset_markers_2.csv")
write.csv(t_or_macrophages_subset.markers[[3]], "t_or_macrophages_subset_markers_3.csv")
write.csv(t_or_macrophages_subset.markers[[4]], "t_or_macrophages_subset_markers_4.csv")
write.csv(t_or_macrophages_subset.markers[[5]], "t_or_macrophages_subset_markers_5.csv")
write.csv(t_or_macrophages_subset.markers[[6]], "t_or_macrophages_subset_markers_6.csv")
write.csv(t_or_macrophages_subset.markers[[7]], "t_or_macrophages_subset_markers_7.csv")
write.csv(t_or_macrophages_subset.markers[[8]], "t_or_macrophages_subset_markers_8.csv")
write.csv(t_or_macrophages_subset.markers[[9]], "t_or_macrophages_subset_markers_9.csv")
write.csv(t_or_macrophages_subset.markers[[10]], "t_or_macrophages_subset_markers_10.csv")
write.csv(t_or_macrophages_subset.markers[[11]], "t_or_macrophages_subset_markers_11.csv")
write.csv(t_or_macrophages_subset.markers[[12]], "t_or_macrophages_subset_markers_12.csv")
write.csv(t_or_macrophages_subset.markers[[13]], "t_or_macrophages_subset_markers_13.csv")
write.csv(t_or_macrophages_subset.markers[[14]], "t_or_macrophages_subset_markers_14.csv")
write.csv(t_or_macrophages_subset.markers[[15]], "t_or_macrophages_subset_markers_15.csv")
write.csv(t_or_macrophages_subset.markers[[16]], "t_or_macrophages_subset_markers_16.csv")
write.csv(t_or_macrophages_subset.markers[[17]], "t_or_macrophages_subset_markers_17.csv")

t_or_macrophages_subset <- RenameIdents(t_or_macrophages_subset, "25" = "19", "24" = "19", "23" = "19", "22" = "19", 
                                        "21" = "19", "20" = "19")
t_or_macrophages_subset <- BuildClusterTree(object = t_or_macrophages_subset, dims = 1:20, reorder = TRUE, 
                                            reorder.numeric = TRUE)
PlotClusterTree(object = t_or_macrophages_subset)

t_or_macrophages_subset$my_clusters <- t_or_macrophages_subset@active.ident

t_or_macrophages_subset_ <- clustify(t_or_macrophages_subset, ref_mouse.rnaseq, cluster_col = "my_clusters", obj_out = TRUE)

t_or_macrophages_subset_subset <- subset(t_or_macrophages_subset_, subset = type == "Monocytes (MO.6C+II-)" | 
                                           type == "T cells (T.Tregs)" | type == "T cells (T.DP.69-)"
                                         | type == "Monocytes (MO)"| type == "Macrophages (MF.RP)" 
                                         | type == "Macrophages (MF.480INT.NAIVE)"| 
                                           type == "Macrophages (MF.480HI.NAIVE)" | type == "Macrophages (MFAR-)"| 
                                           type == "Tgd (Tgd.imm.VG1+)" | type == "T cells (T.DN3A)" 
                                         | type == "unassigned") 
t_or_macrophages_subset_subset <- BuildClusterTree(object = t_or_macrophages_subset_subset, dims = 1:20, reorder = TRUE, 
                                            reorder.numeric = TRUE)
PlotClusterTree(object = t_or_macrophages_subset_subset)

t_or_macrophages_subset_subset.markers <- list()

for (n in 1:18) {
  t_or_macrophages_subset_subset.markers[[n]] <- FindMarkers(t_or_macrophages_subset_subset, ident.1 = n, only.pos = TRUE)
}

write.csv(t_or_macrophages_subset_subset.markers[[1]], "t_or_macrophages_subset_markers_1_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[2]], "t_or_macrophages_subset_markers_2_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[3]], "t_or_macrophages_subset_markers_3_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[4]], "t_or_macrophages_subset_markers_4_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[5]], "t_or_macrophages_subset_markers_5_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[6]], "t_or_macrophages_subset_markers_6_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[7]], "t_or_macrophages_subset_markers_7_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[8]], "t_or_macrophages_subset_markers_8_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[9]], "t_or_macrophages_subset_markers_9_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[10]], "t_or_macrophages_subset_markers_10_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[11]], "t_or_macrophages_subset_markers_11_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[12]], "t_or_macrophages_subset_markers_12_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[13]], "t_or_macrophages_subset_markers_13_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[14]], "t_or_macrophages_subset_markers_14_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[15]], "t_or_macrophages_subset_markers_15_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[16]], "t_or_macrophages_subset_markers_16_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[17]], "t_or_macrophages_subset_markers_17_3.csv")
write.csv(t_or_macrophages_subset_subset.markers[[18]], "t_or_macrophages_subset_markers_18_3.csv")



common_genes <- intersect(rownames(ref_tabula_muris_drop), rownames(ref_tabula_muris_facs))
reff <- merge(ref_tabula_muris_drop[common_genes,], ref_tabula_muris_facs[common_genes,], by = 0)
reff2 <- reff[,-1]
rownames(reff2) <- reff[,1]
reff <- reff2

Senis.big_clusters <- clustify(Senis.big, reff, cluster_col = "RNA_snn_res.3", obj_out = TRUE)
