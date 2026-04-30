library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
library(monocle3)
library(clustree)
library(SeuratWrappers)
library(base)
library(scater)
library(cowplot)
library(patchwork)
library(dbplyr)
library(data.table)
library(EnhancedVolcano)
library(SeuratDisk)
library(ggrepel)
library(clustifyr)
library(clustifyrdata)
library(clusterProfiler)

setwd("E:/Immunaging")
setwd("F:/Immunaging")
setwd("D:/Immunaging")
memory.limit(1024000)#memory in MB

Tongue <- readRDS(file = "Clustred_Senescence_Annotated_Tongue.rds")
Heart_and_Aorta <- readRDS(file = "Clustred_Senescence_Annotated_Heart_and_Aorta.rds")
Marrow <- readRDS(file = "Clustred_Senescence_Annotated_Marrow.rds")
Mammary_Gland <- readRDS(file = "Clustred_Senescence_Annotated_Mammary_Gland.rds")
Fat <- readRDS(file = "Clustred_Senescence_Annotated_Fat.rds")
Kidney <- readRDS(file = "Clustred_Senescence_Annotated_Kidney.rds")
Liver <- readRDS(file = "Clustred_Senescence_Annotated_Liver.rds")
Lung <- readRDS(file = "Clustred_Senescence_Annotated_Lung.rds")
Limb_Muscle <- readRDS(file = "Clustred_Senescence_Annotated_Limb_Muscle.rds")
Pancreas <- readRDS(file = "Clustred_Senescence_Annotated_Pancreas.rds")
Spleen <- readRDS(file = "Clustred_Senescence_Annotated_Spleen.rds")
Thymus <- readRDS(file = "Clustred_Senescence_Annotated_Thymus.rds")
Bladder <- readRDS(file = "Clustred_Senescence_Annotated_Bladder.rds")
Skin <- readRDS(file = "Clustred_Senescence_Annotated_Skin.rds")
Large_Intestine <- readRDS(file = "Clustred_Senescence_Annotated_Large_Intestine.rds")
Trachea <- readRDS(file = "Clustred_Senescence_Annotated_Trachea.rds")

Senis.big <- merge(Tongue, c(Heart_and_Aorta, Marrow, Mammary_Gland, Fat, Kidney, Liver, Lung, Limb_Muscle, Pancreas,
                             Spleen, Thymus, Bladder, Skin, Large_Intestine, Trachea))

saveRDS(Senis.big, file = "TMS_merge.rds")
Senis.big <- readRDS(file = "TMS_merge.rds")

cell_type <- Senis.big@active.ident
Senis.big$cell_type <- Senis.big@active.ident
cell_type <- gsub("^.*.monocyte.*", "monocyte", cell_type)
cell_type <- gsub("^.*.macrophage.*", "monocyte", cell_type)
cell_type <- gsub("-.*", "", cell_type)
Senis.big$cell_type_2 <- cell_type
Idents(Senis.big) <- Senis.big$cell_type_2

saveRDS(Senis.big, file = "TMS_merge_cell_type.rds")
Senis.big <- readRDS(file = "TMS_merge_cell_type.rds")

#Highly variable expressed genes
Senis.big <- FindVariableFeatures(Senis.big, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Senis.big),15)
plot1 <- VariableFeaturePlot(Senis.big)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Senis.big), file = "TMS_merge_Highly_Variable_Expessed_Genes.csv")

saveRDS(Senis.big, file = "TMS_merge_VariableGenes.rds")
Senis.big <- readRDS(file = "TMS_merge_VariableGenes.rds")

all.genes <- rownames(Senis.big)
Senis.big <- ScaleData(Senis.big, featuras = all.genes, vars.to.regress = c("tissue", "mouse.id", "age"))

saveRDS(Senis.big, file = "TMS_merge_ScaledData.rds")
Senis.big <- readRDS(file = "TMS_merge_ScaledData.rds")

Senis.big <- RunPCA(Senis.big, npcs = 109, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(Senis.big, reduction = "pca")

Senis.big <- JackStraw(Senis.big, num.replicate = 100, dims = 109)
Senis.big <- ScoreJackStraw(Senis.big, dims = 1:109)
JackStrawPlot(Senis.big, dims = 1:109)

saveRDS(Senis.big, file = "TMS_merge_PCA.rds")
Senis.big <- readRDS(file = "TMS_merge_PCA.rds")

ElbowPlot(Senis.big, ndims = 109)

Senis.big <- RunTSNE(Senis.big, dims = 1:30, perplexity = 50)
DimPlot(Senis.big, reduction = "tsne")

Senis.big <- RunUMAP(Senis.big, dims = 1:30, perplexity = 50)
DimPlot(Senis.big, reduction = "umap")

Senis.big[["Age_group_2"]] <- plyr::mapvalues(x = Senis.big$age, from = c("1m", "3m", "18m", "21m", "24m", "30m"),
                                            to = c("Young", "Young", "Old", "Old", "Old", "Old"))

write.csv(Senis.big@meta.data, "metedata_TMS.csv")

saveRDS(Senis.big, file = "TMS_merge_TSNE_UMAP.rds")
Senis.big <- readRDS(file = "TMS_merge_TSNE_UMAP.rds")

T_cells <- subset(Senis.big, subset = cell_type_2 == "T cell")
#T_cells <- subset(T_cells, subset = Cd3d > 0 | Cd3e > 0 | Cd3g > 0 | Cd247 > 0)
# using canonical T cell marker for positive selection
#T_cells <- subset(T_cells, subset = Ptprc > 0)

saveRDS(T_cells, file = "TMS_merge_T_cells.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells.rds")

#Highly variable expressed genes
T_cells <- FindVariableFeatures(T_cells, selection.method = "mean.var.plot",                                
                                dispersion.cutoff = c(0.8, Inf), mean.cutoff = c(0.02, Inf))#HVG by TMS
top10 <- head(VariableFeatures(T_cells),15)
plot1 <- VariableFeaturePlot(T_cells)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(T_cells), file = "TMS_merge_T_cells_Highly_Variable_Expessed_Genes.csv")

saveRDS(T_cells, file = "TMS_merge_T_cells_VariableGenes.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_VariableGenes.rds")

all.genes <- rownames(T_cells)
T_cells <- ScaleData(T_cells, featuras = all.genes, vars.to.regress = c("tissue", "age"))

saveRDS(T_cells, file = "TMS_merge_T_cells_ScaledData.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_ScaledData.rds")

T_cells <- RunPCA(T_cells, npcs = 65, ndims.print = 1:65, nfeatures.print = 5)
DimPlot(T_cells, reduction = "pca")

saveRDS(T_cells, file = "TMS_merge_T_cells_PCA.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_PCA.rds")

T_cells <- JackStraw(T_cells, num.replicate = 100, dims = 65)
T_cells <- ScoreJackStraw(T_cells, dims = 1:65)
JackStrawPlot(T_cells, dims = 1:65)

ElbowPlot(T_cells, ndims = 65)

T_cells <- RunTSNE(T_cells, dims = 1:65, perplexity = 50)
DimPlot(T_cells, reduction = "tsne")

T_cells <- RunUMAP(T_cells, dims = 1:65)
DimPlot(T_cells, reduction = "umap")

saveRDS(T_cells, file = "TMS_merge_T_cells_TSNE_UMAP.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_TSNE_UMAP.rds")

T_cells <- FindNeighbors(T_cells, dims = 1:65)

#louvian
T_cells <- FindClusters(T_cells, resolution = seq(0,3,by = 0.1))

clustree(T_cells, prefix = "RNA_snn_res.")
clustree(T_cells, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(T_cells) <- T_cells$RNA_snn_res.1.2

T_cells <- BuildClusterTree(T_cells, reorder.numeric = TRUE, reorder = TRUE, dims = 1:65)
PlotClusterTree(object = T_cells)

T_cells.markers <- FindAllMarkers(T_cells)
m <- T_cells.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(T_cells, features = m$gene)
write.csv(T_cells.markers, "TMS_T_cells_merkers.csv")

T_cells.clusters <- list()

for (n in 1:47) {
  T_cells.clusters[[n]] <- subset(T_cells, idents = n)
  saveRDS(T_cells.clusters[[n]], sprintf("TMS_merge_T_cells_clusters_%s.rds", n))
}

saveRDS(T_cells, file = "TMS_merge_T_cells_Clusters.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_Clusters.rds")

T_cells <- merge(Cluster_1, c(Cluster_2, Cluster_3, Cluster_4, Cluster_5, Cluster_6, Cluster_7, Cluster_8, Cluster_9,
                              Cluster_10, Cluster_11, Cluster_12, Cluster_13, Cluster_14, Cluster_15, Cluster_16,
                              Cluster_17, Cluster_18, Cluster_19, Cluster_20, Cluster_21, Cluster_22, Cluster_23, 
                              Cluster_24, Cluster_25, Cluster_26, Cluster_27, Cluster_28, Cluster_29, Cluster_30, 
                              Cluster_31, Cluster_32, Cluster_33, Cluster_34, Cluster_35, Cluster_36, Cluster_37, 
                              Cluster_38, Cluster_39, Cluster_40, Cluster_41, Cluster_42, Cluster_43, Cluster_44, 
                              Cluster_45, Cluster_46, Cluster_47))

saveRDS(T_cells, file = "TMS_merge_T_cells_Clusters_remerge.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_Clusters_remerge.rds")

T_cells$cell_type_T <- T_cells@active.ident

#Highly variable expressed genes
T_cells <- FindVariableFeatures(T_cells, selection.method = "mean.var.plot",                                
                                dispersion.cutoff = c(0.8, Inf), mean.cutoff = c(0.02, Inf))#HVG by TMS
top10 <- head(VariableFeatures(T_cells),15)
plot1 <- VariableFeaturePlot(T_cells)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(T_cells), file = "TMS_merge_T_cells_Highly_Variable_Expessed_Genes_2.csv")

saveRDS(T_cells, file = "TMS_merge_T_cells_VariableGenes_2.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_VariableGenes_2.rds")

all.genes <- rownames(T_cells)
T_cells <- ScaleData(T_cells, featuras = all.genes, vars.to.regress = c("tissue", "age"))

saveRDS(T_cells, file = "TMS_merge_T_cells_ScaledData_2.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_ScaledData_2.rds")

T_cells <- RunPCA(T_cells, npcs = 65, ndims.print = 1:65, nfeatures.print = 5)
DimPlot(T_cells, reduction = "pca")

saveRDS(T_cells, file = "TMS_merge_T_cells_PCA_2.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_PCA_2.rds")

T_cells <- RunTSNE(T_cells, dims = 1:65, perplexity = 50)
DimPlot(T_cells, reduction = "tsne")

T_cells <- RunUMAP(T_cells, dims = 1:65)
DimPlot(T_cells, reduction = "umap")

saveRDS(T_cells, file = "TMS_merge_T_cells_after_remerging.rds")
T_cells <- readRDS(file = "TMS_merge_T_cells_after_remerging.rds")

CD4_T <- subset(T_cells, subset = cell_type_T == "CD4 T")
CD8_T <- subset(T_cells, subset = cell_type_T == "CD8 T")

saveRDS(CD4_T, "total_CD4_T_TMS.rds")
saveRDS(CD8_T, "total_CD8_T_TMS.rds")
CD4_T <- readRDS("total_CD4_T_TMS.rds")
CD8_T <- readRDS("total_CD8_T_TMS.rds")

#Highly variable exposed genes
CD4_T <- FindVariableFeatures(CD4_T, selection.method = "mean.var.plot", dispersion.cutoff = c(0.55, Inf), 
                                mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(CD4_T),15)
plot1 <- VariableFeaturePlot(CD4_T)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(CD4_T), file = "TMS_merge_CD4_HVG.csv")

saveRDS(CD4_T, file = "TMS_merge_CD4_HVG.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_HVG.rds")

all.genes <- rownames(CD4_T)
CD4_T <- ScaleData(CD4_T, features = all.genes, vars.to.regress = "tissue")#

saveRDS(CD4_T, file = "TMS_merge_CD4_Scaled.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_Scaled.rds")

CD4_T <- RunPCA(CD4_T, npcs = 35, ndims.print = 1:35, nfeatures.print = 10)
DimPlot(CD4_T, reduction = "pca")

CD4_T <- JackStraw(CD4_T, num.replicate = 100, dims = 35)
CD4_T <- ScoreJackStraw(CD4_T, dims = 1:35)
JackStrawPlot(CD4_T, dims = 1:35)

saveRDS(CD4_T, file = "TMS_merge_CD4_T_pca.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_T_pca.rds")

ElbowPlot(CD4_T, ndims = 35)

CD4_T <- RunTSNE(CD4_T, dims = 1:35, perplexity = 30)
DimPlot(CD4_T, reduction = "tsne")

CD4_T <- RunUMAP(CD4_T, dims = 1:35)
DimPlot(CD4_T, reduction = "umap")

saveRDS(CD4_T, file = "TMS_merge_CD4_tsne_umap.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_tsne_umap.rds")

CD4_T <- FindNeighbors(CD4_T, dims = 1:35)

#louvian
CD4_T <- FindClusters(CD4_T, resolution = seq(0, 3, by = 0.1))

clustree(CD4_T, prefix = "RNA_snn_res.")
clustree(CD4_T, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(CD4_T, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "Foxp3", exprs = 'scale.data') + 
  scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(CD4_T) <- CD4_T$RNA_snn_res.2.3

CD4_T <- BuildClusterTree(CD4_T, reorder.numeric = TRUE, reorder = TRUE, dims = 1:35)
PlotClusterTree(object = CD4_T)

CD4_T.markers <- FindAllMarkers(CD4_T, only.pos = TRUE)
m <- CD4_T.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(CD4_T, features = m$gene)

write.csv(CD4_T.markers, "CD4_markers_tms_merge.csv")

CD4_T <- RenameIdents(CD4_T, "24" = "23", "7" = "6", "21" = "20", "9" = "8", "11" = "8", "12" = "8", 
                      "13" = "8", "14" = "8", "15" = "8", "16" = "8", "19" = "6", "5" = "6", "18" = "3", "17" = "3",
                      "4" = "3")
#19 - naive_Isg15, 15- rTregs, 17 - TEM

saveRDS(CD4_T, file = "TMS_merge_CD4_clustered.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_clustered.rds")

CD4_T <- BuildClusterTree(CD4_T, reorder.numeric = TRUE, reorder = TRUE, dims = 1:35)
PlotClusterTree(object = CD4_T)

CD4_T.marker_age <- list()

for (n in 1:9) {
  CD4_T.marker_age[[n]] <- FindConservedMarkers(CD4_T, ident.1 = n, grouping.var = "Age_group_2",
                                                min.cell.group = 0, min.cell.feature = 0)
  write.csv(CD4_T.marker_age[[n]], sprintf("CD4_T_marker_age_T_%s.csv", n))
}

CD4_T.markers <- FindAllMarkers(CD4_T, only.pos = TRUE)
m <- CD4_T.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(CD4_T, features = m$gene)

write.csv(CD4_T.markers, "CD4_markers_tms_merge_mergedclusters.csv")

CD4_T <- RenameIdents(CD4_T, "1" = "Not T", "2" = "Lymphoma", "3" = "Exhausted", "5" = "Naive", 
                      "4" = "aTregs", "6" = "Naive_Isg15", "7" = "rTregs", "8" = "Cytotoxic", "9" = "TEM")

CD4_T$subsets <- CD4_T@active.ident

CD4_T <- subset(CD4_T, subset = subsets != "Lymphoma" & subsets != "Not T")

#Highly variable exposed genes
CD4_T <- FindVariableFeatures(CD4_T, selection.method = "mean.var.plot", dispersion.cutoff = c(0.55, Inf), 
                              mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(CD4_T),15)
plot1 <- VariableFeaturePlot(CD4_T)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(CD4_T), file = "TMS_merge_CD4_HVG_without_Not_T.csv")

saveRDS(CD4_T, file = "TMS_merge_CD4_HVG_without_Not_T.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_HVG_without_Not_T.rds")

all.genes <- rownames(CD4_T)
CD4_T <- ScaleData(CD4_T, features = all.genes, vars.to.regress = "tissue")#

saveRDS(CD4_T, file = "TMS_merge_CD4_Scaled_without_Not_T.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_Scaled_without_Not_T.rds")

CD4_T <- RunPCA(CD4_T, npcs = 35, ndims.print = 1:35, nfeatures.print = 5)
DimPlot(CD4_T, reduction = "pca")

CD4_T <- JackStraw(CD4_T, num.replicate = 100, dims = 35)
CD4_T <- ScoreJackStraw(CD4_T, dims = 1:35)
JackStrawPlot(CD4_T, dims = 1:35)

saveRDS(CD4_T, file = "TMS_merge_CD4_T_pca_without_Not_T.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_T_pca_without_Not_T.rds")

ElbowPlot(CD4_T, ndims = 35)

CD4_T <- RunTSNE(CD4_T, dims = 1:35, perplexity = 30)
DimPlot(CD4_T, reduction = "tsne")

CD4_T <- RunUMAP(CD4_T, dims = 1:35)
DimPlot(CD4_T, reduction = "umap")

saveRDS(CD4_T, file = "TMS_merge_CD4_tsne_umap_without_Not_T.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_tsne_umap_without_Not_T.rds")

saveRDS(CD4_T, file = "TMS_merge_CD4_clustered_without_Not_T.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_clustered_without_Not_T.rds")



CD4_T_y_vs_o <- FindMarkers(CD4_T, ident.1 = "Young", group.by = "Age_group_2", test.use = "roc",
                            logfc.threshold = -Inf)
CD4_T_o_vs_y <- FindMarkers(CD4_T, ident.1 = "Old", group.by = "Age_group_2", test.use = "roc",
                            logfc.threshold = -Inf)
CD4_T_Naive <- FindConservedMarkers(CD4_T, ident.1 = "Naive", logfc.threshold = -Inf, grouping.var = "Age_group_2")
CD4_T_TEM <- FindConservedMarkers(CD4_T, ident.1 = "TEM", logfc.threshold = -Inf, grouping.var = "Age_group_2")
CD4_T_Cytotoxic <- FindConservedMarkers(CD4_T, ident.1 = "Cytotoxic", logfc.threshold = -Inf, 
                                        grouping.var = "Age_group_2")
CD4_T_Naive_Isg15 <- FindConservedMarkers(CD4_T, ident.1 = "Naive_Isg15", logfc.threshold = -Inf, 
                                          grouping.var = "Age_group_2")
CD4_T_Exhausted <- FindConservedMarkers(CD4_T, ident.1 = "Exhausted", logfc.threshold = -Inf, 
                                        grouping.var = "Age_group_2")
CD4_T_rTregs <- FindConservedMarkers(CD4_T, ident.1 = "rTregs", logfc.threshold = -Inf, grouping.var = "Age_group_2")
CD4_T_aTregs <- FindConservedMarkers(CD4_T, ident.1 = "aTregs", logfc.threshold = -Inf, grouping.var = "Age_group_2")

write.csv(CD4_T_y_vs_o, "CD4_T_y_vs_o.csv")
write.csv(CD4_T_o_vs_y, "CD4_T_o_vs_y.csv")
write.csv(CD4_T_Naive, "CD4_T_Naive.csv")
write.csv(CD4_T_TEM, "CD4_T_TEM.csv")
write.csv(CD4_T_Cytotoxic, "CD4_T_Cytotoxic.csv")
write.csv(CD4_T_Naive_Isg15, "CD4_T_Naive_Isg15.csv")
write.csv(CD4_T_Exhausted, "CD4_T_Exhausted.csv")
write.csv(CD4_T_rTregs, "CD4_T_rTregs.csv")
write.csv(CD4_T_aTregs, "CD4_T_aTregs.csv")

# add a column of NAs
CD4_T_Naive$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_T_Naive$diffexpressed[CD4_T_Naive$avg_logFC > 0.6 & CD4_T_Naive$p_val_adj < 0.5] <- "Naive"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_T_Naive$diffexpressed[CD4_T_Naive$avg_logFC < -1 & CD4_T_Naive$p_val_adj < 0.5] <- "Oters"
CD4_T_Naive$delabel <- NA
CD4_T_Naive$delabel[CD4_T_Naive$diffexpressed != "NO"] <- rownames(CD4_T_Naive)[CD4_T_Naive$diffexpressed != "NO"]
ggplot(CD4_T_Naive, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-1, 0.6), col="red") + geom_hline(
    yintercept = -log10(0.000005), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")

CD4_T_Cytotoxic_vs_TEM <- FindMarkers(CD4_T, ident.1 = "Cytotoxic", ident.2 = "TEM", logfc.threshold = -Inf)
CD4_T_Exhausted_vs_TEM <- FindMarkers(CD4_T, ident.1 = "Exhausted", ident.2 = "TEM", logfc.threshold = -Inf)
CD4_T_aTregs_vs_rTregs <- FindMarkers(CD4_T, ident.1 = "aTregs", ident.2 = "rTregs", logfc.threshold = -Inf)


write.csv(CD4_T_Cytotoxic_vs_TEM, "CD4_T_Cytotoxic_vs_TEM.csv")
write.csv(CD4_T_Exhausted_vs_TEM, "CD4_T_Exhausted_vs_TEM.csv")
write.csv(CD4_T_aTregs_vs_rTregs, "CD4_T_aTregs_vs_rTregs.csv")

TEM <- subset(CD4_T, subset = subsets == "TEM")
TEM_o_vs_y <- FindMarkers(TEM, ident.1 = "Old", ident.2 = "Young", logfc.threshold = -Inf, group.by = "Age_group_2")
write.csv(TEM_o_vs_y, "TEM_o_vs_y.csv")

# add a column of NAs
TEM_o_vs_y$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
TEM_o_vs_y$diffexpressed[TEM_o_vs_y$avg_logFC > 0.4 & TEM_o_vs_y$p_val_adj < 0.05] <- "Old"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
TEM_o_vs_y$diffexpressed[TEM_o_vs_y$avg_logFC < -0.4 & TEM_o_vs_y$p_val_adj < 0.05] <- "Young"
TEM_o_vs_y$delabel <- NA
TEM_o_vs_y$delabel[TEM_o_vs_y$diffexpressed != "NO"] <- rownames(TEM_o_vs_y)[TEM_o_vs_y$diffexpressed != "NO"]
ggplot(TEM_o_vs_y, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") + 
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")


# add a column of NAs
CD4_T_Cytotoxic_vs_TEM$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_T_Cytotoxic_vs_TEM$diffexpressed[CD4_T_Cytotoxic_vs_TEM$avg_logFC > 0.4 & 
                                       CD4_T_Cytotoxic_vs_TEM$p_val_adj < 0.05] <- "Cytotoxic"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_T_Cytotoxic_vs_TEM$diffexpressed[CD4_T_Cytotoxic_vs_TEM$avg_logFC < -0.4 & 
                                       CD4_T_Cytotoxic_vs_TEM$p_val_adj < 0.05] <- "TEM"
CD4_T_Cytotoxic_vs_TEM$delabel <- NA
CD4_T_Cytotoxic_vs_TEM$delabel[CD4_T_Cytotoxic_vs_TEM$diffexpressed != "NO"] <- 
  rownames(CD4_T_Cytotoxic_vs_TEM)[CD4_T_Cytotoxic_vs_TEM$diffexpressed != "NO"]
ggplot(CD4_T_Cytotoxic_vs_TEM, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+
  geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") + 
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")

# add a column of NAs
CD4_T_Exhausted_vs_TEM$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_T_Exhausted_vs_TEM$diffexpressed[CD4_T_Exhausted_vs_TEM$avg_logFC > 0.4 & 
                                       CD4_T_Exhausted_vs_TEM$p_val_adj < 0.05] <- "Exhausted"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_T_Exhausted_vs_TEM$diffexpressed[CD4_T_Exhausted_vs_TEM$avg_logFC < -0.4 & 
                                       CD4_T_Exhausted_vs_TEM$p_val_adj < 0.05] <- "TEM"
CD4_T_Exhausted_vs_TEM$delabel <- NA
CD4_T_Exhausted_vs_TEM$delabel[CD4_T_Exhausted_vs_TEM$diffexpressed != "NO"] <- 
  rownames(CD4_T_Exhausted_vs_TEM)[CD4_T_Exhausted_vs_TEM$diffexpressed != "NO"]
ggplot(CD4_T_Exhausted_vs_TEM, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+
  geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") + geom_hline(yintercept = -log10(0.05), col="red")+ 
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")

# add a column of NAs
CD4_T_aTregs_vs_rTregs$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_T_aTregs_vs_rTregs$diffexpressed[CD4_T_aTregs_vs_rTregs$avg_logFC > 0.4 & 
                                       CD4_T_aTregs_vs_rTregs$p_val_adj < 0.05] <- "aTregs"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_T_aTregs_vs_rTregs$diffexpressed[CD4_T_aTregs_vs_rTregs$avg_logFC < -0.4 & 
                                       CD4_T_aTregs_vs_rTregs$p_val_adj < 0.05] <- "rTregs"
CD4_T_aTregs_vs_rTregs$delabel <- NA
CD4_T_aTregs_vs_rTregs$delabel[CD4_T_aTregs_vs_rTregs$diffexpressed != "NO"] <- 
  rownames(CD4_T_aTregs_vs_rTregs)[CD4_T_aTregs_vs_rTregs$diffexpressed != "NO"]
ggplot(CD4_T_aTregs_vs_rTregs, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+
  geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") + geom_hline(yintercept = -log10(0.05), col="red")+ 
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")

saveRDS(CD4_T, file = "TMS_merge_CD4_clustered_names.rds")
CD4_T <- readRDS(file = "TMS_merge_CD4_clustered_names.rds")

write.csv(CD8_T@meta.data, "CD8_T_metedata_TMS.csv")
write.csv(CD4_T@meta.data, "CD4_T_metedata_TMS.csv")

CD4_T_lists <- clustify(CD4_T, ref_mat = ref_mouse.rnaseq, cluster_col = "subsets", obj_out = TRUE, threshold = 0.45)
#T_cells.markers <- FindAllMarkers(T_cells)

#Highly variable exposed genes
CD8_T <- FindVariableFeatures(CD8_T, selection.method = "mean.var.plot", dispersion.cutoff = c(0.6, Inf), 
                              mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(CD8_T),15)
plot1 <- VariableFeaturePlot(CD8_T)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(CD8_T), file = "TMS_merge_CD8_HVG.csv")

saveRDS(CD8_T, file = "TMS_merge_CD8_HVG.rds")
CD8_T <- readRDS(file = "TMS_merge_CD8_HVG.rds")

all.genes <- rownames(CD8_T)
CD8_T <- ScaleData(CD8_T, features = all.genes)

saveRDS(CD8_T, file = "TMS_merge_CD8_Scaled.rds")
CD8_T <- readRDS(file = "TMS_merge_CD8_Scaled.rds")

CD8_T <- RunPCA(CD8_T, npcs = 20, ndims.print = 1:20, nfeatures.print = 5)
DimPlot(CD8_T, reduction = "pca")

CD8_T <- JackStraw(CD8_T, num.replicate = 100, dims = 20)
CD8_T <- ScoreJackStraw(CD8_T, dims = 1:20)
JackStrawPlot(CD8_T, dims = 1:20)

saveRDS(CD8_T, file = "TMS_merge_CD8_T_pca.rds")
CD8_T <- readRDS(file = "TMS_merge_CD8_T_pca.rds")

ElbowPlot(CD8_T, ndims = 20)

CD8_T <- RunTSNE(CD8_T, dims = 1:20, perplexity = 30)
DimPlot(CD8_T, reduction = "tsne")

CD8_T <- RunUMAP(CD8_T, dims = 1:20)
DimPlot(CD8_T, reduction = "umap")

saveRDS(CD8_T, file = "TMS_merge_CD8_tsne_umap.rds")
CD8_T <- readRDS(file = "TMS_merge_CD4_tsne_umap.rds")

CD8_T <- FindNeighbors(CD8_T, dims = 1:20)

#louvian
CD8_T <- FindClusters(CD8_T, resolution = seq(0, 3, by = 0.1))

clustree(CD8_T, prefix = "RNA_snn_res.")
clustree(CD8_T, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(CD8_T) <- CD8_T$RNA_snn_res.0.8

CD8_T <- BuildClusterTree(CD8_T, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = CD8_T)

CD8_T.markers <- FindAllMarkers(CD8_T, only.pos = TRUE)
m <- CD8_T.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(CD4_T, features = m$gene)

write.csv(CD8_T.markers, "CD8_markers_tms_merge.csv")

CD8_T <- RenameIdents(CD8_T, "17" = "rTregs", "9" = "rTregs", "11" = "Tregs", "")

CD8_T <- BuildClusterTree(CD8_T, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = CD8_T)

saveRDS(CD8_T, file = "TMS_merge_CD8_clustered.rds")
CD8_T <- readRDS(file = "TMS_merge_CD8_clustered.rds")

CD8_T.marker_age <- list()

for (n in 1:12) {
  CD8_T.marker_age[[n]] <- FindConservedMarkers(CD8_T, ident.1 = n, grouping.var = "Age_group", min.cell.group = 0, 
                                                min.cell.feature = 0)
  write.csv(CD8_T.marker_age[[n]], sprintf("CD4_T_marker_age_T_%s.csv", n))
}

m <- CD4_T.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(CD4_T, features = m$gene)

CD4_T <- RenameIdents(CD4_T, "8" = "3", "5" = "3", )

CD4_T$subsets <- CD4_T@active.ident

CD4_T <- clustify(CD4_T, ref_mat = ref_mouse.rnaseq, cluster_col = "RNA_snn_res.10", obj_out = TRUE, 
                  threshold = 0.45)
#T_cells.markers <- FindAllMarkers(T_cells)

#monocle3
CDS <- as.cell_data_set(CD4_T)
CDS <- cluster_cells(cds = CDS, reduction_method = "UMAP")
plot_cells(cds = CDS, label_groups_by_cluster = FALSE, cell_size = 1, color_cells_by = "Age_group_2")
CDS <- learn_graph(CDS, use_partition = FALSE)
cell_ids <- colnames(CDS)[CDS$ident ==  "0"]
closest_vertex <- CDS@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(CDS), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(CDS)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]
CDS <- order_cells(cds = CDS, root_pr_nodes = root_pr_nodes)
plot_cells(CDS, color_cells_by = "pseudotime", graph_label_size = 4)
plot_cells(CDS, color_cells_by = "Age_group_2", group_label_size = 5, graph_label_size = 4)
de_res <- graph_test(CDS, neighbor_graph = "principal_graph", cores = 3)




CD4_Cytotoxic <- FindMarkers(CD4_T, ident.1 = "Cytotoxic", logfc.threshold = -Inf, grouping.var = "subsets")
# add a column of NAs
CD4_Cytotoxic$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_Cytotoxic$diffexpressed[CD4_Cytotoxic$avg_logFC > 0.4 & CD4_Cytotoxic$p_val_adj < 0.05] <- "Cytotoxic"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_Cytotoxic$diffexpressed[CD4_Cytotoxic$avg_logFC < -0.5 & 
                              CD4_Cytotoxic$p_val_adj < 0.05] <- "CD4 T cells not Cytotoxic"
CD4_Cytotoxic$delabel <- NA
CD4_Cytotoxic$delabel[CD4_Cytotoxic$diffexpressed != "NO"] <- 
  rownames(CD4_Cytotoxic)[CD4_Cytotoxic$diffexpressed != "NO"]
ggplot(CD4_Cytotoxic, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.5, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(CD4_Cytotoxic, "CD4_Cytotoxic_TMS_DEG.csv")


#monocle3
CDS <- as.cell_data_set(CD4_T)
CDS <- cluster_cells(cds = CDS, reduction_method = "UMAP")
plot_cells(cds = CDS, label_groups_by_cluster = FALSE, cell_size = 1, color_cells_by = "subsets")
CDS <- learn_graph(CDS, use_partition = FALSE)
cell_ids <- colnames(CDS)[CDS$subsets ==  "Naive"]
closest_vertex <- CDS@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(CDS), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(CDS)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]
CDS <- order_cells(cds = CDS, root_pr_nodes = root_pr_nodes)
plot_cells(CDS, color_cells_by = "pseudotime", graph_label_size = 4)
plot_cells(CDS, color_cells_by = "subsets", group_label_size = 5, graph_label_size = 4)
de_res <- graph_test(CDS, neighbor_graph = "principal_graph", cores = 3)

saveRDS(CDS, "CD4_T_TMS_CDS.rds")
