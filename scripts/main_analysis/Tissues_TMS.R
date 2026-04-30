#Libraries announcement
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
library(SeuratWrappers)
library(base)
library(renv)
library(ggplot2)
library(ggrepel)
library(clustree)
library(clustifyr)
library(clustifyrdata)

setwd("E:/Immunaging")
setwd("F:/Immunaging")
setwd("D:/Immunaging")

Tissue.list <- readRDS("Tissue_list.rds")

saveRDS(Tissue.list[["Tongue"]], "Tongue.rds")
saveRDS(Tissue.list[["Heart_and_Aorta"]], "Heart_and_Aorta.rds")
saveRDS(Tissue.list[["Marrow"]], "Marrow.rds")
saveRDS(Tissue.list[["Mammary_Gland"]], "Mammary_Gland.rds")
saveRDS(Tissue.list[["Fat"]], "Fat.rds")
saveRDS(Tissue.list[["Kidney"]], "Kidney.rds")
saveRDS(Tissue.list[["Liver"]], "Liver.rds")
saveRDS(Tissue.list[["Lung"]], "Lung.rds")
saveRDS(Tissue.list[["Limb_Muscle"]], "Limb_Muscle.rds")
saveRDS(Tissue.list[["Pancreas"]], "Pancreas.rds")
saveRDS(Tissue.list[["Spleen"]], "Spleen.rds")
saveRDS(Tissue.list[["Thymus"]], "Thymus.rds")
saveRDS(Tissue.list[["Bladder"]], "Bladder.rds")
saveRDS(Tissue.list[["Skin"]], "Skin.rds")
saveRDS(Tissue.list[["Large_Intestine"]], "Large_Intestine.rds")
saveRDS(Tissue.list[["Trachea"]], "Trachea.rds")

Tongue <- readRDS( file = "Tongue.rds")
Heart_and_Aorta <- readRDS( file = "Heart_and_Aorta.rds")
Marrow <- readRDS( file = "Marrow.rds")
Mammary_Gland <- readRDS( file = "Mammary_Gland.rds")
Fat <- readRDS( file = "Fat.rds")
Kidney <- readRDS( file = "Kidney.rds")
Liver <- readRDS( file = "Liver.rds")
Lung <- readRDS( file = "Lung.rds")
Limb_Muscle <- readRDS( file = "Limb_Muscle.rds")
Pancreas <- readRDS( file = "Pancreas.rds")
Spleen <- readRDS( file = "Spleen.rds")
Thymus <- readRDS( file = "Thymus.rds")
Bladder <- readRDS( file = "Bladder.rds")
Skin <- readRDS( file = "Skin.rds")
Large_Intestine <- readRDS( file = "Large_Intestine.rds")
Trachea <- readRDS( file = "Trachea.rds")

res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)

#Tongue
#Highly variable expressed genes
Tongue <- FindVariableFeatures(Tongue, selection.method = "mean.var.plot", dispersion.cutoff = c(0.6, Inf), 
                               mean.cutoff = c(0, 10))#HVG by TMS
top10 <- head(VariableFeatures(Tongue),15)
plot1 <- VariableFeaturePlot(Tongue)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Tongue), file = "Highly_Variable_Expessed_Genes_Tongue.csv")

saveRDS(Tongue, file = "VariableGenes_Tongue.rds")
Tongue <- readRDS(file = "VariableGenes_Tongue.rds")

all.genes <- rownames(Tongue)
Tongue <- ScaleData(Tongue, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Tongue, file = "ScaledData_Tongue.rds")
Tongue <- readRDS(file = "ScaledData_Tongue.rds")

Tongue <- RunPCA(Tongue, npcs = 40, ndims.print = 1:40, nfeatures.print = 5)
DimPlot(Tongue, reduction = "pca")

saveRDS(Tongue, file = "PCA_Tongue.rds")
Tongue <- readRDS(file = "PCA_Tongue.rds")

Tongue <- JackStraw(Tongue, num.replicate = 100, dims = 40)
Tongue <- ScoreJackStraw(Tongue, dims = 1:40)
JackStrawPlot(Tongue, dims = 1:40)

saveRDS(Tongue, file = "JackStraw_Tongue.rds")
Tongue <- readRDS(file = "JackStraw_Tongue.rds")

ElbowPlot(Tongue, ndims = 40)

Tongue <- FindNeighbors(Tongue, dims = 1:40)

#louvian
Tongue <- FindClusters(Tongue, resolution = res)

clustree(Tongue, prefix = "RNA_snn_res.")

Idents(Tongue) <- Tongue$RNA_snn_res.0.8

Tongue <- BuildClusterTree(Tongue, reorder.numeric = TRUE, reorder = TRUE, dims = 1:40)
PlotClusterTree(object = Tongue)

Tongue <- RunTSNE(Tongue, dims = 1:40, perplexity = 50)
DimPlot(Tongue, reduction = "tsne")

Tongue <- RunUMAP(Tongue, dims = 1:40)
DimPlot(Tongue, reduction = "umap")

Tongue.markers <- FindAllMarkers(Tongue)
write.csv(Tongue.markers, "Tongue_markers.csv")

Tongue <- RenameIdents(Tongue, "2" = "T cell")

Tongue$clusters_reordered <- Tongue@active.ident

Tongue <- clustify(Tongue, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", seurat_output = TRUE, 
                   threshold = 0.65)

Tongue_lists <- clustify(Tongue, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                         seurat_out = FALSE, threshold = 0.65)

plot_cor_heatmap(Tongue_lists)

Idents(Tongue) <- Tongue$type
Tongue <- RenameIdents(Tongue, "DN1 thymic pro-T cell-Thymus-CLASH!" = "T cell-Tongue")

saveRDS(Tongue, file = "Clustred_Tongue.rds")
Tongue <- readRDS(file = "Clustred_Tongue.rds")

#Heart and aorta
#Highly variable expressed genes
Heart_and_Aorta <- FindVariableFeatures(Heart_and_Aorta, selection.method = "mean.var.plot", 
                                        dispersion.cutoff = c(0, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Heart_and_Aorta),15)
plot1 <- VariableFeaturePlot(Heart_and_Aorta)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Heart_and_Aorta), file = "Highly_Variable_Expessed_Genes_Heart_and_Aorta.csv")

saveRDS(Heart_and_Aorta, file = "VariableGenes_Heart_and_Aorta.rds")
Heart_and_Aorta <- readRDS(file = "VariableGenes_Heart_and_Aorta.rds")

all.genes <- rownames(Heart_and_Aorta)
Heart_and_Aorta <- ScaleData(Heart_and_Aorta, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Heart_and_Aorta, file = "ScaledData_Heart_and_Aorta.rds")
Heart_and_Aorta <- readRDS(file = "ScaledData_Heart_and_Aorta.rds")

Heart_and_Aorta <- RunPCA(Heart_and_Aorta, npcs = 40, ndims.print = 1:40, nfeatures.print = 5)
DimPlot(Heart_and_Aorta, reduction = "pca")

saveRDS(Heart_and_Aorta, file = "PCA_Heart_and_Aorta.rds")
Heart_and_Aorta <- readRDS(file = "PCA_Heart_and_Aorta.rds")

Heart_and_Aorta <- JackStraw(Heart_and_Aorta, num.replicate = 100, dims = 40)
Heart_and_Aorta <- ScoreJackStraw(Heart_and_Aorta, dims = 1:40)
JackStrawPlot(Heart_and_Aorta, dims = 1:40)

saveRDS(Heart_and_Aorta, file = "JackStraw_Heart_and_Aorta.rds")
Heart_and_Aorta <- readRDS(file = "JackStraw_Heart_and_Aorta.rds")

ElbowPlot(Heart_and_Aorta, ndims = 40)

Heart_and_Aorta <- FindNeighbors(Heart_and_Aorta, dims = 1:40)

#louvian
Heart_and_Aorta <- FindClusters(Heart_and_Aorta, resolution = res)

clustree(Heart_and_Aorta, prefix = "RNA_snn_res.")

Idents(Heart_and_Aorta) <- Heart_and_Aorta$RNA_snn_res.0.6

Heart_and_Aorta <- BuildClusterTree(Heart_and_Aorta, reorder.numeric = TRUE, reorder = TRUE, dims = 1:40)
PlotClusterTree(object = Heart_and_Aorta)

Heart_and_Aorta <- RunTSNE(Heart_and_Aorta, dims = 1:40, perplexity = 50)
DimPlot(Heart_and_Aorta, reduction = "tsne")

Heart_and_Aorta <- RunUMAP(Heart_and_Aorta, dims = 1:40)
DimPlot(Heart_and_Aorta, reduction = "umap")

Heart_and_Aorta.markers <- FindAllMarkers(Heart_and_Aorta)
write.csv(Heart_and_Aorta.markers, "Heart_and_Aorta_markers.csv")

Heart_and_Aorta <- RenameIdents(Heart_and_Aorta, "7" = "1", "6" = "1")

Heart_and_Aorta$clusters_reordered <- Heart_and_Aorta@active.ident

Heart_and_Aorta <- clustify(Heart_and_Aorta, ref_mat = ref_tabula_muris_facs, cluster_col = "clusters_reordered", 
                            seurat_output = TRUE, threshold = 0.45)

Heart_and_Aorta_lists <- clustify(Heart_and_Aorta, ref_mat = ref_tabula_muris_facs, 
                                  cluster_col = "clusters_reordered", seurat_out = FALSE, threshold = 0.45)

plot_cor_heatmap(Heart_and_Aorta_lists)
Idents(Heart_and_Aorta) <- Heart_and_Aorta$type
Heart_and_Aorta <- RenameIdents(Heart_and_Aorta,
                                "luminal epithelial cell of mammary gland-Mammary" = "endothelial cell-Heart", 
                                "lung endothelial cell-Lung" = "endothelial cell-Heart", 
                                "leukocyte-Heart" = "monocyte-Heart", 
                                "brain pericyte-Brain" = "pericyte-Heart",
                                "endothelial cell-Fat" = "endothelial cell-Heart")

saveRDS(Heart_and_Aorta, file = "Clustred_Heart_and_Aorta.rds")
Heart_and_Aorta <- readRDS(file = "Clustred_Heart_and_Aorta.rds")

#Marrow
#Highly variable expressed genes
Marrow <- FindVariableFeatures(Marrow, selection.method = "mean.var.plot", dispersion.cutoff = c(0, Inf), 
                               mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Marrow),15)
plot1 <- VariableFeaturePlot(Marrow)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Marrow), file = "Highly_Variable_Expessed_Genes_Marrow.csv")

saveRDS(Marrow, file = "VariableGenes_Marrow.rds")
Marrow <- readRDS(file = "VariableGenes_Marrow.rds")

all.genes <- rownames(Marrow)
Marrow <- ScaleData(Marrow, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Marrow, file = "ScaledData_Marrow.rds")
Marrow <- readRDS(file = "ScaledData_Marrow.rds")

Marrow <- RunPCA(Marrow, npcs = 50, ndims.print = 1:50, nfeatures.print = 5)
DimPlot(Marrow, reduction = "pca")

saveRDS(Marrow, file = "PCA_Marrow.rds")
Marrow <- readRDS(file = "PCA_Marrow.rds")

Marrow <- JackStraw(Marrow, num.replicate = 100, dims = 50)
Marrow <- ScoreJackStraw(Marrow, dims = 1:50)
JackStrawPlot(Marrow, dims = 1:50)

saveRDS(Marrow, file = "JackStraw_Marrow.rds")
Marrow <- readRDS(file = "JackStraw_Marrow.rds")

ElbowPlot(Marrow, ndims = 50)

Marrow <- FindNeighbors(Marrow, dims = 1:50)

#louvian
Marrow <- FindClusters(Marrow, resolution = res)

clustree(Marrow, prefix = "RNA_snn_res.")

Idents(Marrow) <- Marrow$RNA_snn_res.0.6

Marrow <- RunTSNE(Marrow, dims = 1:50, perplexity = 50)
DimPlot(Marrow, reduction = "tsne")

Marrow <- RunUMAP(Marrow, dims = 1:50)
DimPlot(Marrow, reduction = "umap")

Marrow <- BuildClusterTree(Marrow, reorder.numeric = TRUE, reorder = TRUE, dims = 1:50)
PlotClusterTree(object = Marrow)

Marrow.markers <- FindAllMarkers(Marrow)
write.csv(Marrow.markers, "Marrow_markers.csv")

Marrow <- RenameIdents(Marrow, "24" = "T cell")

Marrow$clusters_reordered <- Marrow@active.ident

Marrow <- clustify(Marrow, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", seurat_output = TRUE,
                   threshold = 0.7)

Marrow_lists <- clustify(Marrow, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                         seurat_out = FALSE, threshold = 0.7)

plot_cor_heatmap(Marrow_lists)

Idents(Marrow) <- Marrow$type
Marrow <- RenameIdents(Marrow, "B cell-Spleen" = "B cell-Marrow")

saveRDS(Marrow, file = "Clustred_Marrow.rds")
Marrow <- readRDS(file = "Clustred_Marrow.rds")

#Mammary_Gland
#Highly variable expressed genes
Mammary_Gland <- FindVariableFeatures(Mammary_Gland, selection.method = "mean.var.plot", 
                                      dispersion.cutoff = c(0.2, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Mammary_Gland),15)
plot1 <- VariableFeaturePlot(Mammary_Gland)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Mammary_Gland), file = "Highly_Variable_Expessed_Genes_Mammary_Gland.csv")

saveRDS(Mammary_Gland, file = "VariableGenes_Mammary_Gland.rds")
Mammary_Gland <- readRDS(file = "VariableGenes_Mammary_Gland.rds")

all.genes <- rownames(Mammary_Gland)
Mammary_Gland <- ScaleData(Mammary_Gland, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Mammary_Gland, file = "ScaledData_Mammary_Gland.rds")
Mammary_Gland <- readRDS(file = "ScaledData_Mammary_Gland.rds")

Mammary_Gland <- RunPCA(Mammary_Gland, npcs = 50, ndims.print = 1:50, nfeatures.print = 5)
DimPlot(Mammary_Gland, reduction = "pca")

saveRDS(Mammary_Gland, file = "PCA_Mammary_Gland.rds")
Mammary_Gland <- readRDS(file = "PCA_Mammary_Gland.rds")

Mammary_Gland <- JackStraw(Mammary_Gland, num.replicate = 100, dims = 50)
Mammary_Gland <- ScoreJackStraw(Mammary_Gland, dims = 1:50)
JackStrawPlot(Mammary_Gland, dims = 1:50)

saveRDS(Mammary_Gland, file = "JackStraw_Mammary_Gland.rds")
Mammary_Gland <- readRDS(file = "JackStraw_Mammary_Gland.rds")

ElbowPlot(Mammary_Gland, ndims = 50)

Mammary_Gland <- FindNeighbors(Mammary_Gland, dims = 1:50)

#louvian
Mammary_Gland <- FindClusters(Mammary_Gland, resolution = res)

clustree(Mammary_Gland, prefix = "RNA_snn_res.")

Idents(Mammary_Gland) <- Mammary_Gland$RNA_snn_res.0.6

Mammary_Gland <- BuildClusterTree(Mammary_Gland, reorder.numeric = TRUE, reorder = TRUE, dims = 1:50)
PlotClusterTree(object = Mammary_Gland)

Mammary_Gland.markers <- FindAllMarkers(Mammary_Gland)
write.csv(Mammary_Gland.markers, "Mammary_Gland_markers.csv")

Mammary_Gland <- RunTSNE(Mammary_Gland, dims = 1:50, perplexity = 50)
DimPlot(Mammary_Gland, reduction = "tsne")

Mammary_Gland <- RunUMAP(Mammary_Gland, dims = 1:50)
DimPlot(Mammary_Gland, reduction = "umap")

Mammary_Gland <- RenameIdents(Mammary_Gland, "11" = "2", "10" = "2", "9" = "2", "8" = "2", "7" = "2", "6" = "2",
                              "3" = "2")

Mammary_Gland$clusters_reordered <- Mammary_Gland@active.ident

Mammary_Gland <- clustify(Mammary_Gland, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                          seurat_output = TRUE, threshold = 0.65)

Mammary_Gland_lists <- clustify(Mammary_Gland, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                                seurat_out = FALSE, threshold = 0.65)

plot_cor_heatmap(Mammary_Gland_lists)

Idents(Mammary_Gland) <- Mammary_Gland$type
Mammary_Gland <- RenameIdents(Mammary_Gland, "macrophage-Marrow" = "macrophage-Mammary",
                              "mesenchymal stem cell-Limb" = "mesenchymal stem cell-Mammary",
                              "unknown-Limb" = "smooth muscle cell-Mammary")

saveRDS(Mammary_Gland, file = "Clustred_Mammary_Gland.rds")
Mammary_Gland <- readRDS(file = "Clustred_Mammary_Gland.rds")

#Fat
#Highly variable expressed genes
Fat <- FindVariableFeatures(Fat, selection.method = "mean.var.plot", 
                            dispersion.cutoff = c(0, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Fat),15)
plot1 <- VariableFeaturePlot(Fat)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Fat), file = "Highly_Variable_Expessed_Genes_Fat.csv")

saveRDS(Fat, file = "VariableGenes_Fat.rds")
Fat <- readRDS(file = "VariableGenes_Fat.rds")

all.genes <- rownames(Fat)
Fat <- ScaleData(Fat, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Fat, file = "ScaledData_Fat.rds")
Fat <- readRDS(file = "ScaledData_Fat.rds")

Fat <- RunPCA(Fat, npcs = 50, ndims.print = 1:50, nfeatures.print = 5)
DimPlot(Fat, reduction = "pca")

saveRDS(Fat, file = "PCA_Fat.rds")
Fat <- readRDS(file = "PCA_Fat.rds")

Fat <- JackStraw(Fat, num.replicate = 100, dims = 50)
Fat <- ScoreJackStraw(Fat, dims = 1:50)
JackStrawPlot(Fat, dims = 1:50)

saveRDS(Fat, file = "JackStraw_Fat.rds")
Fat <- readRDS(file = "JackStraw_Fat.rds")

ElbowPlot(Fat, ndims = 50)

Fat <- FindNeighbors(Fat, dims = 1:50)

#louvian
Fat <- FindClusters(Fat, resolution = res)

clustree(Fat, prefix = "RNA_snn_res.")

Idents(Fat) <- Fat$RNA_snn_res.0.4

Fat <- BuildClusterTree(Fat, reorder.numeric = TRUE, reorder = TRUE, dims = 1:50)
PlotClusterTree(object = Fat)

Fat <- RunTSNE(Fat, dims = 1:50, perplexity = 50)
DimPlot(Fat, reduction = "tsne")

Fat <- RunUMAP(Fat, dims = 1:50)
DimPlot(Fat, reduction = "umap")

Fat.markers <- FindAllMarkers(Fat)
write.csv(Fat.markers, "Fat_markers.csv")

Fat <- RenameIdents(Fat, "16" = "15")

Fat$clusters_reordered <- Fat@active.ident

Fat <- clustify(Fat, ref_mat = ref_tabula_muris_facs, cluster_col = "clusters_reordered", seurat_output = TRUE, 
                threshold = 0.5)

Fat_lists <- clustify(Fat, ref_mat = ref_tabula_muris_facs, cluster_col = "clusters_reordered", 
                      seurat_out = FALSE, threshold = 0.5)

plot_cor_heatmap(Fat_lists)

Idents(Fat) <- Fat$type
Fat <- RenameIdents(Fat, "mesenchymal stem cell of adipose-Fat-CLASH!" = "mesenchymal stem cell of adipose-Fat",
                    "myeloid cell-Fat" = "monocyte-Fat",
                    "stromal cell-Mammary" = "stromal cell-Fat",
                    "leukocyte-Heart" = "monocyte-Fat",
                    "mesenchymal cell-Trachea" = "mesenchymal cell-Fat",
                    "luminal epithelial cell of mammary gland-Mammary" = "epithelial cell-Fat", 
                    "granulocyte-Marrow" = "granulocyte-Fat")

saveRDS(Fat, file = "Clustred_Fat.rds")
Fat <- readRDS(file = "Clustred_Fat.rds")

#Kidney
#Highly variable expressed genes
Kidney <- FindVariableFeatures(Kidney, selection.method = "mean.var.plot", 
                               dispersion.cutoff = c(0.3, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Kidney),15)
plot1 <- VariableFeaturePlot(Kidney)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Kidney), file = "Highly_Variable_Expessed_Genes_Kidney.csv")

saveRDS(Kidney, file = "VariableGenes_Kidney.rds")
Kidney <- readRDS(file = "VariableGenes_Kidney.rds")

all.genes <- rownames(Kidney)
Kidney <- ScaleData(Kidney, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Kidney, file = "ScaledData_Kidney.rds")
Kidney <- readRDS(file = "ScaledData_Kidney.rds")

Kidney <- RunPCA(Kidney, npcs = 80)
DimPlot(Kidney, reduction = "pca")

saveRDS(Kidney, file = "PCA_Kidney.rds")
Kidney <- readRDS(file = "PCA_Kidney.rds")

Kidney <- JackStraw(Kidney, num.replicate = 100, dims = 80)
Kidney <- ScoreJackStraw(Kidney, dims = 1:80)
JackStrawPlot(Kidney, dims = 1:80)

saveRDS(Kidney, file = "JackStraw_Kidney.rds")
Kidney <- readRDS(file = "JackStraw_Kidney.rds")

ElbowPlot(Kidney, ndims = 80)

Kidney <- FindNeighbors(Kidney, dims = 1:80)

#louvian
Kidney <- FindClusters(Kidney, resolution = res)

clustree(Kidney, prefix = "RNA_snn_res.")

Idents(Kidney) <- Kidney$RNA_snn_res.0.2

Kidney <- BuildClusterTree(Kidney, reorder.numeric = TRUE, reorder = TRUE, dims = 1:80)
PlotClusterTree(object = Kidney)

Kidney.markers <- FindAllMarkers(Kidney)
write.csv(Kidney.markers, "Kidney_markers.csv")

Kidney <- RunTSNE(Kidney, dims = 1:80, perplexity = 50)
DimPlot(Kidney, reduction = "tsne")

Kidney <- RunUMAP(Kidney, dims = 1:80)
DimPlot(Kidney, reduction = "umap")

Kidney <- RenameIdents(Kidney, "30" = "28", "29" = "28")

Kidney$clusters_reordered <- Kidney@active.ident

Kidney <- clustify(Kidney, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", seurat_output = TRUE, 
                   threshold = 0.45)

Kidney_lists <- clustify(Kidney, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                         seurat_out = FALSE, threshold = 0.45)

plot_cor_heatmap(Kidney_lists)

Idents(Kidney) <- Kidney$type
Kidney <- RenameIdents(Kidney, "B cell-Spleen" = "B cell-Kidney",
                       "macrophage-Spleen" = "macrophage-Kidney",
                       "granulocyte-Marrow" = "granulocyte-Kidney",
                       "T cell-Spleen" = "T cell-Kidney",
                       "proerythroblast-Marrow" = "proerythroblast-Kidney")

saveRDS(Kidney, file = "Clustred_Kidney.rds")
Kidney <- readRDS(file = "Clustred_Kidney.rds")

#Liver
#Highly variable expressed genes
Liver <- FindVariableFeatures(Liver, selection.method = "mean.var.plot", 
                              dispersion.cutoff = c(0.3, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Liver),15)
plot1 <- VariableFeaturePlot(Liver)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Liver), file = "Highly_Variable_Expessed_Genes_Liver.csv")

saveRDS(Liver, file = "VariableGenes_Liver.rds")
Liver <- readRDS(file = "VariableGenes_Liver.rds")

all.genes <- rownames(Liver)
Liver <- ScaleData(Liver, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Liver, file = "ScaledData_Liver.rds")
Liver <- readRDS(file = "ScaledData_Liver.rds")

Liver <- RunPCA(Liver, npcs = 40)
DimPlot(Liver, reduction = "pca")

saveRDS(Liver, file = "PCA_Liver.rds")
Liver <- readRDS(file = "PCA_Liver.rds")

Liver <- JackStraw(Liver, num.replicate = 100, dims = 40)
Liver <- ScoreJackStraw(Liver, dims = 1:40)
JackStrawPlot(Liver, dims = 1:40)

saveRDS(Liver, file = "JackStraw_Liver.rds")
Liver <- readRDS(file = "JackStraw_Liver.rds")

ElbowPlot(Liver, ndims = 40)

Liver <- FindNeighbors(Liver, dims = 1:40)

#louvian
Liver <- FindClusters(Liver, resolution = res)

clustree(Liver, prefix = "RNA_snn_res.")

Idents(Liver) <- Liver$RNA_snn_res.0.5

Liver <- BuildClusterTree(Liver, reorder.numeric = TRUE, reorder = TRUE, dims = 1:40)
PlotClusterTree(object = Liver)

Liver.markers <- FindAllMarkers(Liver)
write.csv(Liver.markers, "Liver_markers.csv")

Liver <- RunTSNE(Liver, dims = 1:40, perplexity = 50)
DimPlot(Liver, reduction = "tsne")

Liver <- RunUMAP(Liver, dims = 1:40)
DimPlot(Liver, reduction = "umap")

Liver <- RenameIdents(Liver, "28" = "26", "27" = "26")

Liver$clusters_reordered <- Liver@active.ident

Liver <- clustify(Liver, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", seurat_output = TRUE, 
                  threshold = 0.4)

Liver_lists <- clustify(Liver, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                         seurat_out = FALSE, threshold = 0.45)

plot_cor_heatmap(Liver_lists)

Idents(Liver) <- Liver$type
Liver <- RenameIdents(Liver, "leukocyte-Liver" = "monocyte-Liver",
                      "macrophage-Limb" = "macrophage-Liver",
                      "macrophage-Spleen" = "macrophage-Liver",
                      "macrophage-Marrow" = "macrophage-Liver",
                      "T cell-Spleen" = "T cell-Liver",
                      "B cell-Spleen" = "B cell-Liver",
                      "neuroendocrine cell-Trachea" = "neuroendocrine cell-Liver",
                      "midlobular male-Liver" = "hepatocyte-Liver",
                      "macrophage-Kidney" = "macrophage-Liver",
                      "monocyte-Marrow" = "monocyte-Liver")

saveRDS(Liver, file = "Clustred_Liver.rds")
Liver <- readRDS(file = "Clustred_Liver.rds")

#Lung
#Highly variable expressed genes
Lung <- FindVariableFeatures(Lung, selection.method = "mean.var.plot",                                
                             dispersion.cutoff = c(0, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Lung),15)
plot1 <- VariableFeaturePlot(Lung)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Lung), file = "Highly_Variable_Expessed_Genes_Lung.csv")

saveRDS(Lung, file = "VariableGenes_Lung.rds")
Lung <- readRDS(file = "VariableGenes_Lung.rds")

all.genes <- rownames(Lung)
Lung <- ScaleData(Lung, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Lung, file = "ScaledData_Lung.rds")
Lung <- readRDS(file = "ScaledData_Lung.rds")

Lung <- RunPCA(Lung, npcs = 40)
DimPlot(Lung, reduction = "pca")

saveRDS(Lung, file = "PCA_Lung.rds")
Lung <- readRDS(file = "PCA_Lung.rds")

Lung <- JackStraw(Lung, num.replicate = 100, dims = 40)
Lung <- ScoreJackStraw(Lung, dims = 1:40)
JackStrawPlot(Lung, dims = 1:40)

saveRDS(Lung, file = "JackStraw_Lung.rds")
Lung <- readRDS(file = "JackStraw_Lung.rds")

ElbowPlot(Lung, ndims = 40)

Lung <- FindNeighbors(Lung, dims = 1:40)

#louvian
Lung <- FindClusters(Lung, resolution = res)

clustree(Lung, prefix = "RNA_snn_res.")

Idents(Lung) <- Lung$RNA_snn_res.0.8

Lung <- BuildClusterTree(Lung, reorder.numeric = TRUE, reorder = TRUE, dims = 1:40)
PlotClusterTree(object = Lung)

Lung <- RunTSNE(Lung, dims = 1:40, perplexity = 50)
DimPlot(Lung, reduction = "tsne")

Lung <- RunUMAP(Lung, dims = 1:40)
DimPlot(Lung, reduction = "umap")

Lung.markers <- FindAllMarkers(Lung)
write.csv(Lung.markers, "Lung_markers.csv")

Lung <- RenameIdents(Lung, "46" = "10", "45" = "10", "44" = "10", "43" = "10", "42" = "10", "41" = "10",
                     "11" = "10")

Lung$clusters_reordered <- Lung@active.ident

Lung <- clustify(Lung, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", seurat_output = TRUE, 
                 threshold = 0.6)

Lung_lists <- clustify(Lung, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                        seurat_out = FALSE, threshold = 0.6)

plot_cor_heatmap(Lung_lists)

Idents(Lung) <- Lung$type
Lung <- RenameIdents(Lung, "non-classical monocyte-Lung-CLASH!" = "non-classical monocyte-Lung",
                     "myeloid cell-Lung-CLASH!" = "dendritic cell-Lung",
                     "classical monocyte-Lung-CLASH!" = "classical monocyte-Lung",
                     "monocyte-Marrow" = "monocyte-Lung",
                     "type II pneumocyte-Lung-CLASH!" = "type II pneumocyte-Lung",
                     "ciliated columnar cell of tracheobronchial tree-Lung-CLASH!" = "ciliated columnar cell of tracheobronchial tree-Lung",
                     "leukocyte-Lung" = "monocyte-Lung",
                     "B cell-Spleen" = "B cell-Lung", 
                     "basophil-Marrow" = "basophil-Lung",
                     "promonocyte-Marrow" = "monocyte-Lung",
                     "macrophage-Kidney" = "macrophage-Lung")

saveRDS(Lung, file = "Clustred_Lung.rds")
Lung <- readRDS(file = "Clustred_Lung.rds")

#Limb_Muscle
#Highly variable expressed genes
Limb_Muscle <- FindVariableFeatures(Limb_Muscle, selection.method = "mean.var.plot",                                
                                    dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Limb_Muscle),15)
plot1 <- VariableFeaturePlot(Limb_Muscle)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Limb_Muscle), file = "Highly_Variable_Expessed_Genes_Limb_Muscle.csv")

saveRDS(Limb_Muscle, file = "VariableGenes_Limb_Muscle.rds")
Limb_Muscle <- readRDS(file = "VariableGenes_Limb_Muscle.rds")

all.genes <- rownames(Limb_Muscle)
Limb_Muscle <- ScaleData(Limb_Muscle, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Limb_Muscle, file = "ScaledData_Limb_Muscle.rds")
Limb_Muscle <- readRDS(file = "ScaledData_Limb_Muscle.rds")

Limb_Muscle <- RunPCA(Limb_Muscle, npcs = 80)
DimPlot(Limb_Muscle, reduction = "pca")

saveRDS(Limb_Muscle, file = "PCA_Limb_Muscle.rds")
Limb_Muscle <- readRDS(file = "PCA_Limb_Muscle.rds")

Limb_Muscle <- JackStraw(Limb_Muscle, num.replicate = 100, dims = 80)
Limb_Muscle <- ScoreJackStraw(Limb_Muscle, dims = 1:80)
JackStrawPlot(Limb_Muscle, dims = 1:80)

saveRDS(Limb_Muscle, file = "JackStraw_Limb_Muscle.rds")
Limb_Muscle <- readRDS(file = "JackStraw_Limb_Muscle.rds")

ElbowPlot(Limb_Muscle, ndims = 80)

Limb_Muscle <- FindNeighbors(Limb_Muscle, dims = 1:80)

#louvian
Limb_Muscle <- FindClusters(Limb_Muscle, resolution = res)

clustree(Limb_Muscle, prefix = "RNA_snn_res.")

Idents(Limb_Muscle) <- Limb_Muscle$RNA_snn_res.0.6

Limb_Muscle <- BuildClusterTree(Limb_Muscle, reorder.numeric = TRUE, reorder = TRUE, dims = 1:80)
PlotClusterTree(object = Limb_Muscle)

Limb_Muscle.markers <- FindAllMarkers(Limb_Muscle)
write.csv(Limb_Muscle.markers, "Limb_Muscle_markers.csv")

Limb_Muscle <- RunTSNE(Limb_Muscle, dims = 1:80, perplexity = 50)
DimPlot(Limb_Muscle, reduction = "tsne")

Limb_Muscle <- RunUMAP(Limb_Muscle, dims = 1:80)
DimPlot(Limb_Muscle, reduction = "umap")

Limb_Muscle <- RenameIdents(Limb_Muscle, "6" = "4", "5" = "4")

Limb_Muscle$clusters_reordered <- Limb_Muscle@active.ident

Limb_Muscle <- clustify(Limb_Muscle, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                        seurat_output = TRUE, threshold = 0.6)

Limb_Muscle_lists <- clustify(Limb_Muscle, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                       seurat_out = FALSE, threshold = 0.6)

plot_cor_heatmap(Limb_Muscle_lists)

Idents(Limb_Muscle) <- Limb_Muscle$type
Limb_Muscle <- RenameIdents(Limb_Muscle, "endothelial cell-Bladder" = "endothelial cell-Limb",
                            "leukocyte-Lung" = "macrophage-Limb",
                            "unknown-Limb" = "smooth muscle cell-Limb",
                            "neuroendocrine cell-Trachea" = "neuroendocrine cell-Limb",
                            "endothelial cell-Mammary" = "endothelial cell-Limb")

saveRDS(Limb_Muscle, file = "Clustred_Limb_Muscle.rds")
Limb_Muscle <- readRDS(file = "Clustred_Limb_Muscle.rds")

#Pancreas
#Highly variable expressed genes
Pancreas <- FindVariableFeatures(Pancreas, selection.method = "mean.var.plot",                                
                                 dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Pancreas),15)
plot1 <- VariableFeaturePlot(Pancreas)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Pancreas), file = "Highly_Variable_Expessed_Genes_Pancreas.csv")

saveRDS(Pancreas, file = "VariableGenes_Pancreas.rds")
Pancreas <- readRDS(file = "VariableGenes_Pancreas.rds")

all.genes <- rownames(Pancreas)
Pancreas <- ScaleData(Pancreas, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Pancreas, file = "ScaledData_Pancreas.rds")
Pancreas <- readRDS(file = "ScaledData_Pancreas.rds")

Pancreas <- RunPCA(Pancreas, npcs = 30)
DimPlot(Pancreas, reduction = "pca")

saveRDS(Pancreas, file = "PCA_Pancreas.rds")
Pancreas <- readRDS(file = "PCA_Pancreas.rds")

Pancreas <- JackStraw(Pancreas, num.replicate = 100, dims = 30)
Pancreas <- ScoreJackStraw(Pancreas, dims = 1:30)
JackStrawPlot(Pancreas, dims = 1:30)

saveRDS(Pancreas, file = "JackStraw_Pancreas.rds")
Pancreas <- readRDS(file = "JackStraw_Pancreas.rds")

ElbowPlot(Pancreas, ndims = 30)

Pancreas <- FindNeighbors(Pancreas, dims = 1:30)

#louvian
Pancreas <- FindClusters(Pancreas, resolution = res)

clustree(Pancreas, prefix = "RNA_snn_res.")

Idents(Pancreas) <- Pancreas$RNA_snn_res.0.6

Pancreas <- BuildClusterTree(Pancreas, reorder.numeric = TRUE, reorder = TRUE, dims = 1:30)
PlotClusterTree(object = Pancreas)

Pancreas.markers <- FindAllMarkers(Pancreas)
write.csv(Pancreas.markers, "Pancreas_markers.csv")

Pancreas <- RunTSNE(Pancreas, dims = 1:30, perplexity = 50)
DimPlot(Pancreas, reduction = "tsne")

Pancreas <- RunUMAP(Pancreas, dims = 1:30)
DimPlot(Pancreas, reduction = "umap")

Pancreas <- RenameIdents(Pancreas, "20" = "19")

Pancreas$clusters_reordered <- Pancreas@active.ident

Pancreas <- clustify(Pancreas, ref_mat = ref_tabula_muris_facs, cluster_col = "clusters_reordered", 
                     seurat_output = TRUE, threshold = 0.55)

Pancreas_lists <- clustify(Pancreas, ref_mat = ref_tabula_muris_facs, cluster_col = "clusters_reordered", 
                              seurat_out = FALSE, threshold = 0.55)

plot_cor_heatmap(Pancreas_lists)

Idents(Pancreas) <- Pancreas$type
Pancreas <- RenameIdents(Pancreas, "type B pancreatic cell-Pancreas-CLASH!" = "type B pancreatic cell-Pancreas",
                         "pancreatic D cell-Pancreas-CLASH!" = "pancreatic D cell-Pancreas",
                         "pancreatic A cell-Pancreas-CLASH!" = "pancreatic A cell-Pancreas",
                         "pancreatic acinar cell-Pancreas-CLASH!" = "pancreatic acinar cell-Pancreas",
                         "pancreatic PP cell-Pancreas-CLASH!" = "pancreatic PP cell-Pancreas",
                         "myeloid cell-Fat" = "macrophage-Pancreas",
                         "pancreatic stellate cell-Pancreas-CLASH!" = "pancreatic stellate cell-Pancreas",
                         "immature B cell-Marrow" = "B cell-Pancreas",
                         "T cell-Fat" = "T cell-Pancreas",
                         "pancreatic ductal cell-Pancreas-CLASH!" = "pancreatic ductal cell-Pancreas")

saveRDS(Pancreas, file = "Clustred_Pancreas.rds")
Pancreas <- readRDS(file = "Clustred_Pancreas.rds")

#Spleen
#Highly variable expessed genes
Spleen <- FindVariableFeatures(Spleen, selection.method = "mean.var.plot",                                
                               dispersion.cutoff = c(0, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Spleen),15)
plot1 <- VariableFeaturePlot(Spleen)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Spleen), file = "Highly_Variable_Expessed_Genes_Spleen.csv")

saveRDS(Spleen, file = "VariableGenes_Spleen.rds")
Spleen <- readRDS(file = "VariableGenes_Spleen.rds")

all.genes <- rownames(Spleen)
Spleen <- ScaleData(Spleen, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Spleen, file = "ScaledData_Spleen.rds")
Spleen <- readRDS(file = "ScaledData_Spleen.rds")

Spleen <- RunPCA(Spleen, npcs = 60)
DimPlot(Spleen, reduction = "pca")

saveRDS(Spleen, file = "PCA_Spleen.rds")
Spleen <- readRDS(file = "PCA_Spleen.rds")

Spleen <- JackStraw(Spleen, num.replicate = 100, dims = 60)
Spleen <- ScoreJackStraw(Spleen, dims = 1:60)
JackStrawPlot(Spleen, dims = 1:60)

saveRDS(Spleen, file = "JackStraw_Spleen.rds")
Spleen <- readRDS(file = "JackStraw_Spleen.rds")

ElbowPlot(Spleen, ndims = 60)

Spleen <- FindNeighbors(Spleen, dims = 1:60)

#louvian
Spleen <- FindClusters(Spleen, resolution = res)

clustree(Spleen, prefix = "RNA_snn_res.")

Idents(Spleen) <- Spleen$RNA_snn_res.0.6

Spleen <- BuildClusterTree(Spleen, reorder.numeric = TRUE, reorder = TRUE, dims = 1:60)
PlotClusterTree(object = Spleen)

Spleen.markers <- FindAllMarkers(Spleen)
write.csv(Spleen.markers, "Spleen_markers.csv")

Spleen <- RunTSNE(Spleen, dims = 1:60, perplexity = 50)
DimPlot(Spleen, reduction = "tsne")

Spleen <- RunUMAP(Spleen, dims = 1:60)
DimPlot(Spleen, reduction = "umap")

Spleen <- RenameIdents(Spleen, "14" = "6", "13" = "6", "12" = "6", "11" = "6", "10" = "6", "9" = "6", "8" = "6",
                       "7" = "6")

Spleen$clusters_reordered <- Spleen@active.ident

Spleen <- clustify(Spleen, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                   seurat_output = TRUE, threshold = 0.65)

Spleen_lists <- clustify(Spleen, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                       seurat_out = FALSE, threshold = 0.65)

plot_cor_heatmap(Spleen_lists)

Idents(Spleen) <- Spleen$type
Spleen <- RenameIdents(Spleen, "proerythroblast-Marrow" = "proerythroblast-Spleen",
                       "hematopoietic precursor cell-Marrow" = "hematopoietic precursor cell-Spleen",
                       "granulocyte-Marrow" = "granulocyte-Spleen",
                       "late pro-B cell-Marrow" = "B cell-Spleen",
                       "promonocyte-Marrow" = "monocyte-Spleen",
                       "macrophage-Marrow" = "macrophage-Spleen")

saveRDS(Spleen, file = "Clustred_Spleen.rds")
Spleen <- readRDS(file = "Clustred_Spleen.rds")

#Thymus
#Highly variable expressed genes
Thymus <- FindVariableFeatures(Thymus, selection.method = "mean.var.plot",                                
                               dispersion.cutoff = c(0, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Thymus),15)
plot1 <- VariableFeaturePlot(Thymus)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Thymus), file = "Highly_Variable_Expessed_Genes_Thymus.csv")

saveRDS(Thymus, file = "VariableGenes_Thymus.rds")
Thymus <- readRDS(file = "VariableGenes_Thymus.rds")

all.genes <- rownames(Thymus)
Thymus <- ScaleData(Thymus, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Thymus, file = "ScaledData_Thymus.rds")
Thymus <- readRDS(file = "ScaledData_Thymus.rds")

Thymus <- RunPCA(Thymus, npcs = 30)
DimPlot(Thymus, reduction = "pca")

saveRDS(Thymus, file = "PCA_Thymus.rds")
Thymus <- readRDS(file = "PCA_Thymus.rds")

Thymus <- JackStraw(Thymus, num.replicate = 100, dims = 30)
Thymus <- ScoreJackStraw(Thymus, dims = 1:30)
JackStrawPlot(Thymus, dims = 1:30)

saveRDS(Thymus, file = "JackStraw_Thymus.rds")
Thymus <- readRDS(file = "JackStraw_Thymus.rds")

ElbowPlot(Thymus, ndims = 30)

Thymus <- FindNeighbors(Thymus, dims = 1:30)

#louvian
Thymus <- FindClusters(Thymus, resolution = res)

clustree(Thymus, prefix = "RNA_snn_res.")

Idents(Thymus) <- Thymus$RNA_snn_res.1.2

Thymus <- BuildClusterTree(Thymus, reorder.numeric = TRUE, reorder = TRUE, dims = 1:30)
PlotClusterTree(object = Thymus)

Thymus.markers <- FindAllMarkers(Thymus)
write.csv(Thymus.markers, "Thymus_markers.csv")

Thymus <- RunTSNE(Thymus, dims = 1:30, perplexity = 50)
DimPlot(Thymus, reduction = "tsne")

Thymus <- RunUMAP(Thymus, dims = 1:30)
DimPlot(Thymus, reduction = "umap")

Thymus <- RenameIdents(Thymus, "26" = "7", "25" = "7", "24" = "7", "22" = "7", "21" = "7", "20" = "7", "19" = "7", 
                       "18" = "7", "13" = "7", "12" = "7", "11" = "7", "10" = "7", "9" = "7", "8" = "7", 
                       "17" = "14", "16" = "14", "15" = "14")

Thymus$clusters_reordered <- Thymus@active.ident

Thymus <- clustify(Thymus, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                   seurat_output = TRUE, threshold = 0.7)

Thymus_lists <- clustify(Thymus, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                         seurat_out = FALSE)

plot_cor_heatmap(Thymus_lists)

Idents(Thymus) <- Thymus$type
Thymus <- RenameIdents(Thymus, "B cell-Spleen" = "B cell-Thymus", 
                       "immature T cell-Thymus" = "T cell-Thymus",
                       "myeloid cell-Lung-CLASH!" = "monocyte-Thymus",
                       "DN4-DP transition Cd69 negative rapidly dividing thymocytes-Thymus" = "T cell-Thymus", 
                       "macrophage-Spleen" = "macrophage-Thymus")

saveRDS(Thymus, file = "Clustred_Thymus.rds")
Thymus <- readRDS(file = "Clustred_Thymus.rds")

#Bladder
#Highly variable expressed genes
Bladder <- FindVariableFeatures(Bladder, selection.method = "mean.var.plot",                                
                                dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Bladder),15)
plot1 <- VariableFeaturePlot(Bladder)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Bladder), file = "Highly_Variable_Expessed_Genes_Bladder.csv")

saveRDS(Bladder, file = "VariableGenes_Bladder.rds")
Bladder <- readRDS(file = "VariableGenes_Bladder.rds")

all.genes <- rownames(Bladder)
Bladder <- ScaleData(Bladder, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Bladder, file = "ScaledData_Bladder.rds")
Bladder <- readRDS(file = "ScaledData_Bladder.rds")

Bladder <- RunPCA(Bladder, npcs = 30)
DimPlot(Bladder, reduction = "pca")

saveRDS(Bladder, file = "PCA_Bladder.rds")
Bladder <- readRDS(file = "PCA_Bladder.rds")

Bladder <- JackStraw(Bladder, num.replicate = 100, dims = 30)
Bladder <- ScoreJackStraw(Bladder, dims = 1:30)
JackStrawPlot(Bladder, dims = 1:30)

saveRDS(Bladder, file = "JackStraw_Bladder.rds")
Bladder <- readRDS(file = "JackStraw_Bladder.rds")

ElbowPlot(Bladder, ndims = 30)

Bladder <- FindNeighbors(Bladder, dims = 1:30)

#louvian
Bladder <- FindClusters(Bladder, resolution = res)

clustree(Bladder, prefix = "RNA_snn_res.")

Idents(Bladder) <- Bladder$RNA_snn_res.0.7

Bladder <- BuildClusterTree(Bladder, reorder.numeric = TRUE, reorder = TRUE, dims = 1:30)
PlotClusterTree(object = Bladder)

Bladder.markers <- FindAllMarkers(Bladder)
write.csv(Bladder.markers, "Bladder_markers.csv")

Bladder <- RunTSNE(Bladder, dims = 1:30, perplexity = 50)
DimPlot(Bladder, reduction = "tsne")

Bladder <- RunUMAP(Bladder, dims = 1:30)
DimPlot(Bladder, reduction = "umap")

Bladder <- RenameIdents(Bladder, "3" = "T cell")

Bladder$clusters_reordered <- Bladder@active.ident

Bladder <- clustify(Bladder, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                    seurat_output = TRUE, threshold = 0.6)

Bladder_lists <- clustify(Bladder, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                         seurat_out = FALSE, threshold = 0.6)

plot_cor_heatmap(Bladder_lists)

Idents(Bladder) <- Bladder$type
Bladder <- RenameIdents(Bladder, "bladder cell-Bladder-CLASH!" = "bladder cell-Bladder",
                       "T cell-Mammary" = "T cell-Bladder",
                       "leukocyte-Bladder-CLASH!" = "monocyte-Bladder",
                       "neuroendocrine cell-Trachea" = "neuroendocrine cell-Bladder")


saveRDS(Bladder, file = "Clustred_Bladder.rds")
Bladder <- readRDS(file = "Clustred_Bladder.rds")

#Skin
#Highly variable expressed genes
Skin <- FindVariableFeatures(Skin, selection.method = "mean.var.plot",                                
                             dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0, 10))#HVG by TMS
top10 <- head(VariableFeatures(Skin),15)
plot1 <- VariableFeaturePlot(Skin)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Skin), file = "Highly_Variable_Expessed_Genes_Skin.csv")

saveRDS(Skin, file = "VariableGenes_Skin.rds")
Skin <- readRDS(file = "VariableGenes_Skin.rds")

all.genes <- rownames(Skin)
Skin <- ScaleData(Skin, features = all.genes, vars.to.regress = c("mouse.id", "age"))

saveRDS(Skin, file = "ScaledData_Skin.rds")
Skin <- readRDS(file = "ScaledData_Skin.rds")

Skin <- RunPCA(Skin, npcs = 30)
DimPlot(Skin, reduction = "pca")

saveRDS(Skin, file = "PCA_Skin.rds")
Skin <- readRDS(file = "PCA_Skin.rds")

Skin <- JackStraw(Skin, num.replicate = 100, dims = 30)
Skin <- ScoreJackStraw(Skin, dims = 1:30)
JackStrawPlot(Skin, dims = 1:30)

saveRDS(Skin, file = "JackStraw_Skin.rds")
Skin <- readRDS(file = "JackStraw_Skin.rds")

ElbowPlot(Skin, ndims = 30)

Skin <- FindNeighbors(Skin, dims = 1:30)

#louvian
Skin <- FindClusters(Skin, resolution = res)

clustree(Skin, prefix = "RNA_snn_res.")

Idents(Skin) <- Skin$RNA_snn_res.0.8

Skin <- BuildClusterTree(Skin, reorder.numeric = TRUE, reorder = TRUE, dims = 1:30)
PlotClusterTree(object = Skin)

Skin.markers <- FindAllMarkers(Skin)
write.csv(Skin.markers, "Skin_markers.csv")

Skin <- RunTSNE(Skin, dims = 1:30, perplexity = 50)
DimPlot(Skin, reduction = "tsne")

Skin <- RunUMAP(Skin, dims = 1:30)
DimPlot(Skin, reduction = "umap")

Skin <- RenameIdents(Skin, "19" = "1")

Skin$clusters_reordered <- Skin@active.ident

Skin <- clustify(Skin, ref_mat = ref_tabula_muris_facs, cluster_col = "clusters_reordered", seurat_output = TRUE, 
                 threshold = 0.6)

Skin_lists <- clustify(Skin, ref_mat = ref_tabula_muris_facs, cluster_col = "clusters_reordered", 
                          seurat_out = FALSE, threshold = 0.6)

plot_cor_heatmap(Skin_lists)

Idents(Skin) <- Skin$type
Skin <- RenameIdents(Skin, "basal cell of epidermis-Skin-CLASH!" = "basal cell of epidermis-Skin",
                        "epidermal cell-Skin-CLASH!" = "epidermal cell-Skin",
                        "T cell-Fat" = "T cell-Skin")

saveRDS(Skin, file = "Clustred_Skin.rds")
Skin <- readRDS(file = "Clustred_Skin.rds")


#Large_Intestine
#Highly variable expressed genes
Large_Intestine <- FindVariableFeatures(Large_Intestine, selection.method = "mean.var.plot",                                
                                        dispersion.cutoff = c(0, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Large_Intestine),15)
plot1 <- VariableFeaturePlot(Large_Intestine)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Large_Intestine), file = "Highly_Variable_Expessed_Genes_Large_Intestine.csv")

saveRDS(Large_Intestine, file = "VariableGenes_Large_Intestine.rds")
Large_Intestine <- readRDS(file = "VariableGenes_Large_Intestine.rds")

all.genes <- rownames(Large_Intestine)
Large_Intestine <- ScaleData(Large_Intestine, features = all.genes, vars.to.regress = "mouse.id")

saveRDS(Large_Intestine, file = "ScaledData_Large_Intestine.rds")
Large_Intestine <- readRDS(file = "ScaledData_Large_Intestine.rds")

Large_Intestine <- RunPCA(Large_Intestine, npcs = 20)
DimPlot(Large_Intestine, reduction = "pca")

saveRDS(Large_Intestine, file = "PCA_Large_Intestine.rds")
Large_Intestine <- readRDS(file = "PCA_Large_Intestine.rds")

Large_Intestine <- JackStraw(Large_Intestine, num.replicate = 100, dims = 20)
Large_Intestine <- ScoreJackStraw(Large_Intestine, dims = 1:20)
JackStrawPlot(Large_Intestine, dims = 1:20)

saveRDS(Large_Intestine, file = "JackStraw_Large_Intestine.rds")
Large_Intestine <- readRDS(file = "JackStraw_Large_Intestine.rds")

ElbowPlot(Large_Intestine, ndims = 20)

Large_Intestine <- FindNeighbors(Large_Intestine, dims = 1:20)

#louvian
Large_Intestine <- FindClusters(Large_Intestine, resolution = res)

clustree(Large_Intestine, prefix = "RNA_snn_res.")

Idents(Large_Intestine) <- Large_Intestine$RNA_snn_res.0.5

Large_Intestine <- BuildClusterTree(Large_Intestine, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = Large_Intestine)

Large_Intestine.markers <- FindAllMarkers(Large_Intestine)
write.csv(Large_Intestine.markers, "Large_Intestine_markers.csv")

Large_Intestine <- RunTSNE(Large_Intestine, dims = 1:20, perplexity = 50)
DimPlot(Large_Intestine, reduction = "tsne")

Large_Intestine <- RunUMAP(Large_Intestine, dims = 1:20)
DimPlot(Large_Intestine, reduction = "umap")

Large_Intestine <- RenameIdents(Large_Intestine, "1" = "T cell")

Large_Intestine$clusters_reordered <- Large_Intestine@active.ident

Large_Intestine <- clustify(Large_Intestine, ref_mat = ref_tabula_muris_facs, cluster_col = "clusters_reordered", 
                            seurat_output = TRUE)

Large_Intestine_lists <- clustify(Large_Intestine, ref_mat = ref_tabula_muris_facs, 
                                  cluster_col = "clusters_reordered", seurat_out = FALSE)

plot_cor_heatmap(Large_Intestine_lists)

Idents(Large_Intestine) <- Large_Intestine$type
Large_Intestine <- RenameIdents(Large_Intestine, "T cell-Fat" = "T cell-Large",
                                "Lgr5- amplifying undifferentiated cell-Large" = "Lgr5(neg) amplifying undifferentiated cell-Large",
                                "Lgr5+ undifferentiated cell (Distal)-Large" = "Lgr5(pos) undifferentiated cell (Distal)-Large")

saveRDS(Large_Intestine, file = "Clustred_Large_Intestine.rds")
Large_Intestine <- readRDS(file = "Clustred_Large_Intestine.rds")

#Trachea
#Highly variable expressed genes
Trachea <- FindVariableFeatures(Trachea, selection.method = "mean.var.plot",                                
                                dispersion.cutoff = c(0, Inf), mean.cutoff = c(0.0125, 10))#HVG by TMS
top10 <- head(VariableFeatures(Trachea),15)
plot1 <- VariableFeaturePlot(Trachea)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Trachea), file = "Highly_Variable_Expessed_Genes_Trachea.csv")

saveRDS(Trachea, file = "VariableGenes_Trachea.rds")
Trachea <- readRDS(file = "VariableGenes_Trachea.rds")

all.genes <- rownames(Trachea)
Trachea <- ScaleData(Trachea, features = all.genes, vars.to.regress = "mouse.id")

saveRDS(Trachea, file = "ScaledData_Trachea.rds")
Trachea <- readRDS(file = "ScaledData_Trachea.rds")

Trachea <- RunPCA(Trachea, npcs = 50)
DimPlot(Trachea, reduction = "pca")

saveRDS(Trachea, file = "PCA_Trachea.rds")
Trachea <- readRDS(file = "PCA_Trachea.rds")

Trachea <- JackStraw(Trachea, num.replicate = 100, dims = 50)
Trachea <- ScoreJackStraw(Trachea, dims = 1:50)
JackStrawPlot(Trachea, dims = 1:50)

saveRDS(Trachea, file = "JackStraw_Trachea.rds")
Trachea <- readRDS(file = "JackStraw_Trachea.rds")

ElbowPlot(Trachea, ndims = 50)

Trachea <- FindNeighbors(Trachea, dims = 1:50)

#louvian
Trachea <- FindClusters(Trachea, resolution = res)

clustree(Trachea, prefix = "RNA_snn_res.")

Idents(Trachea) <- Trachea$RNA_snn_res.0.7

Trachea <- BuildClusterTree(Trachea, reorder.numeric = TRUE, reorder = TRUE, dims = 1:50)
PlotClusterTree(object = Trachea)

Trachea.markers <- FindAllMarkers(Trachea)
write.csv(Trachea.markers, "Trachea_markers.csv")

Trachea <- RunTSNE(Trachea, dims = 1:50, perplexity = 50)
DimPlot(Trachea, reduction = "tsne")

Trachea <- RunUMAP(Trachea, dims = 1:50)
DimPlot(Trachea, reduction = "umap")

Trachea <- RenameIdents(Trachea, "4" = "T cell")

Trachea$clusters_reordered <- Trachea@active.ident

Trachea <- clustify(Trachea, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                    seurat_output = TRUE, threshold = 0.55)

Trachea_lists <- clustify(Trachea, ref_mat = ref_tabula_muris_drop, cluster_col = "clusters_reordered", 
                                  seurat_out = FALSE, threshold = 0.55)

plot_cor_heatmap(Trachea_lists)

Idents(Trachea) <- Trachea$type
Trachea <- RenameIdents(Trachea, "T cell-Mammary" = "T cell-Trachea",
                        "blood cell-Trachea" = "macrophage-Trachea",
                        "ciliated columnar cell of tracheobronchial tree-Lung-CLASH!" = "ciliated columnar cell of tracheobronchial tree-Trachea")

saveRDS(Trachea, file = "Clustred_Trachea.rds")
Trachea <- readRDS(file = "Clustred_Trachea.rds")

