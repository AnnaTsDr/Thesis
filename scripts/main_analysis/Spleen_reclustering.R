Spleen <- readRDS(file = "Clustred_Spleen.rds")

Spleen_T <- subset(Spleen, subset = type == "T cell-Spleen")

#Spleen
#Highly variable expressed genes
Spleen_T <- FindVariableFeatures(Spleen_T, selection.method = "mean.var.plot", dispersion.cutoff = c(0.65, Inf), 
                               mean.cutoff = c(0.0125, 2))#HVG by TMS
top10 <- head(VariableFeatures(Spleen_T),15)
plot1 <- VariableFeaturePlot(Spleen_T)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Spleen_T), file = "Highly_Variable_Expessed_Genes_Spleen_T_cells.csv")

saveRDS(Spleen_T, file = "VariableGenes_Spleen_T_cells.rds")
Spleen_T <- readRDS(file = "VariableGenes_Spleen_T_cells.rds")

all.genes <- rownames(Spleen_T)
Spleen_T <- ScaleData(Spleen_T, features = all.genes, vars.to.regress = "mouse.id")

saveRDS(Spleen_T, file = "ScaledData_Spleen_T_cells.rds")
Spleen_T <- readRDS(file = "ScaledData_Spleen_T_cells.rds")

Spleen_T <- RunPCA(Spleen_T, npcs = 30)
DimPlot(Spleen_T, reduction = "pca")

saveRDS(Spleen_T, file = "PCA_Spleen_T_cells.rds")
Spleen_T <- readRDS(file = "PCA_Spleen_T_cells.rds")

Spleen_T <- JackStraw(Spleen_T, num.replicate = 100, dims = 30)
Spleen_T <- ScoreJackStraw(Spleen_T, dims = 1:30)
JackStrawPlot(Spleen_T, dims = 1:30)

saveRDS(Spleen_T, file = "JackStraw_Spleen_T_cells.rds")
Spleen_T <- readRDS(file = "JackStraw_Spleen_T_cells.rds")

ElbowPlot(Spleen_T, ndims = 30)

Spleen_T <- RunTSNE(Spleen_T, dims = 1:30, perplexity = 30)
DimPlot(Spleen_T, reduction = "tsne")

Spleen_T <- RunUMAP(Spleen_T, dims = 1:30)
DimPlot(Spleen_T, reduction = "umap")

Spleen_T <- FindNeighbors(Spleen_T, dims = 1:30)

#louvian
Spleen_T <- FindClusters(Spleen_T, resolution = 3)

clustree(Spleen_T)

Idents(Spleen_T) <- Spleen_T$RNA_snn_res.0.4

Spleen_T <- BuildClusterTree(Spleen_T, dims = 1:30, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(Spleen_T)

saveRDS(Spleen_T, file = "Clustred_Spleen_T_cells.rds")
Spleen_T <- readRDS("Clustred_Spleen_T_cells.rds")

Spleen_T.marker_all <- FindAllMarkers(Spleen_T)
write.csv(Spleen_T.marker_all, "Spleen_marker_all_T.csv")

Spleen_T.marker_age <- list()

for (n in 1:12) {
  Spleen_T.marker_age[[n]] <- FindConservedMarkers(Spleen_T, ident.1 = n, grouping.var = "Age_group", min.cell.group = 0, min.cell.feature = 0)
  write.csv(Spleen_T.marker_age[[n]], sprintf("Spleen_marker_age_T_%s.csv", n))
}

top10 <- Spleen_T.marker_all %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(Spleen_T, features = top10$gene)
Spleen_T$my_clusters <- Spleen_T@active.ident

Spleen_T <- clustify(Spleen_T, ref_mat = ref_mouse.rnaseq, cluster_col = "RNA_snn_res.10", obj_out = TRUE, 
                   threshold = 0.4)

saveRDS(Spleen_T, "Clustred_Spleen_T_cells_type_clustify.rds")
Spleen_T <- readRDS("Clustred_Spleen_T_cells_type_clustify.rds")

Idents(Spleen_T) <- Spleen_T$type.clustify
Spleen_T <- RenameIdents(Spleen_T, "NK cells" = "NK T cell-Spleen", "T cells" = "T cell-Spleen", "B cells" = "CD3+ B cell-Spleen")

Spleen_T <- subset(Spleen_T, subset = type.clustify == "T cells")

#Highly variable expessed genes
Spleen <- FindVariableFeatures(Spleen, selection.method = "mean.var.plot", dispersion.cutoff = c(0.5, 10), mean.cutoff = c(0.0125, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Spleen),15)
plot1 <- VariableFeaturePlot(Spleen)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Spleen), file = "Highly_Variable_Expessed_Genes_Spleen_T_cells.csv")

Spleen <- ScaleData(Spleen)

Spleen <- RunPCA(Spleen, npcs = 20)
DimPlot(Spleen, reduction = "pca")

ElbowPlot(Spleen, ndims = 20)

Spleen <- RunTSNE(Spleen, dims = 1:20, perplexity = 30)
DimPlot(Spleen, reduction = "tsne")

Spleen <- RunUMAP(Spleen, dims = 1:20, perplexity = 30)
DimPlot(Spleen, reduction = "umap")

Spleen <- FindNeighbors(Spleen, dims = 1:20)

#louvian
Spleen <- FindClusters(Spleen, resolution = 1)



Spleen_TMS_CD4 <- subset(Spleen, subset = RNA_snn_res.1 == "0" | RNA_snn_res.1 == "2" | RNA_snn_res.1 == "3" | RNA_snn_res.1 == "5" | 
                           RNA_snn_res.1 == "8" | RNA_snn_res.1 == "10" | RNA_snn_res.1 == "16")

#recluster
#Highly variable expessed genes
Spleen_TMS_CD4 <- FindVariableFeatures(Spleen_TMS_CD4, selection.method = "mean.var.plot", dispersion.cutoff = c(0.5, Inf), 
                                       mean.cutoff = c(0.0125, 4))#HVG by TMS
top10 <- head(VariableFeatures(Spleen_TMS_CD4),15)
plot1 <- VariableFeaturePlot(Spleen_TMS_CD4)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Spleen_TMS_CD4), file = "Highly_Variable_Expessed_Genes_Spleen_TMS_CD4.csv")

saveRDS(Spleen_TMS_CD4, file = "VariableGenes_Spleen_TMS_CD4.rds")
Spleen_TMS_CD4 <- readRDS(file = "VariableGenes_Spleen_TMS_CD4.rds")

Spleen_TMS_CD4 <- ScaleData(Spleen_TMS_CD4)

saveRDS(Spleen_TMS_CD4, file = "ScaledData_Spleen_TMS_CD4.rds")
Spleen_TMS_CD4 <- readRDS(file = "ScaledData_Spleen_TMS_CD4.rds")

Spleen_TMS_CD4 <- RunPCA(Spleen_TMS_CD4, npcs = 20)
DimPlot(Spleen_TMS_CD4, reduction = "pca")

saveRDS(Spleen_TMS_CD4, file = "PCA_Spleen_TMS_CD4.rds")
Spleen_TMS_CD4 <- readRDS(file = "PCA_Spleen_TMS_CD4.rds")

Spleen_TMS_CD4 <- JackStraw(Spleen_TMS_CD4, num.replicate = 100, dims = 20)
Spleen_TMS_CD4 <- ScoreJackStraw(Spleen_TMS_CD4, dims = 1:20)
JackStrawPlot(Spleen_TMS_CD4, dims = 1:20)

saveRDS(Spleen_TMS_CD4, file = "JackStraw_Spleen_TMS_CD4.rds")
Spleen_TMS_CD4 <- readRDS(file = "JackStraw_Spleen_TMS_CD4.rds")

ElbowPlot(Spleen_TMS_CD4, ndims = 20)

Spleen_TMS_CD4 <- FindNeighbors(Spleen_TMS_CD4, dims = 1:20)

#louvian
Spleen_TMS_CD4 <- FindClusters(Spleen_TMS_CD4, resolution = c(0.6, 0.8, 1.2))

Spleen_TMS_CD4 <- RunTSNE(Spleen_TMS_CD4, dims = 1:20, perplexity = 30)
DimPlot(Spleen_TMS_CD4, reduction = "tsne")

Spleen_TMS_CD4 <- RunUMAP(Spleen_TMS_CD4, dims = 1:20, perplexity = 30)
DimPlot(Spleen_TMS_CD4, reduction = "umap")

Spleen_TMS_CD4.markers <- FindAllMarkers(Spleen_TMS_CD4)

Spleen_TMS_CD4 <- RenameIdents(Spleen_TMS_CD4, "0" = "Naive", "1" = "Exhausted", "2" = "Naive", "3" = "Exhausted", 
                               "4" = "aTreg", "5" = "TEM", "6" = "Exhausted",
                               "7" = "TEM", "8" = "Cytotoxic", "9" = "TEM", "10" = "TEM")
Spleen_TMS_CD4$my_clusters <- Spleen_TMS_CD4@active.ident

Spleen_TMS_CD4_lists <- clustify(Spleen_TMS_CD4, ref_mat = ref_aging, cluster_col = "RNA_snn_res.1.2", obj_out = FALSE)
plot_cor_heatmap(Spleen_TMS_CD4_lists)

Spleen_TMS_CD4 <- clustify(Spleen_TMS_CD4, ref_mat = ref_aging, cluster_col = "RNA_snn_res.1.2", obj_out = TRUE)

clu_markers <- findmarkergenes(object = Spleen, species = 'Mouse')
clu_annotation <- scCATCH(object = clu_markers$clu_markers,species = 'Mouse',tissue = 'Spleen')
