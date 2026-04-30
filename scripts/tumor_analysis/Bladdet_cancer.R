Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_seurat_object.rds")
Anti_PD_L1_B_tumor <- readRDS(file = "Anti_PD_L1_B_tumor_seurat_object.rds")
Anti_PD_L1_C_tumor <- readRDS(file = "Anti_PD_L1_C_tumor_seurat_object.rds")
Anti_PD_L1_D_tumor <- readRDS(file = "Anti_PD_L1_D_tumor_seurat_object.rds")
Chemo_tumor <- readRDS(file = "Chemo_tumor_seurat_object.rds")
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_seurat_object.rds")
Untreated_B_tumor <- readRDS(file = "Untreated_B_tumor_seurat_object.rds")

Anti_PD_L1_A_tumor$id <- "Anti_PD_L1_A_tumor"
Anti_PD_L1_B_tumor$id <- "Anti_PD_L1_B_tumor"
Anti_PD_L1_C_tumor$id <- "Anti_PD_L1_C_tumor"
Anti_PD_L1_D_tumor$id <- "Anti_PD_L1_D_tumor"
Chemo_tumor$id <- "Chemo_tumor"
Untreated_A_tumor$id <- "Untreated_A_tumor"
Untreated_B_tumor$id <- "Untreated_B_tumor"


Bladdet_cancer <- merge(Anti_PD_L1_A_tumor, c(Anti_PD_L1_B_tumor, Anti_PD_L1_C_tumor, Anti_PD_L1_D_tumor, Chemo_tumor,
                                              Untreated_A_tumor, Untreated_B_tumor))

saveRDS(Bladdet_cancer, file = "Bladdet_cancer.rds")
Bladdet_cancer <- readRDS(file = "Bladdet_cancer.rds")

Bladdet_cancer[["percent.mt"]] <- PercentageFeatureSet(Bladdet_cancer, pattern = "^MT-")
Bladdet_cancer[['percent.ribo']] <- PercentageFeatureSet(Bladdet_cancer, pattern = "^RP[SL]")
VlnPlot(Bladdet_cancer, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(Bladdet_cancer, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(Bladdet_cancer, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Bladdet_cancer <- NormalizeData(Bladdet_cancer, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(Bladdet_cancer, file = "Bladdet_cancer_Normalized.rds")
Bladdet_cancer <- readRDS(file = "Bladdet_cancer_Normalized.rds")

Bladdet_cancer <- subset(Bladdet_cancer, subset = CD8A == 0)
Bladdet_cancer <- subset(Bladdet_cancer, subset = CD79A == 0)
Bladdet_cancer <- subset(Bladdet_cancer, subset = LYZ == 0)
Bladdet_cancer <- subset(Bladdet_cancer, subset = CST3 == 0)

Bladdet_cancer <- FindVariableFeatures(Bladdet_cancer, selection.method = "mean.var.plot", 
                                           dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, 2))
top10 <- head(VariableFeatures(Bladdet_cancer), 10)
plot1 <- VariableFeaturePlot(Bladdet_cancer)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

write.csv(VariableFeatures(Bladdet_cancer), "Bladdet_cancer_HVG.csv")

saveRDS(Bladdet_cancer, file = "Bladdet_cancer_HVG.rds")
Bladdet_cancer <- readRDS(file = "Bladdet_cancer_HVG.rds")

all.genes <- rownames(Bladdet_cancer)
Bladdet_cancer <- ScaleData(Bladdet_cancer, features = all.genes)
Bladdet_cancer <- RunPCA(Bladdet_cancer, npcs = 15)
DimPlot(Bladdet_cancer, reduction = "pca")
ElbowPlot(Bladdet_cancer, ndims = 15)

saveRDS(Bladdet_cancer, file = "Bladdet_cancer_PCA.rds")
Bladdet_cancer <- readRDS(file = "Bladdet_cancer_PCA.rds")

Bladdet_cancer <- RunTSNE(object = Bladdet_cancer, dims = 1:15)
TSNEPlot(object = Bladdet_cancer, label = TRUE, pt.size = 0.5)
FeaturePlot(Bladdet_cancer, features = "nCount_RNA")

Bladdet_cancer <- RunUMAP(object = Bladdet_cancer, dims = 1:15)
DimPlot(Bladdet_cancer, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(Bladdet_cancer, file = "Bladdet_cancer_TSNE_UMAP.rds")
Bladdet_cancer <- readRDS(file = "Bladdet_cancer_TSNE_UMAP.rds")

Bladdet_cancer <- FindNeighbors(Bladdet_cancer, dims = 1:15)
res <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2)
Bladdet_cancer <- FindClusters(Bladdet_cancer, resolution = res)

clustree(Bladdet_cancer, prefix = "RNA_snn_res.")

saveRDS(Bladdet_cancer, "Bladdet_cancer_clustered.rds")
Bladdet_cancer <- readRDS(file = "Bladdet_cancer_clustered.rds")

Idents(Bladdet_cancer) <- Bladdet_cancer$RNA_snn_res.1.2

Bladdet_cancer <- BuildClusterTree(object = Bladdet_cancer, dims = 1:15, reorder = TRUE, reorder.numeric = TRUE)
PlotClusterTree(object = Bladdet_cancer)

Bladdet_cancer.markers <- FindAllMarkers(Bladdet_cancer)
m <- Bladdet_cancer.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Bladdet_cancer, features = m$gene)
write.csv(Bladdet_cancer.markers, "Bladdet_cancer_all_markers_before_names.csv")

Bladdet_cancer.table_idents <- as.data.frame(table(Idents(Bladdet_cancer)))
write.csv(Bladdet_cancer.table_idents, "Bladdet_cancer_table_idents_before_names.csv")

Bladdet_cancer <- RenameIdents(Bladdet_cancer, "10" = "8", "9" = "8")