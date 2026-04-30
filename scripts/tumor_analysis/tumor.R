
#Anti_PD_L1_A_tumor
Anti_PD_L1_A_tumor <- readRDS(file = "Anti_PD_L1_A_tumor_seurat_object.rds")
Anti_PD_L1_A_tumor.metadata <- read.csv("GSE149652/Anti-PD-L1_A_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Anti_PD_L1_A_tumor.metadata) <- Anti_PD_L1_A_tumor.metadata[,1]
Anti_PD_L1_A_tumor.index <- rownames(Anti_PD_L1_A_tumor.metadata)
rownames(Anti_PD_L1_A_tumor.metadata) <- gsub("-", ".", Anti_PD_L1_A_tumor.index)

Anti_PD_L1_A_tumor <- AddMetaData(Anti_PD_L1_A_tumor, metadata = Anti_PD_L1_A_tumor.metadata)

#Anti_PD_L1_A_tumor$nCount_RNA <- Anti_PD_L1_A_tumor$n_counts
#Anti_PD_L1_A_tumor$nFeature_RNA <- Anti_PD_L1_A_tumor$n_genes

#Anti_PD_L1_B_tumor
Anti_PD_L1_B_tumor <- readRDS(file = "Anti_PD_L1_B_tumor_seurat_object.rds")
Anti_PD_L1_B_tumor.metadata <- read.csv("GSE149652/Anti-PD-L1_B_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Anti_PD_L1_B_tumor.metadata) <- Anti_PD_L1_B_tumor.metadata[,1]
Anti_PD_L1_B_tumor.index <- rownames(Anti_PD_L1_B_tumor.metadata)
rownames(Anti_PD_L1_B_tumor.metadata) <- gsub("-", ".", Anti_PD_L1_B_tumor.index)

Anti_PD_L1_B_tumor <- AddMetaData(Anti_PD_L1_B_tumor, metadata = Anti_PD_L1_B_tumor.metadata)

#Anti_PD_L1_B_tumor$nCount_RNA <- Anti_PD_L1_B_tumor$n_counts
#Anti_PD_L1_B_tumor$nFeature_RNA <- Anti_PD_L1_B_tumor$n_genes

#Anti_PD_L1_C_tumor
Anti_PD_L1_C_tumor <- readRDS(file = "Anti_PD_L1_C_tumor_seurat_object.rds")
Anti_PD_L1_C_tumor.metadata <- read.csv("GSE149652/Anti-PD-L1_C_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Anti_PD_L1_C_tumor.metadata) <- Anti_PD_L1_C_tumor.metadata[,1]
Anti_PD_L1_C_tumor.index <- rownames(Anti_PD_L1_C_tumor.metadata)
rownames(Anti_PD_L1_C_tumor.metadata) <- gsub("-", ".", Anti_PD_L1_C_tumor.index)

Anti_PD_L1_C_tumor <- AddMetaData(Anti_PD_L1_C_tumor, metadata = Anti_PD_L1_C_tumor.metadata)

#Anti_PD_L1_C_tumor$nCount_RNA <- Anti_PD_L1_C_tumor$n_counts
#Anti_PD_L1_C_tumor$nFeature_RNA <- Anti_PD_L1_C_tumor$n_genes

#Anti_PD_L1_D_tumor
Anti_PD_L1_D_tumor <- readRDS(file = "Anti_PD_L1_D_tumor_seurat_object.rds")
Anti_PD_L1_D_tumor.metadata <- read.csv("GSE149652/Anti-PD-L1_D_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Anti_PD_L1_D_tumor.metadata) <- Anti_PD_L1_D_tumor.metadata[,1]
Anti_PD_L1_D_tumor.index <- rownames(Anti_PD_L1_D_tumor.metadata)
rownames(Anti_PD_L1_D_tumor.metadata) <- gsub("-", ".", Anti_PD_L1_D_tumor.index)

Anti_PD_L1_D_tumor <- AddMetaData(Anti_PD_L1_D_tumor, metadata = Anti_PD_L1_D_tumor.metadata)

#Anti_PD_L1_D_tumor$nCount_RNA <- Anti_PD_L1_D_tumor$n_counts
#Anti_PD_L1_D_tumor$nFeature_RNA <- Anti_PD_L1_D_tumor$n_genes

#Chemo_tumor
Chemo_tumor <- readRDS(file = "Chemo_tumor_seurat_object.rds")
Chemo_tumor.metadata <- read.csv("GSE149652/Chemo_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Chemo_tumor.metadata) <- Chemo_tumor.metadata[,1]
Chemo_tumor.index <- rownames(Chemo_tumor.metadata)
rownames(Chemo_tumor.metadata) <- gsub("-", ".", Chemo_tumor.index)

Chemo_tumor <- AddMetaData(Chemo_tumor, metadata = Chemo_tumor.metadata)

#Chemo_tumor$nCount_RNA <- Chemo_tumor$n_counts
#Chemo_tumor$nFeature_RNA <- Chemo_tumor$n_genes

#Untreated_A_tumor
Untreated_A_tumor <- readRDS(file = "Untreated_A_tumor_seurat_object.rds")
Untreated_A_tumor.metadata <- read.csv("GSE149652/Untreated_A_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Untreated_A_tumor.metadata) <- Untreated_A_tumor.metadata[,1]
Untreated_A_tumor.index <- rownames(Untreated_A_tumor.metadata)
rownames(Untreated_A_tumor.metadata) <- gsub("-", ".", Untreated_A_tumor.index)

Untreated_A_tumor <- AddMetaData(Untreated_A_tumor, metadata = Untreated_A_tumor.metadata)

#Untreated_A_tumor$nCount_RNA <- Untreated_A_tumor$n_counts
#Untreated_A_tumor$nFeature_RNA <- Untreated_A_tumor$n_genes

#Untreated_B_tumor
Untreated_B_tumor <- readRDS(file = "Untreated_B_tumor_seurat_object.rds")
Untreated_B_tumor.metadata <- read.csv("GSE149652/Untreated_B_tumor_CD4_droplet_cellinfo_matrice.csv")

rownames(Untreated_B_tumor.metadata) <- Untreated_B_tumor.metadata[,1]
Untreated_B_tumor.index <- rownames(Untreated_B_tumor.metadata)
rownames(Untreated_B_tumor.metadata) <- gsub("-", ".", Untreated_B_tumor.index)

Untreated_B_tumor <- AddMetaData(Untreated_B_tumor, metadata = Untreated_B_tumor.metadata)

#Untreated_B_tumor$nCount_RNA <- Untreated_B_tumor$n_counts
#Untreated_B_tumor$nFeature_RNA <- Untreated_B_tumor$n_genes

tumor <- merge(Anti_PD_L1_A_tumor, c(Anti_PD_L1_B_tumor, Anti_PD_L1_C_tumor, Anti_PD_L1_D_tumor, Chemo_tumor, 
                                     Untreated_A_tumor, Untreated_B_tumor))

saveRDS(tumor, file = "tumor_Seurat_Object.rds")
tumor <- readRDS(file = "tumor_Seurat_Object.rds")

tumor[["percent.mt"]] <- PercentageFeatureSet(tumor, pattern = "^MT-")
VlnPlot(tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))#, "percent.mt", "percent.ribo"
FeatureScatter(tumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(tumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(tumor, feature1 = "nCount_RNA", feature2 = "n_counts")
FeatureScatter(tumor, feature1 = "nFeature_RNA", feature2 = "n_genes")
FeatureScatter(tumor, feature1 = "percent.mt", feature2 = "percent_mito")

tumor <- NormalizeData(tumor, normalization.method = "LogNormalize", scale.factor = 10000)

saveRDS(tumor, file = "tumor_Normalized.rds")
tumor <- readRDS(file = "tumor_Normalized.rds")

all.genes <- rownames(tumor)
tumor <- ScaleData(tumor, features = all.genes, vars.to.regress = c("orig.ident", "nCount_RNA", "percent.mt"))

tumor <- FindVariableFeatures(tumor, selection.method = "mean.var.plot", 
                                           dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0.0125, 3))
top10 <- head(VariableFeatures(tumor), 10)
plot1 <- VariableFeaturePlot(tumor)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

write.csv(VariableFeatures(tumor), "tumor_HVG.csv")

saveRDS(tumor, file = "tumor_HVG.rds")
tumor <- readRDS(file = "tumor_HVG.rds")

tumor <- RunPCA(tumor, npcs = 20)
DimPlot(tumor, reduction = "pca")
ElbowPlot(tumor, ndims = 20)

saveRDS(tumor, file = "tumor_PCA.rds")
tumor <- readRDS(file = "tumor_PCA.rds")

tumor <- RunTSNE(object = tumor, dims = 1:20)
TSNEPlot(object = tumor, label = TRUE, pt.size = 0.5)
FeaturePlot(tumor, features = "nCount_RNA")

tumor <- RunUMAP(object = tumor, dims = 1:20)
DimPlot(tumor, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(tumor, file = "tumor_TSNE_UMAP.rds")
tumor <- readRDS(file = "tumor_TSNE_UMAP.rds")

tumor <- FindNeighbors(tumor, dims = 1:20)
tumor <- FindClusters(tumor, resolution = seq(0,1.5, by = 0.1))

clustree(tumor, prefix = "RNA_snn_res.")
clustree(tumor, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(tumor) <- tumor$RNA_snn_res.1.5

tumor <- BuildClusterTree(tumor, reorder.numeric = TRUE, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = tumor)

tumor.markers <- FindAllMarkers(tumor)
m <- tumor.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(tumor, features = m$gene)
write.csv(tumor.markers, "tumor_markers.csv")

tumor <- RenameIdents(tumor, "2" = "1", "8" = "1", "10" = "1", "9" = "1", "16" = "11", "13" = "11",
                      "12" = "11", "7" = "4", "5" = "3", "6" = "3", "14" = "4")

tumor <- RenameIdents(tumor, "1" = "Cytotoxic", "3" = "TEM", "4" = "TCM", "11" = "aTregs", "15" = "rTregs")

tumor <- BuildClusterTree(tumor, reorder = TRUE, dims = 1:20)
PlotClusterTree(object = tumor)

tumor$cell_type_Anna <- tumor@active.ident

saveRDS(tumor, file = "tumor_Clusters.rds")
tumor <- readRDS(file = "tumor_Clusters.rds")

tumor_CTL <- subset(tumor, subset = cell_type_Anna == "Cytotoxic")

saveRDS(tumor_CTL, file = "tumor_CTL_TSNE_UMAP.rds")
tumor_CTL <- readRDS(file = "tumor_CTL_TSNE_UMAP.rds")

#monocle3
CDS <- as.cell_data_set(tumor)
CDS <- cluster_cells(cds = CDS, reduction_method = "UMAP")
plot_cells(cds = CDS, label_groups_by_cluster = FALSE, cell_size = 1, color_cells_by = "cell_type_Anna")
CDS <- learn_graph(CDS, use_partition = FALSE)
cell_ids <- colnames(CDS)[CDS$cell_type_Anna ==  "TCM"]
closest_vertex <- CDS@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_vertex <- as.matrix(closest_vertex[colnames(CDS), ])
closest_vertex <- closest_vertex[cell_ids, ]
closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
mst <- principal_graph(CDS)$UMAP
root_pr_nodes <- igraph::V(mst)$name[closest_vertex]
CDS <- order_cells(cds = CDS, root_pr_nodes = root_pr_nodes)
plot_cells(CDS, color_cells_by = "pseudotime", graph_label_size = 4)
plot_cells(CDS, color_cells_by = "cell_type_Anna", group_label_size = 5, graph_label_size = 4)
de_res <- graph_test(CDS, neighbor_graph = "principal_graph", cores = 3)

saveRDS(CDS, "tumor_CDS.rds")
