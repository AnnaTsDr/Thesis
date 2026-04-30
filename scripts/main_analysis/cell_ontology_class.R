T_cells_2 <- subset(Senis.big, subset = cell_ontology_class == "CD4-positive, alpha-beta T cell" |
                    cell_ontology_class == "CD8-positive, alpha-beta T cell" |
                    cell_ontology_class == "immature NKT cell"|
                    cell_ontology_class == "leukocyte" | cell_ontology_class == "lymphocyte" |
                    cell_ontology_class == "mature NK T cell" | cell_ontology_class == "naive T cell" |
                    cell_ontology_class == "NK cell" | cell_ontology_class == "regulatory T cell" |
                    cell_ontology_class == "T cell")

#cell_ontology_class == "DN3 thymocyte" | cell_ontology_class == "DN4 thymocyte" | | cell_ontology_class == "thymocyte" |
#cell_ontology_class == "double negative T cell" "blood cell" | cell_ontology_class == "immature T cell" 

T_cells_3 <- subset(T_cells_2, subset = free_annotation == "leukocyte" | free_annotation == "nan" |
                      free_annotation == "T cell" | free_annotation == "CD45" | free_annotation == "CD45    T cell" |
                      free_annotation == "CD45    NK cell" | free_annotation == "NK/T" | 
                      free_annotation == "Natural Killer" | free_annotation == "CD8+ T" | 
                      free_annotation == "Natural Killer T" | free_annotation == "Proliferating T" | 
                      free_annotation == "Ly6g5b+ T" | free_annotation == "CD4+ T" |
                      free_annotation == "Regulatory T" | free_annotation == "Proliferating NK" |
                      free_annotation == "DN4" | free_annotation == "unknown" )

#free_annotation == "DN to DP transition, dividing (more DN)" | free_annotation == "Double negative thymocyte, DN3 (Cd8-, Cd4-), some express pre TCR alpha" |
#free_annotation == "DN to DP transition, dividing (some are Cd8+/ Cd4+, some undergoing VDJ recombination)" |
#free_annotation == "DN to DP transition (some are Cd8+/ Cd4+, some undergoing VDJ recombination  )" |
#free_annotation == "Double negative Thymocyte, DN4 (Cd8-, Cd4-), some undergoing VDJ recombination " |
  

T_cells_3 <- FindVariableFeatures(T_cells_3, selection.method = "mean.var.plot",                                
                                dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0.02, 2.5))#HVG by TMS
top10 <- head(VariableFeatures(T_cells_3),15)
plot1 <- VariableFeaturePlot(T_cells_3)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

T_cells_3 <- ScaleData(T_cells_3)

T_cells_3 <- RunPCA(T_cells_3, npcs = 70, ndims.print = 1:20, nfeatures.print = 5)
DimPlot(T_cells_3, reduction = "pca")

ElbowPlot(T_cells_3, ndims = 70)

T_cells_3 <- RunTSNE(T_cells_3, dims = 1:70, perplexity = 50)
DimPlot(T_cells_3, reduction = "tsne")

T_cells_3 <- RunUMAP(T_cells_3, dims = 1:70)
DimPlot(T_cells_3, reduction = "umap")

T_cells <- subset(T_cells_3, subset = cell_type == "T cell" )

T_cells <- FindVariableFeatures(T_cells, selection.method = "mean.var.plot",                                
                                dispersion.cutoff = c(0.5, Inf), mean.cutoff = c(0.02, Inf))#HVG by TMS
top10 <- head(VariableFeatures(T_cells),15)
plot1 <- VariableFeaturePlot(T_cells)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

T_cells <- ScaleData(T_cells)

T_cells <- RunPCA(T_cells, npcs = 57, ndims.print = 1:20, nfeatures.print = 5)
DimPlot(T_cells, reduction = "pca")

T_cells <- JackStraw(T_cells, num.replicate = 100, dims = 57)
T_cells <- ScoreJackStraw(T_cells, dims = 1:57)
JackStrawPlot(T_cells, dims = 1:57)

ElbowPlot(T_cells, ndims = 57)

T_cells <- RunTSNE(T_cells, dims = 1:57, perplexity = 50)
DimPlot(T_cells, reduction = "tsne")

T_cells <- RunUMAP(T_cells, dims = 1:57)
DimPlot(T_cells, reduction = "umap")

T_cells <- FindNeighbors(T_cells, dims = 1:57)

#louvian
T_cells <- FindClusters(T_cells, resolution = seq(0, 1.5, by = 0.1))

clustree(T_cells, prefix = "RNA_snn_res.")
clustree(T_cells, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)

Idents(T_cells) <- T_cells$RNA_snn_res.1.1

T_cells <- BuildClusterTree(T_cells, reorder.numeric = TRUE, reorder = TRUE, dims = 1:57)
PlotClusterTree(object = T_cells)

T_cells.markers <- FindAllMarkers(T_cells)

T_cells_ <- clustify(T_cells, ref_mat = ref_immgen, cluster_col = "RNA_snn_res.1.5", obj_out = TRUE, 
                  threshold = 0.45)

T_cells_2 <- FindVariableFeatures(T_cells_2, selection.method = "mean.var.plot",                                
                                  dispersion.cutoff = c(0.3, Inf), mean.cutoff = c(0.0125, 3))#HVG by TMS
top10 <- head(VariableFeatures(T_cells_2),15)
plot1 <- VariableFeaturePlot(T_cells_2)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

T_cells_2 <- ScaleData(T_cells_2)

T_cells_2 <- RunPCA(T_cells_2, npcs = 52, ndims.print = 1:20, nfeatures.print = 5)
DimPlot(T_cells_2, reduction = "pca")

ElbowPlot(T_cells_2, ndims = 52)

T_cells_2 <- RunTSNE(T_cells_2, dims = 1:52, perplexity = 50)
DimPlot(T_cells_2, reduction = "tsne")

T_cells_2 <- RunUMAP(T_cells_2, dims = 1:52)
DimPlot(T_cells_2, reduction = "umap")

T_cells <- subset(T_cells_2, subset = cell_type == "T cell" )
