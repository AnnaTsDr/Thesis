#Load Senescence list
Senescence.list <- readRDS("Senescence_list.rds")

Senescence_age_group.list <- SplitObject(Senescence.list[["Senescence.pos"]], split.by = "Age_group")

#
Senescence.list <- lapply(X = Senescence.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)# 
})

features <- SelectIntegrationFeatures(object.list = Senescence.list)
Senescence.list <- lapply(X = Senescence.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs = 50)
})

ElbowPlot(Senescence.list[["Senescence.pos"]], ndims = 100)

anchors <- FindIntegrationAnchors(object.list = Senescence.list, reduction = "rpca", 
                                  dims = 1:50)
Senescence.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

Senescence.integrated <- ScaleData(Senescence.integrated, verbose = FALSE)
Senescence.integrated <- RunPCA(Senescence.integrated, npcs = 50, verbose = FALSE)
Senescence.integrated <- RunUMAP(Senescence.integrated, dims = 1:50)

DimPlot(Senescence.integrated, group.by = "orig.ident")

Idents(Senescence.integrated) <- Senescence.integrated$Senescence.group

Senescence.integrated.markers_Sen_ves_all <- FindMarkers(Senescence.integrated, ident.1 = "Senescence.pos")

#Senescence.integrated.markers_Sen_vs_all_age_group <- FindMarkers(Senescence.integrated, ident.1 = "Senescence.pos", 
 #                                                                  group.by = "Age_group")



#
Senescence_age_group.list <- lapply(X = Senescence_age_group.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)# 
})

features <- SelectIntegrationFeatures(object.list = Senescence_age_group.list)
Senescence_age_group.list <- lapply(X = Senescence_age_group.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs = 50)
})

ElbowPlot(Senescence_age_group.list[["Old"]], ndims = 50)

anchors_age_groups <- FindIntegrationAnchors(object.list = Senescence_age_group.list, reduction = "rpca", 
                                  dims = 1:50)
Senescence_age_group.integrated <- IntegrateData(anchorset = anchors_age_groups, dims = 1:50)

Senescence_age_group.integrated <- ScaleData(Senescence_age_group.integrated, verbose = FALSE)
Senescence_age_group.integrated <- RunPCA(Senescence_age_group.integrated, npcs = 50, verbose = FALSE)
Senescence_age_group.integrated <- RunUMAP(Senescence_age_group.integrated, dims = 1:50)
Senescence_age_group.integrated <- RunTSNE(Senescence_age_group.integrated, dims = 1:50, perplexity = 50)

DimPlot(Senescence_age_group.integrated, group.by = "Age_group")

Idents(Senescence_age_group.integrated) <- Senescence_age_group.integrated$Age_group

Senescence_age_group.integrated.markers <- FindAllMarkers(Senescence_age_group.integrated)

Senescence_age_group.integrated.markers_Young_vs_Old <- FindMarkers(Senescence_age_group.integrated, ident.1 = "Young", 
                                                                    ident.2 = "Old", logfc.threshold = 0)
Senescence_age_group.integrated.markers_Young_vs_Supercentenarians <- FindMarkers(Senescence_age_group.integrated, ident.1 = "Young", 
                                                                    ident.2 = "supercentenarians", logfc.threshold = 0)
Senescence_age_group.integrated.markers_Old_vs_Supercentenarians <- FindMarkers(Senescence_age_group.integrated, ident.1 = "Old", 
                                                                    ident.2 = "supercentenarians", logfc.threshold = 0)
EnhancedVolcano(Senescence_age_group.integrated.markers_Young_vs_Old, 
                lab = rownames(Senescence_age_group.integrated.markers_Young_vs_Old), x = "avg_logFC", y = "p_val_adj",
                FCcutoff = 0.25)

ggplot(Senescence_age_group.integrated.markers_Young_vs_Old, aes(x = avg_logFC, y = -log10(p_val_adj)))



#nonsenescent cells
Senescence_neg.list <- SplitObject(Senescence.list[["Senescence.neg"]], split.by = )
Senescence_neg.list <- lapply(X = Senescence.list[["Senescence.neg"]], FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)# 
})

features <- SelectIntegrationFeatures(object.list = Senescence_neg.list)
Senescence_neg.list <- lapply(X = Senescence_neg.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, npcs = 50)
})

ElbowPlot(Senescence.list[["Senescence.neg"]], ndims = 100)

anchors <- FindIntegrationAnchors(object.list = Senescence_neg.list, reduction = "rpca", 
                                  dims = 1:50)
Senescence_neg.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

Senescence_neg.integrated <- ScaleData(Senescence_neg.integrated, verbose = FALSE)
Senescence_neg.integrated <- RunPCA(Senescence_neg.integrated, npcs = 50, verbose = FALSE)
Senescence_neg.integrated <- RunUMAP(Senescence_neg.integrated, dims = 1:50)

DimPlot(Senescence_neg.integrated, group.by = "orig.ident")



#split senescence by tissue

Senescence_tissue.list <- SplitObject(Senescence.list[["Senescence.pos"]], split.by = "tissue")

Senescence_tissue_Tongue.list <- SplitObject(Senescence_tissue.list[["Tongue"]], split.by = "Age_group")
Senescence_tissue_Bladder.list <- SplitObject(Senescence_tissue.list[["Bladder"]], split.by = "Age_group")
Senescence_tissue_Heart_and_Aorta.list <- SplitObject(Senescence_tissue.list[["Heart_and_Aorta"]], split.by = "Age_group")
Senescence_tissue_Kidney.list <- SplitObject(Senescence_tissue.list[["Kidney"]], split.by = "Age_group")
Senescence_tissue_Liver.list <- SplitObject(Senescence_tissue.list[["Liver"]], split.by = "Age_group")
Senescence_tissue_Lung.list <- SplitObject(Senescence_tissue.list[["Lung"]], split.by = "Age_group")
Senescence_tissue_Marrow.list <- SplitObject(Senescence_tissue.list[["Marrow"]], split.by = "Age_group")
Senescence_tissue_Spleen.list <- SplitObject(Senescence_tissue.list[["Spleen"]], split.by = "Age_group")
Senescence_tissue_Limb_Muscle.list <- SplitObject(Senescence_tissue.list[["Limb_Muscle"]], split.by = "Age_group")

#Tongue
Senescence_tissue_Tongue.list <- lapply(X = Senescence_tissue_Tongue.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "mean.var.plot", dispersion.cutoff = c(0.5, 5), 
                            mean.cutoff = c(0.0125, 1))# 
})

top10 <- head(VariableFeatures(Senescence_tissue_Tongue.list[["Young"]]),15)
plot1 <- VariableFeaturePlot(Senescence_tissue_Tongue.list[["Young"]])
LabelPoints(plot = plot1, points = top10, repel = TRUE)


features_Tongue <- SelectIntegrationFeatures(object.list = Senescence_tissue_Tongue.list)
Senescence_tissue_Tongue.list <- lapply(X = Senescence_tissue_Tongue.list, FUN = function(x) {
  x <- ScaleData(x, features = features_Tongue, verbose = FALSE)
  x <- RunPCA(x, features = features_Tongue, verbose = FALSE, npcs = 20)
})

ElbowPlot(Senescence_tissue_Tongue.list[["Young"]], ndims = 20)

anchors_Tongue <- FindIntegrationAnchors(object.list = Senescence_tissue_Tongue.list, reduction = "rpca", 
                                  dims = 1:20)
Senescence_tissue_Tongue.integrated <- IntegrateData(anchorset = anchors_Tongue, dims = 1:20)

Senescence_tissue_Tongue.integrated <- ScaleData(Senescence_tissue_Tongue.integrated, verbose = FALSE)
Senescence_tissue_Tongue.integrated <- RunPCA(Senescence_tissue_Tongue.integrated, npcs = 20, verbose = FALSE)
Senescence_tissue_Tongue.integrated <- RunUMAP(Senescence_tissue_Tongue.integrated, dims = 1:20)

DimPlot(Senescence_tissue_Tongue.integrated, group.by = "orig.ident")

Idents(Senescence_tissue_Tongue.integrated) <- Senescence_tissue_Tongue.integrated$Age_group

Senescence_tissue_Tongue.integrated.markers_y_vs_o <- FindMarkers(Senescence_tissue_Tongue.integrated, 
                                                                  ident.1 = "Young", ident.2 = "Old")

#Bladder
Senescence_tissue_Bladder.list <- lapply(X = Senescence_tissue_Bladder.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)# 
})

features_Bladder <- SelectIntegrationFeatures(object.list = Senescence_tissue_Bladder.list)
Senescence_tissue_Bladder.list <- lapply(X = Senescence_tissue_Bladder.list, FUN = function(x) {
  x <- ScaleData(x, features = features_Bladder, verbose = FALSE)
  x <- RunPCA(x, features = features_Bladder, verbose = FALSE, npcs = 10)
})

anchors_Bladder <- FindIntegrationAnchors(object.list = Senescence_tissue_Bladder.list, reduction = "rpca", 
                                         dims = 1:10)
Senescence_tissue_Bladder.integrated <- IntegrateData(anchorset = anchors_Bladder, dims = 1:10)

Senescence_tissue_Bladder.integrated <- ScaleData(Senescence_tissue_Bladder.integrated, verbose = FALSE)
Senescence_tissue_Bladder.integrated <- RunPCA(Senescence_tissue_Bladder.integrated, npcs = 10, verbose = FALSE)
Senescence_tissue_Bladder.integrated <- RunUMAP(Senescence_tissue_Bladder.integrated, dims = 1:10)

DimPlot(Senescence_tissue_Bladder.integrated, group.by = "orig.ident")

Idents(Senescence_tissue_Bladder.integrated) <- Senescence_tissue_Bladder.integrated$Age_group

Senescence_tissue_Bladder.integrated.markers_y_vs_o <- FindMarkers(Senescence_tissue_Bladder.integrated, 
                                                                  ident.1 = "Young", ident.2 = "Old")

#Heart & Aorta





mat <- Seurat::GetAssayData(Senescence.list[["Senescence.pos"]], assay = "RNA", slot = "scale.data")
pca <- Senescence.list[["Senescence.pos"]][["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance*100

