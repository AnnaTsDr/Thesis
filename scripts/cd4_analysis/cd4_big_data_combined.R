cd4_big_data.list <- list()
v <- c("cd4t0n6", "cd4t6", "cd4t24", "Supercentenerians", "Anti_PD_L1_A_tumor", "Anti_PD_L1_B_tumor", 
       "Anti_PD_L1_C_tumor", "Anti_PD_L1_D_tumor", "Chemo_tumor", "Untreated_A_tumor", "Untreated_B_tumor", "Aging_mouse")
#for (n in length(v)){
#cd4_big_data.list[[n]] <- 
#}

Covid_data <- merge(cd4t0n6, c(cd4t6, cd4t24))

saveRDS(Covid_data, file = "Covid_data.rds")
Covid_data <- readRDS(file = "Covid_data.rds")

tumor_data <- merge(Anti_PD_L1_A_tumor, c(Anti_PD_L1_B_tumor, Anti_PD_L1_C_tumor, Anti_PD_L1_D_tumor, 
                                          Untreated_A_tumor, Untreated_B_tumor, Chemo_tumor))

saveRDS(tumor_data, file = "tumor_data.rds")
tumor_data <- readRDS(file = "tumor_data.rds")

cd4_big_data <- merge(Supercentenerians, c(tumor_data, Covid_data, agingdata))

saveRDS(cd4_big_data, file = "cd4_big_data_seurat_object.rds")
cd4_big_data <- readRDS(file = "cd4_big_data_seurat_object.rds")

cd4_big_data[["percent.mt"]] <- PercentageFeatureSet(cd4_big_data, pattern = "^MT-")
cd4_big_data[['percent.ribo']] <- PercentageFeatureSet(cd4_big_data, pattern = "^RP[SL]")
VlnPlot(cd4_big_data, features = c("nFeature_RNA", "nCount_RNA", "percent.ribo", "percent.mt"))#, "percent.mt"
FeatureScatter(cd4_big_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(cd4_big_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

cd4_big_data.list <- SplitObject(cd4_big_data, split.by = "orig.ident")



cd4_big_data.list <- lapply(X = cd4_big_data.list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "mean.var.plot", dispersion.cutoff = c(0.9, Inf), 
                            mean.cutoff = c(0.0125, 2), verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = cd4_big_data.list)
cd4_big_data.list <- lapply(X = cd4_big_data.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE, dims = 1:20)
})

reference_dataset <- which(names(cd4_big_data.list) == "Aging_mouse")

cd4_big_data.anchors <- FindIntegrationAnchors(object.list = cd4_big_data.list, reduction = "rpca", 
                                  dims = 1:20)

saveRDS(cd4_big_data.list, file = "cd4_big_data_list.rds")
cd4_big_data.list <- readRDS(file = "cd4_big_data_list.rds")
saveRDS(features, file = "cd4_big_data_features.rds")
features <- readRDS(file = "cd4_big_data_features.rds")
saveRDS(cd4_big_data.anchors, file = "cd4_big_data_Anchors_rPCA.rds")
cd4_big_data.anchors <- readRDS(file = "cd4_big_data_Anchors_rPCA.rds")

saveRDS(cd4t0n6.markers, file = "cd4t0n6_markers.rds")
saveRDS(cd4t6.markers, file = "cd4t6_markers.rds")
saveRDS(cd4t24.markers, file = "cd4t24_markers.rds")
saveRDS(Supercentenerians, file = "Supercentenerians_markers.rds")
saveRDS(Anti_PD_L1_A_tumor.markers, file = "Anti_PD_L1_A_tumor_markers.rds")
saveRDS(Anti_PD_L1_B_tumor.markers, file = "Anti_PD_L1_B_tumor_markers.rds")
saveRDS(Anti_PD_L1_C_tumor.markers, file = "Anti_PD_L1_C_tumor_markers.rds")
saveRDS(Anti_PD_L1_D_tumor.markers, file = "Anti_PD_L1_D_tumor_markers.rds")
saveRDS(Chemo_tumor.markers, file = "Chemo_tumor_markers.rds")
saveRDS(Untreated_A_tumor.markers, file = "Untreated_A_tumor_markers.rds")
saveRDS(Untreated_B_tumor.markers, file = "Untreated_B_tumor_markers.rds")

cd4_big_data.integrated <- IntegrateData(anchorset = cd4_big_data.anchors, dims = 1:20)

DefaultAssay(cd4_big_data.integrated) <- "integrated"

cd4_big_data.integrated <- ScaleData(cd4_big_data.integrated, verbose = FALSE)
cd4_big_data.integrated <- RunPCA(cd4_big_data.integrated, verbose = FALSE)
cd4_big_data.integrated <- RunUMAP(cd4_big_data.integrated, dims = 1:20)
