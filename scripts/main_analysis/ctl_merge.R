#CTL_CD4_Subset_agingdata <- readRDS("CTL_CD4_Subset_agingdata.rds")

CTL_CD4_Subset_agingdata <- readRDS("AgingDataCytotoxic.rds")
Supercentenarians_Cytotoxic_CD4 <- readRDS("Supercentenarians_Cytotoxic_CD4.rds")
tumor_CTL <- readRDS(file = "tumor_CTL_TSNE_UMAP.rds")

Supercentenarians_Cytotoxic_CD4$Species <- "Human"
tumor_CTL$Species <- "Human"
CTL_CD4_Subset_agingdata$Species <- "Mouse"
tumor_CTL$Condition <- "Tumor"
Supercentenarians_Cytotoxic_CD4[["Condition"]] <- plyr::mapvalues(x = Supercentenarians_Cytotoxic_CD4$V3, 
                                                                  from = c("SC", "CT"),
                                                                  to = c("Old", "Young"))
CTL_CD4_Subset_agingdata$Condition <- CTL_CD4_Subset_agingdata$Age_group
tumor_CTL$State <- "Activated"
Supercentenarians_Cytotoxic_CD4$State <- "Non-activated"
CTL_CD4_Subset_agingdata$State <- "Non-activated"
Supercentenarians_Cytotoxic_CD4$orig.data <- "Human_Aging"
tumor_CTL$orig.data <- "Bladder_Cancer"
CTL_CD4_Subset_agingdata$orig.data <- "Mouse_Aging"
Supercentenarians_Cytotoxic_CD4$treatment <- "no-tratment_aging_human"
tumor_CTL$treatment <- plyr::mapvalues(x = tumor_CTL$orig.ident, 
                                       from = c("Anti_PD_L1_A_tumor", "Anti_PD_L1_B_tumor", "Anti_PD_L1_C_tumor", 
                                                "Anti_PD_L1_D_tumor", "Chemo_tumor", "Untreated_A_tumor", 
                                                "Untreated_B_tumor"),
                                       to = c("Anti_PD_L1", "Anti_PD_L1", "Anti_PD_L1", "Anti_PD_L1", "Chemo",
                                              "Untreated", "Untreated"))
CTL_CD4_Subset_agingdata$treatment <- "no-tratment_aging_mouse"

CD4_CTL <- merge(tumor_CTL, c(Supercentenarians_Cytotoxic_CD4, CTL_CD4_Subset_agingdata))

CD4_CTL <- FindVariableFeatures(CD4_CTL, selection.method = "mean.var.plot", 
                              dispersion.cutoff = c(0.9, Inf), mean.cutoff = c(0.0125, Inf))
top10 <- head(VariableFeatures(CD4_CTL), 10)
plot1 <- VariableFeaturePlot(CD4_CTL)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

write.csv(VariableFeatures(CD4_CTL), "CD4_CTL_HVG.csv")

saveRDS(CD4_CTL, file = "CD4_CTL_HVG.rds")
CD4_CTL <- readRDS(file = "CD4_CTL_HVG.rds")

all.genes <- rownames(CD4_CTL)
CD4_CTL <- ScaleData(CD4_CTL, features = all.genes, vars.to.regress = "orig.ident")
CD4_CTL <- RunPCA(CD4_CTL, npcs = 15)
DimPlot(CD4_CTL, reduction = "pca")
ElbowPlot(CD4_CTL, ndims = 15)

saveRDS(CD4_CTL, file = "CD4_CTL_PCA.rds")
CD4_CTL <- readRDS(file = "CD4_CTL_PCA.rds")

CD4_CTL <- RunTSNE(object = CD4_CTL, dims = 1:15)
TSNEPlot(object = CD4_CTL, label = TRUE, pt.size = 0.5)
FeaturePlot(CD4_CTL, features = "nCount_RNA")

CD4_CTL <- RunUMAP(object = CD4_CTL, dims = 1:15)
DimPlot(CD4_CTL, label = TRUE, pt.size = 0.5, reduction = "umap")

saveRDS(CD4_CTL, file = "CD4_CTL_TSNE_UMAP.rds")
CD4_CTL <- readRDS(file = "CD4_CTL_TSNE_UMAP.rds")

CD4_CTL.markers_Super_vs_tomor <- FindMarkers(CD4_CTL, ident.1 = "Human_Aging", ident.2 = "Bladder_Cancer", 
                                              group.by = "orig.data", logfc.threshold = -Inf)
# add a column of NAs
CD4_CTL.markers_Super_vs_tomor$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_CTL.markers_Super_vs_tomor$diffexpressed[CD4_CTL.markers_Super_vs_tomor$avg_logFC > 0.6 & CD4_CTL.markers_Super_vs_tomor$p_val_adj < 0.5] <- "Human_Aging"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_CTL.markers_Super_vs_tomor$diffexpressed[CD4_CTL.markers_Super_vs_tomor$avg_logFC < -1 & CD4_CTL.markers_Super_vs_tomor$p_val_adj < 0.5] <- "Bladder_Cancer"
CD4_CTL.markers_Super_vs_tomor$delabel <- NA
CD4_CTL.markers_Super_vs_tomor$delabel[CD4_CTL.markers_Super_vs_tomor$diffexpressed != "NO"] <- rownames(CD4_CTL.markers_Super_vs_tomor)[CD4_CTL.markers_Super_vs_tomor$diffexpressed != "NO"]
ggplot(CD4_CTL.markers_Super_vs_tomor, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-1, 0.6), col="red") +
  geom_hline(yintercept = -log10(0.5), col="red")+ scale_color_manual(values=c("blue", "red", "grey")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(CD4_CTL.markers_Super_vs_tomor, "CD4_CTL_markers_Super_vs_tomor.csv")


CD4_CTL.markers_Super_vs_aging <- FindMarkers(CD4_CTL, ident.1 = "Human_Aging", ident.2 = "Mouse_Aging", 
                                              group.by = "orig.data", logfc.threshold = -Inf)
# add a column of NAs
CD4_CTL.markers_Super_vs_aging$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_CTL.markers_Super_vs_aging$diffexpressed[CD4_CTL.markers_Super_vs_aging$avg_logFC > 2 & CD4_CTL.markers_Super_vs_aging$p_val_adj < 0.05] <- "Human_Aging"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_CTL.markers_Super_vs_aging$diffexpressed[CD4_CTL.markers_Super_vs_aging$avg_logFC < -1 & CD4_CTL.markers_Super_vs_aging$p_val_adj < 0.05] <- "Mouse_Aging"
CD4_CTL.markers_Super_vs_aging$delabel <- NA
CD4_CTL.markers_Super_vs_aging$delabel[CD4_CTL.markers_Super_vs_aging$diffexpressed != "NO"] <- rownames(CD4_CTL.markers_Super_vs_aging)[CD4_CTL.markers_Super_vs_aging$diffexpressed != "NO"]

ggplot(CD4_CTL.markers_Super_vs_aging, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-1, 2), col="red") + 
  geom_hline(yintercept= -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "red", "grey")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(CD4_CTL.markers_Super_vs_aging, "CD4_CTL_markers_Super_vs_aging.csv")

CD4_CTL.markers_aging_vs_yuong <- FindMarkers(CD4_CTL, ident.1 = "Old", ident.2 = "Young", 
                                              group.by = "Condition", logfc.threshold = -Inf)
# add a column of NAs
CD4_CTL.markers_aging_vs_yuong$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_CTL.markers_aging_vs_yuong$diffexpressed[CD4_CTL.markers_aging_vs_yuong$avg_logFC > 0.25 & CD4_CTL.markers_aging_vs_yuong$p_val_adj < 0.05] <- "Old"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_CTL.markers_aging_vs_yuong$diffexpressed[CD4_CTL.markers_aging_vs_yuong$avg_logFC < -0.25 & CD4_CTL.markers_aging_vs_yuong$p_val_adj < 0.05] <- "Young"
CD4_CTL.markers_aging_vs_yuong$delabel <- NA
CD4_CTL.markers_aging_vs_yuong$delabel[CD4_CTL.markers_aging_vs_yuong$diffexpressed != "NO"] <- rownames(CD4_CTL.markers_aging_vs_yuong)[CD4_CTL.markers_aging_vs_yuong$diffexpressed != "NO"]

ggplot(CD4_CTL.markers_aging_vs_yuong, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.25, 0.25), col="red") + 
  geom_hline(yintercept= -log10(0.05), col="red")+ scale_color_manual(values=c("grey", "red", "blue")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(CD4_CTL.markers_aging_vs_yuong, "CD4_CTL_markers_aging_vs_yuong.csv")

CD4_CTL.markers_aging_vs_tomor <- FindMarkers(CD4_CTL, ident.1 = "Mouse_Aging", ident.2 = "Bladder_Cancer", 
                                              group.by = "orig.data", logfc.threshold = -Inf)
# add a column of NAs
CD4_CTL.markers_aging_vs_tomor$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_CTL.markers_aging_vs_tomor$diffexpressed[CD4_CTL.markers_aging_vs_tomor$avg_logFC > 1 & CD4_CTL.markers_aging_vs_tomor$p_val_adj < 0.05] <- "Mouse_Aging"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_CTL.markers_aging_vs_tomor$diffexpressed[CD4_CTL.markers_aging_vs_tomor$avg_logFC < -2.5 & CD4_CTL.markers_aging_vs_tomor$p_val_adj < 0.05] <- "Bladder_Cancer"
CD4_CTL.markers_aging_vs_tomor$delabel[CD4_CTL.markers_aging_vs_tomor$diffexpressed != "NO"] <- rownames(CD4_CTL.markers_aging_vs_tomor)[CD4_CTL.markers_aging_vs_tomor$diffexpressed != "NO"]

ggplot(CD4_CTL.markers_aging_vs_tomor, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-2.5, 1), col="red") + geom_hline(yintercept= -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "red", "grey")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(CD4_CTL.markers_aging_vs_tomor, "CD4_CTL_markers_aging_vs_tomor.csv")

CD4_CTL.markers_mouse_vs_human <- FindMarkers(CD4_CTL, ident.1 = "Mouse_Aging", ident.2 = 
                                                c("Bladder_Cancer", "Human_Aging"), 
                                              group.by = "orig.data", logfc.threshold = -Inf)
# add a column of NAs
CD4_CTL.markers_mouse_vs_human$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_CTL.markers_mouse_vs_human$diffexpressed[CD4_CTL.markers_mouse_vs_human$avg_logFC > 0.4 & CD4_CTL.markers_mouse_vs_human$p_val_adj < 0.05] <- "Mouse"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_CTL.markers_mouse_vs_human$diffexpressed[CD4_CTL.markers_mouse_vs_human$avg_logFC < -2 & CD4_CTL.markers_mouse_vs_human$p_val_adj < 0.05] <- "Human"
CD4_CTL.markers_mouse_vs_human$delabel[CD4_CTL.markers_mouse_vs_human$diffexpressed != "NO"] <- rownames(CD4_CTL.markers_mouse_vs_human)[CD4_CTL.markers_mouse_vs_human$diffexpressed != "NO"]

ggplot(CD4_CTL.markers_mouse_vs_human, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-2, 0.4), col="red") + geom_hline(yintercept= -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "red", "grey")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(CD4_CTL.markers_mouse_vs_human, "CD4_CTL_markers_mouse_vs_human.csv")

CD4_CTL.markers_Super_vs_Untreated <- FindMarkers(CD4_CTL, ident.1 = "no-tratment_aging_human", ident.2 = "Untreated", 
                                              group.by = "treatment", logfc.threshold = -Inf)
# add a column of NAs
CD4_CTL.markers_Super_vs_Untreated$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
CD4_CTL.markers_Super_vs_Untreated$diffexpressed[CD4_CTL.markers_Super_vs_Untreated$avg_logFC > 0.6 & CD4_CTL.markers_Super_vs_Untreated$p_val_adj < 0.5] <- "Human_Aging"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
CD4_CTL.markers_Super_vs_Untreated$diffexpressed[CD4_CTL.markers_Super_vs_Untreated$avg_logFC < -1 & CD4_CTL.markers_Super_vs_Untreated$p_val_adj < 0.5] <- "Untreated_Bladder_Cancer"
CD4_CTL.markers_Super_vs_Untreated$delabel <- NA
CD4_CTL.markers_Super_vs_Untreated$delabel[CD4_CTL.markers_Super_vs_Untreated$diffexpressed != "NO"] <- rownames(CD4_CTL.markers_Super_vs_Untreated)[CD4_CTL.markers_Super_vs_Untreated$diffexpressed != "NO"]
ggplot(CD4_CTL.markers_Super_vs_Untreated, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-1, 0.6), col="red") + geom_hline(yintercept = -log10(0.000005), col="red")+ scale_color_manual(values=c("blue", "red", "grey")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(CD4_CTL.markers_Super_vs_Untreated, "CD4_CTL_markers_Super_vs_Untreated.csv")


CD4_CTL <- FindNeighbors(CD4_CTL, dims = 1:15)

#louvian
CD4_CTL <- FindClusters(CD4_CTL, resolution = seq(0,1.2,by = 0.1))

saveRDS(CD4_CTL, file = "CD4_CTL_Clusters.rds")
CD4_CTL <- readRDS(file = "CD4_CTL_Clusters.rds")

Idents(CD4_CTL) <- CD4_CTL$RNA_snn_res.0.8

DimPlot(CD4_CTL, label = TRUE, pt.size = 0.5, reduction = "umap")

CD4_CTL <- BuildClusterTree(CD4_CTL, reorder.numeric = TRUE, reorder = TRUE, dims = 1:15)
PlotClusterTree(object = CD4_CTL)

saveRDS(CD4_CTL, file = "CD4_CTL_Clusters_reorder.rds")
CD4_CTL <- readRDS(file = "CD4_CTL_Clusters_reorder.rds")

CD4_CTL.markers <- FindAllMarkers(CD4_CTL)
m <- CD4_CTL.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(CD4_CTL, features = m$gene)
write.csv(CD4_CTL.markers, "CD4_CTL_markers.csv")

# subsetting rows with p-value 0:
pval_0 <- subset(CD4_CTL.markers_Super_vs_tomor, CD4_CTL.markers_Super_vs_tomor$p_val_adj == 0)
# creating a random list of numbers around 1e-300. 
pval_jitter <- round(rnorm(nrow(pval_0), mean =295, sd = 5),0) 
pval_jitter <- 1*10^-pval_jitter
# replace old p-values by new values
pval_0$p_val_adj <- pval_jitter
# creating a new data set for visualisation 
resA_1 <- subset(CD4_CTL.markers_Super_vs_tomor,CD4_CTL.markers_Super_vs_tomor$p_val_adj!=0)
resA_2 <- rbind(resA_1, pval_0)

