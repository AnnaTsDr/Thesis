Senescence.pos <- readRDS(file = "Senescence_subset.rds")

#Highly variable expressed genes
Senescence.pos <- FindVariableFeatures(Senescence.pos, selection.method = "mean.var.plot",                                
                                dispersion.cutoff = c(0.8, Inf), mean.cutoff = c(0.02, Inf))#HVG by TMS
top10 <- head(VariableFeatures(Senescence.pos),15)
plot1 <- VariableFeaturePlot(Senescence.pos)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Senescence.pos), file = "TMS_merge_Senescence_pos_Highly_Variable_Expessed_Genes.csv")

saveRDS(Senescence.pos, file = "TMS_merge_Senescence_pos_VariableGenes.rds")
Senescence.pos <- readRDS(file = "TMS_merge_Senescence_pos_VariableGenes.rds")

all.genes <- rownames(Senescence.pos)
Senescence.pos <- ScaleData(Senescence.pos, features = all.genes, vars.to.regress = c("tissue", "age"))

saveRDS(Senescence.pos, file = "TMS_merge_Senescence_pos_ScaledData.rds")
Senescence.pos <- readRDS(file = "TMS_merge_Senescence_pos_ScaledData.rds")

Senescence.pos <- RunPCA(Senescence.pos, npcs = 50, ndims.print = 1:20, nfeatures.print = 5)
DimPlot(Senescence.pos, reduction = "pca")

saveRDS(Senescence.pos, file = "TMS_merge_Senescence_pos_PCA.rds")
Senescence.pos <- readRDS(file = "TMS_merge_Senescence_pos_PCA.rds")

Senescence.pos <- JackStraw(Senescence.pos, num.replicate = 100, dims = 50)
Senescence.pos <- ScoreJackStraw(Senescence.pos, dims = 1:50)
JackStrawPlot(Senescence.pos, dims = 1:50)

ElbowPlot(Senescence.pos, ndims = 50)

Senescence.pos <- RunTSNE(Senescence.pos, dims = 1:50, perplexity = 30)
DimPlot(Senescence.pos, reduction = "tsne")

Senescence.pos <- RunUMAP(Senescence.pos, dims = 1:50)
DimPlot(Senescence.pos, reduction = "umap")

saveRDS(Senescence.pos, file = "TMS_merge_Senescence_pos_TSNE_UMAP.rds")
Senescence.pos <- readRDS(file = "TMS_merge_Senescence_pos_TSNE_UMAP.rds")

Senescence.pos.markers_o_vs_y <- FindMarkers(Senescence.pos, ident.1 = "Old", ident.2 = "Young", 
                                              group.by = "Age_group", logfc.threshold = -Inf)
# add a column of NAs
Senescence.pos.markers_o_vs_y$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Senescence.pos.markers_o_vs_y$diffexpressed[Senescence.pos.markers_o_vs_y$avg_logFC > 0.4 & Senescence.pos.markers_o_vs_y$p_val_adj < 0.05] <- "Old"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Senescence.pos.markers_o_vs_y$diffexpressed[Senescence.pos.markers_o_vs_y$avg_logFC < -0.4 & Senescence.pos.markers_o_vs_y$p_val_adj < 0.05] <- "Young"
Senescence.pos.markers_o_vs_y$delabel <- NA
Senescence.pos.markers_o_vs_y$delabel[Senescence.pos.markers_o_vs_y$diffexpressed != "NO"] <- rownames(Senescence.pos.markers_o_vs_y)[Senescence.pos.markers_o_vs_y$diffexpressed != "NO"]
ggplot(Senescence.pos.markers_o_vs_y, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel()

VlnPlot(Senescence.pos, "Ccl5", group.by = "Age_group")
VlnPlot(Senescence.pos, "Ccl2", group.by = "Age_group")
VlnPlot(Senescence.pos, "Il6", group.by = "Age_group")
VlnPlot(Senescence.pos, "Il1b", group.by = "Age_group")

write.csv(Senescence.pos.markers_o_vs_y, "Senescence_pos_markers_o_vs_y.csv")

Old_sen <- subset(Senis.big, subset = Age_group == "Old")
Old_sen_vs_nonsen <- FindMarkers(Old_sen, ident.1 = "Senescence.pos", ident.2 = "Senescence.neg", 
                                 group.by = "Senescence.group", logfc.threshold = -Inf)
saveRDS(Old_sen_vs_nonsen, file = "Old_sen_vs_nonsen.rds")

# add a column of NAs
Old_sen_vs_nonsen$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Old_sen_vs_nonsen$diffexpressed[Old_sen_vs_nonsen$avg_logFC > 0.4 & Old_sen_vs_nonsen$p_val_adj < 0.05] <- "Old"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Old_sen_vs_nonsen$diffexpressed[Old_sen_vs_nonsen$avg_logFC < -0.4 & Old_sen_vs_nonsen$p_val_adj < 0.05] <- "Young"
Old_sen_vs_nonsen$delabel <- NA
Old_sen_vs_nonsen$delabel[Old_sen_vs_nonsen$diffexpressed != "NO"] <- rownames(Old_sen_vs_nonsen)[Old_sen_vs_nonsen$diffexpressed != "NO"]
ggplot(Old_sen_vs_nonsen, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel()

Old <- subset(Senis.big, subset = Age_group == "Old")
Young <- subset(Senis.big, subset = Age_group == "Young")

Old_senescence_vs_non <- FindMarkers(Old, ident.1 = "Senescence.pos", logfc.threshold = -Inf, 
                                     group.by = "Senescence.group")
Young_senescence_vs_non <- FindMarkers(Young, ident.1 = "Senescence.pos", logfc.threshold = -Inf, 
                                       group.by = "Senescence.group")

write.csv(Old_senescence_vs_non, "Old_senescence_vs_non.csv")
write.csv(Young_senescence_vs_non, "Young_senescence_vs_non.csv")
saveRDS(Old, "Old_TMS.rds")
saveRDS(Young, "Young_TMS.rds")

# add a column of NAs
Old_senescence_vs_non$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Old_senescence_vs_non$diffexpressed[Old_senescence_vs_non$avg_logFC > 0.4 & Old_senescence_vs_non$p_val_adj < 0.001] <- "Senescence.pos"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Old_senescence_vs_non$diffexpressed[Old_senescence_vs_non$avg_logFC < -0.6 & Old_senescence_vs_non$p_val_adj < 0.001] <- "Senescence.neg"
Old_senescence_vs_non$delabel <- NA
Old_senescence_vs_non$delabel[Old_senescence_vs_non$diffexpressed != "NO"] <- rownames(Old_senescence_vs_non)[Old_senescence_vs_non$diffexpressed != "NO"]
ggplot(Old_senescence_vs_non, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.6, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.001), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(Old_senescence_vs_non, "Old_senescence_vs_non_with_diffexpressed.csv")

# add a column of NAs
Young_senescence_vs_non$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Young_senescence_vs_non$diffexpressed[Young_senescence_vs_non$avg_logFC > 0.4 & Young_senescence_vs_non$p_val_adj < 0.05] <- "Senescence.pos"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Young_senescence_vs_non$diffexpressed[Young_senescence_vs_non$avg_logFC < -0.5 & Young_senescence_vs_non$p_val_adj < 0.05] <- "Senescence.neg"
Young_senescence_vs_non$delabel <- NA
Young_senescence_vs_non$delabel[Young_senescence_vs_non$diffexpressed != "NO"] <- rownames(Young_senescence_vs_non)[Young_senescence_vs_non$diffexpressed != "NO"]
ggplot(Young_senescence_vs_non, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.5, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(Young_senescence_vs_non, "Young_senescence_vs_non_with_diffexpressed.csv")
