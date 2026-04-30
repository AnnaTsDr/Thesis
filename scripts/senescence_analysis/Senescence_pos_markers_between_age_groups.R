Senescence.pos <- readRDS("Senescence_subset.rds")

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
Senescence.pos <- ScaleData(Senescence.pos, features = all.genes, vars.to.regress = c("tissue", "age", "mouse.id"))

saveRDS(Senescence.pos, file = "TMS_merge_Senescence_pos_ScaledData.rds")
Senescence.pos <- readRDS(file = "TMS_merge_Senescence_pos_ScaledData.rds")

Senescence.pos <- RunPCA(Senescence.pos, ndims.print = 1:20, nfeatures.print = 5)
DimPlot(Senescence.pos, reduction = "pca")

saveRDS(Senescence.pos, file = "TMS_merge_Senescence_pos_PCA.rds")
Senescence.pos <- readRDS(file = "TMS_merge_Senescence_pos_PCA.rds")

Senescence.pos <- JackStraw(Senescence.pos, num.replicate = 100, dims = 20)
Senescence.pos <- ScoreJackStraw(Senescence.pos, dims = 1:20)
JackStrawPlot(Senescence.pos, dims = 1:20)

ElbowPlot(Senescence.pos, ndims = 20)

Senescence.pos <- RunTSNE(Senescence.pos, dims = 1:20, perplexity = 30)
DimPlot(Senescence.pos, reduction = "tsne")

Senescence.pos <- RunUMAP(Senescence.pos, dims = 1:20)
DimPlot(Senescence.pos, reduction = "umap")

Senescence.pos[["Age_group_2"]] <- plyr::mapvalues(x = Senescence.pos$age, from = c("1m", "3m", "18m", "21m", "24m", "30m"),
                                              to = c("Young", "Young", "Old", "Old", "Old", "Old"))

saveRDS(Senescence.pos, file = "TMS_merge_Senescence_pos_TSNE_UMAP.rds")
Senescence.pos <- readRDS(file = "TMS_merge_Senescence_pos_TSNE_UMAP.rds")


Senescence.pos.markers_o_vs_y <- FindMarkers(Senescence.pos, ident.1 = "Old", ident.2 = "Young", 
                                             group.by = "Age_group", logfc.threshold = -Inf)

saveRDS(Senescence.pos.markers_o_vs_y, file = "Senescence_pos_markers_old_vs_young_age_group.rds")

# add a column of NAs
Senescence.pos.markers_o_vs_y$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Senescence.pos.markers_o_vs_y$diffexpressed[Senescence.pos.markers_o_vs_y$avg_logFC > 0.4 & 
                                              Senescence.pos.markers_o_vs_y$p_val_adj < 0.05] <- "Old"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Senescence.pos.markers_o_vs_y$diffexpressed[Senescence.pos.markers_o_vs_y$avg_logFC < -0.4 & 
                                              Senescence.pos.markers_o_vs_y$p_val_adj < 0.05] <- "Young"
Senescence.pos.markers_o_vs_y$delabel <- NA
Senescence.pos.markers_o_vs_y$delabel[Senescence.pos.markers_o_vs_y$diffexpressed != "NO"] <- 
  rownames(Senescence.pos.markers_o_vs_y)[Senescence.pos.markers_o_vs_y$diffexpressed != "NO"]
ggplot(Senescence.pos.markers_o_vs_y, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+
  geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("grey", "#996035", "#65EFFF")) +
  geom_text_repel() + xlab("ln(fold change)")

Senescence.pos.markers_o_vs_s <- FindMarkers(Senescence.pos, ident.1 = "Old", ident.2 = "Supercentenarian", 
                                             group.by = "Age_group", logfc.threshold = -Inf)

saveRDS(Senescence.pos.markers_o_vs_s, file = "Senescence_pos_markers_old_vs_Supercentenarian_age_group.rds")

# add a column of NAs
Senescence.pos.markers_o_vs_s$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Senescence.pos.markers_o_vs_s$diffexpressed[Senescence.pos.markers_o_vs_s$avg_logFC > 0.4 & 
                                              Senescence.pos.markers_o_vs_s$p_val_adj < 0.05] <- "Old"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Senescence.pos.markers_o_vs_s$diffexpressed[Senescence.pos.markers_o_vs_s$avg_logFC < -0.4 & 
                                              Senescence.pos.markers_o_vs_s$p_val_adj < 0.05] <- "Supercentenarian"
Senescence.pos.markers_o_vs_s$delabel <- NA
Senescence.pos.markers_o_vs_s$delabel[Senescence.pos.markers_o_vs_s$diffexpressed != "NO"] <- 
  rownames(Senescence.pos.markers_o_vs_s)[Senescence.pos.markers_o_vs_s$diffexpressed != "NO"]
ggplot(Senescence.pos.markers_o_vs_s, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, 
                                          label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("grey", "#996035", "#331900")) +
  geom_text_repel() + xlab("ln(fold change)")

Senescence.pos.markers_o_vs_y_2 <- FindMarkers(Senescence.pos, ident.1 = "Old", ident.2 = "Young", 
                                             group.by = "Age_group_2", logfc.threshold = -Inf)

saveRDS(Senescence.pos.markers_o_vs_y_2, file = "Senescence_pos_markers_old_vs_young_age_group_2.rds")

# add a column of NAs
Senescence.pos.markers_o_vs_y_2$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Senescence.pos.markers_o_vs_y_2$diffexpressed[Senescence.pos.markers_o_vs_y_2$avg_logFC > 0.4 & 
                                                Senescence.pos.markers_o_vs_y_2$p_val_adj < 0.05] <- "Old"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Senescence.pos.markers_o_vs_y_2$diffexpressed[Senescence.pos.markers_o_vs_y_2$avg_logFC < -0.4 & 
                                              Senescence.pos.markers_o_vs_y_2$p_val_adj < 0.05] <- "Young"
Senescence.pos.markers_o_vs_y_2$delabel <- NA
Senescence.pos.markers_o_vs_y_2$delabel[Senescence.pos.markers_o_vs_y_2$diffexpressed != "NO"] <- 
  rownames(Senescence.pos.markers_o_vs_y_2)[Senescence.pos.markers_o_vs_y_2$diffexpressed != "NO"]
ggplot(Senescence.pos.markers_o_vs_y_2, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+
  geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("grey", "#996035", "#65EFFF")) +
  geom_text_repel() + xlab("ln(fold change)")

Senescence.pos.markers_y_vs_s <- FindMarkers(Senescence.pos, ident.1 = "Young", ident.2 = "Supercentenarian", 
                                             group.by = "Age_group", logfc.threshold = -Inf)

saveRDS(Senescence.pos.markers_y_vs_s, file = "Senescence_pos_markers_young_vs_Supercentenarian_age_group.rds")

# add a column of NAs
Senescence.pos.markers_y_vs_s$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Senescence.pos.markers_y_vs_s$diffexpressed[Senescence.pos.markers_y_vs_s$avg_logFC > 0.4 & 
                                              Senescence.pos.markers_y_vs_s$p_val_adj < 0.05] <- "Young"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Senescence.pos.markers_y_vs_s$diffexpressed[Senescence.pos.markers_y_vs_s$avg_logFC < -0.4 & 
                                              Senescence.pos.markers_y_vs_s$p_val_adj < 0.05] <- "Supercentenarian"
Senescence.pos.markers_y_vs_s$delabel <- NA
Senescence.pos.markers_y_vs_s$delabel[Senescence.pos.markers_y_vs_s$diffexpressed != "NO"] <- 
  rownames(Senescence.pos.markers_y_vs_s)[Senescence.pos.markers_y_vs_s$diffexpressed != "NO"]
ggplot(Senescence.pos.markers_y_vs_s, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, 
                                          label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("grey", "#65EFFF", "#331900")) +
  geom_text_repel() + xlab("ln(fold change)")

Senescence.pos.markers_s_vs_o <- FindMarkers(Senescence.pos, ident.1 = "Supercentenarian", ident.2 = "Old", 
                                             group.by = "Age_group", logfc.threshold = -Inf)

saveRDS(Senescence.pos.markers_s_vs_o, file = "Senescence_pos_markers_Supercentenarian_vs_old_age_group.rds")

# add a column of NAs
Senescence.pos.markers_s_vs_o$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Senescence.pos.markers_s_vs_o$diffexpressed[Senescence.pos.markers_s_vs_o$avg_logFC > 0.4 & 
                                              Senescence.pos.markers_s_vs_o$p_val_adj < 0.05] <- "Supercentenarian"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Senescence.pos.markers_s_vs_o$diffexpressed[Senescence.pos.markers_s_vs_o$avg_logFC < -0.4 & 
                                              Senescence.pos.markers_s_vs_o$p_val_adj < 0.05] <- "Old"
Senescence.pos.markers_s_vs_o$delabel <- NA
Senescence.pos.markers_s_vs_o$delabel[Senescence.pos.markers_s_vs_o$diffexpressed != "NO"] <- 
  rownames(Senescence.pos.markers_s_vs_o)[Senescence.pos.markers_s_vs_o$diffexpressed != "NO"]

#GO and KEGG
#gene set enrichment in senescence positive cells s_vs_o
mm <- org.Mm.eg.db
#my.symbols <- Senescence.pos.markers_s_vs_o$delabel[
#  Senescence.pos.markers_s_vs_o$diffexpressed == "Supercentenarian"]

my.symbols <- Senescence.pos.markers_s_vs_o$delabel

Senescence.markers_genes_s_vs_o_for_gseGO_and_gseKEGG <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                         keytype = "SYMBOL")
df3 <- Senescence.pos.markers_s_vs_o[Senescence.pos.markers_s_vs_o$delabel %in% 
                                       Senescence.markers_genes_s_vs_o_for_gseGO_and_gseKEGG$SYMBOL,]
df3$id <- Senescence.markers_genes_s_vs_o_for_gseGO_and_gseKEGG$ENTREZID
Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG <- df3$avg_logFC
names(Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG) <- df3$id
Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG <- 
  na.omit(Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG)
Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG <- 
  sort(Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG, decreasing = TRUE)
Senescence_s_vs_o_gsea_kegg <- gseKEGG(Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG, "mmu", 
                                       minGSSize = 3, pAdjustMethod = "none")
Senescence_s_vs_o_gsea_kegg <- setReadable(Senescence_s_vs_o_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(Senescence_s_vs_o_gsea_kegg, color = "p.adjust")
emapplot(Senescence_s_vs_o_gsea_kegg)
cnetplot(Senescence_s_vs_o_gsea_kegg, foldChange=Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG)
browseKEGG(Senescence_s_vs_o_gsea_kegg, "mmu00071")
library("pathview")
mmu00071 <- pathview(gene.data  = Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG,
                     pathway.id = "mmu00071", species    = "mmu",
                     limit      = list(gene=max(abs(Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG)), cpd=1))

Senescence_s_vs_o_gsea_go <- gseGO(Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG, OrgDb = mm, minGSSize = 3,
                            ont = "All")
Senescence_s_vs_o_gsea_go <- setReadable(Senescence_s_vs_o_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(Senescence_s_vs_o_gsea_go)
emapplot(Senescence_s_vs_o_gsea_go)
cnetplot(Senescence_s_vs_o_gsea_go, foldChange=Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG)
gseaplot(Senescence_s_vs_o_gsea_go, by = "all", title = Senescence_s_vs_o_gsea_go$Description[12], geneSetID = 1)
pmcplot(terms, 2016:2022, proportion=FALSE)


#over expression analysis of senescence cells s_vs_o
Senescence_over_enrich_kegg <- enrichKEGG(Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG$ENTREZID, "mmu")
Senescence_over_enrich_kegg <- setReadable(Senescence_over_enrich_kegg, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_kegg)
emapplot(Senescence_over_enrich_kegg)
cnetplot(Senescence_over_enrich_kegg)
upsetplot(Senescence_over_enrich_kegg)
Senescence_over_enrich_go <- enrichGO(Senescence.markers_s_vs_o_genes_list_for_gseGO_and_gseKEGG$ENTREZID, OrgDb = mm)
Senescence_over_enrich_go <- setReadable(Senescence_over_enrich_go, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_go)
emapplot(Senescence_over_enrich_go)
cnetplot(Senescence_over_enrich_go)
upsetplot(Senescence_over_enrich_go)

#GO and KEGG
#gene set enrichment in senescence positive cells o_vs_y_2
mm <- org.Mm.eg.db
#my.symbols <- Senescence.pos.markers_o_vs_y_2$delabel[
#  Senescence.pos.markers_o_vs_y_2$diffexpressed == "Old"]

my.symbols <- Senescence.pos.markers_o_vs_y_2$delabel

Senescence.markers_genes_o_vs_y_2_for_gseGO_and_gseKEGG <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                                keytype = "SYMBOL")
df4 <- Senescence.pos.markers_o_vs_y_2[Senescence.pos.markers_o_vs_y_2$delabel %in% 
                                         Senescence.markers_genes_o_vs_y_2_for_gseGO_and_gseKEGG$SYMBOL,]
df4$id <- Senescence.markers_genes_o_vs_y_2_for_gseGO_and_gseKEGG$ENTREZID
Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG <- df4$avg_logFC
names(Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG) <- df4$id
Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG <- 
  na.omit(Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG)
Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG <- 
  sort(Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG, decreasing = TRUE)
Senescence_o_vs_y_2_gsea_kegg <- gseKEGG(Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG, "mmu", 
                                       minGSSize = 3, pAdjustMethod = "none")
Senescence_o_vs_y_2_gsea_kegg <- setReadable(Senescence_o_vs_y_2_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(Senescence_o_vs_y_2_gsea_kegg, color = "p.adjust")
emapplot(Senescence_o_vs_y_2_gsea_kegg)
cnetplot(Senescence_o_vs_y_2_gsea_kegg, foldChange=Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG)
browseKEGG(Senescence_o_vs_y_2_gsea_kegg, "mmu00071")
library("pathview")
mmu00071 <- pathview(gene.data  = Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG,
                     pathway.id = "mmu00071", species    = "mmu",
                     limit      = list(gene=max(abs(Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG)), 
                                       cpd=1))

Senescence_o_vs_y_2_gsea_go <- gseGO(Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG, OrgDb = mm, 
                                   minGSSize = 3, ont = "All")
Senescence_o_vs_y_2_gsea_go <- setReadable(Senescence_o_vs_y_2_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(Senescence_o_vs_y_2_gsea_go)
emapplot(Senescence_o_vs_y_2_gsea_go)
cnetplot(Senescence_o_vs_y_2_gsea_go, foldChange=Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG)
gseaplot(Senescence_o_vs_y_2_gsea_go, by = "all", title = Senescence_o_vs_y_2_gsea_go$Description[12], geneSetID = 1)
pmcplot(terms, 2016:2022, proportion=FALSE)


#over expression analysis of senescence cells o_vs_y_2
Senescence_o_vs_y_2_over_enrich_kegg <- enrichKEGG(Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG$ENTREZID, "mmu")
Senescence_o_vs_y_2_over_enrich_kegg <- setReadable(Senescence_o_vs_y_2_over_enrich_kegg, OrgDb = mm, "ENTREZID")
dotplot(Senescence_o_vs_y_2_over_enrich_kegg)
emapplot(Senescence_o_vs_y_2_over_enrich_kegg)
cnetplot(Senescence_o_vs_y_2_over_enrich_kegg)
upsetplot(Senescence_o_vs_y_2_over_enrich_kegg)
Senescence_o_vs_y_2_over_enrich_go <- enrichGO(Senescence.markers_o_vs_y_2_genes_list_for_gseGO_and_gseKEGG$ENTREZID, OrgDb = mm)
Senescence_o_vs_y_2_over_enrich_go <- setReadable(Senescence_o_vs_y_2_over_enrich_go, OrgDb = mm, "ENTREZID")
dotplot(Senescence_o_vs_y_2_over_enrich_go)
emapplot(Senescence_o_vs_y_2_over_enrich_go)
cnetplot(Senescence_o_vs_y_2_over_enrich_go)
upsetplot(Senescence_o_vs_y_2_over_enrich_go)

