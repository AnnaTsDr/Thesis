library(org.Mm.eg.db)

Senescence.markers_conserved_in_age_group <- FindConservedMarkers(Senis.big, ident.1 = "Senescence.pos", ident.2 = "Senescence.neg",
                                                     grouping.var = "Age_group")
Senescence.markers <- FindMarkers(Senis.big, ident.1 = "Senescence.pos", logfc.threshold = -Inf)

# add a column of NAs
Senescence.markers$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Senescence.markers$diffexpressed[Senescence.markers$avg_logFC > 0.4 & Senescence.markers$p_val_adj < 0.001] <- "Senescence.pos"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Senescence.markers$diffexpressed[Senescence.markers$avg_logFC < -0.6 & Senescence.markers$p_val_adj < 0.001] <- "Senescence.neg"
Senescence.markers$delabel <- NA
Senescence.markers$delabel[Senescence.markers$diffexpressed != "NO"] <- rownames(Senescence.markers)[Senescence.markers$diffexpressed != "NO"]
ggplot(Senescence.markers, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.6, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.001), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(Senescence.markers, "Senescence_markers_sen_pos_vs_sen_neg.csv")

write.csv(Senescence.markers_conserved_in_age_group, "Senescence_markers_conserved_in_age_group.csv")

Senescence_pos <- subset(Senis.big, subset = Senescence.group == "Senescence.pos")

Senescence_pos <- subset(Senis.big_without_defect_mice, subset = Senescence.group == "Senescence.pos")

Senescence.markers_old_vs_young <- FindMarkers(Senescence_pos, ident.1 = "Old", logfc.threshold = -Inf, 
                                               group.by = "Age_group")

# add a column of NAs
Senescence.markers_old_vs_young$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Senescence.markers_old_vs_young$diffexpressed[Senescence.markers_old_vs_young$avg_logFC > 0.4 & Senescence.markers_old_vs_young$p_val_adj < 0.05] <- "Old"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Senescence.markers_old_vs_young$diffexpressed[Senescence.markers_old_vs_young$avg_logFC < -0.4 & Senescence.markers_old_vs_young$p_val_adj < 0.05] <- "Young"
Senescence.markers_old_vs_young$delabel <- NA
Senescence.markers_old_vs_young$delabel[Senescence.markers_old_vs_young$diffexpressed != "NO"] <- rownames(Senescence.markers_old_vs_young)[Senescence.markers_old_vs_young$diffexpressed != "NO"]
ggplot(Senescence.markers_old_vs_young, aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("blue", "black", "red")) +
  geom_text_repel() + xlab("ln(fold change)")
write.csv(Senescence.markers_old_vs_young, "Senescence_markers_old_vs_young.csv")

#GO and KEGG
#gene set enrichment in senescence positive cells
Senescence.markers <- read.csv("Senescence_markers_sen_pos_vs_sen_neg.csv")
mm <- org.Mm.eg.db
my.symbols <- Senescence.markers$X
Senescence.markers_genes_for_gseGO_and_gseKEGG <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                         keytype = "SYMBOL")
df2 <- Senescence.markers[Senescence.markers$X %in% Senescence.markers_genes_for_gseGO_and_gseKEGG$SYMBOL,]
df2$id <- Senescence.markers_genes_for_gseGO_and_gseKEGG$ENTREZID
Senescence.markers_genes_list_for_gseGO_and_gseKEGG <- df2$avg_logFC
names(Senescence.markers_genes_list_for_gseGO_and_gseKEGG) <- df2$id
Senescence.markers_genes_list_for_gseGO_and_gseKEGG <- na.omit(Senescence.markers_genes_list_for_gseGO_and_gseKEGG)
Senescence.markers_genes_list_for_gseGO_and_gseKEGG <- sort(Senescence.markers_genes_list_for_gseGO_and_gseKEGG, decreasing = TRUE)
Senescence_gsea_kegg <- gseKEGG(Senescence.markers_genes_list_for_gseGO_and_gseKEGG, "mmu")
Senescence_gsea_kegg <- setReadable(Senescence_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(Senescence_gsea_kegg)
emapplot(Senescence_gsea_kegg)
cnetplot(Senescence_gsea_kegg)
Senescence_gsea_go <- gseGO(Senescence.markers_genes_list_for_gseGO_and_gseKEGG, OrgDb = mm)
Senescence_gsea_go <- setReadable(Senescence_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(Senescence_gsea_go)
emapplot(Senescence_gsea_go)
cnetplot(Senescence_gsea_go)

#gene set enrichment in old vs young senescence cells
Senescence.markers_old_vs_young <- read.csv("Senescence_markers_old_vs_young.csv")
mm <- org.Mm.eg.db
my.symbols <- Senescence.markers_old_vs_young$X
Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                         keytype = "SYMBOL")
df2 <- Senescence.markers_old_vs_young[Senescence.markers_old_vs_young$X %in% Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young$SYMBOL,]
df2$id <- Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young$ENTREZID
Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young <- df2$avg_logFC
names(Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young) <- df2$id
Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young <- na.omit(Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young)
Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young <- sort(Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young, decreasing = TRUE)
Senescence_gsea_kegg_old_vs_young <- gseKEGG(Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young, "mmu")
Senescence_gsea_kegg_old_vs_young <- setReadable(Senescence_gsea_kegg_old_vs_young, OrgDb = mm, "ENTREZID")
dotplot(Senescence_gsea_kegg_old_vs_young)
emapplot(Senescence_gsea_kegg_old_vs_young)
cnetplot(Senescence_gsea_kegg_old_vs_young)
Senescence_gsea_go_old_vs_young <- gseGO(Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young, OrgDb = mm)
Senescence_gsea_go_old_vs_young <- setReadable(Senescence_gsea_go_old_vs_young, OrgDb = mm, "ENTREZID")
dotplot(Senescence_gsea_go_old_vs_young)
emapplot(Senescence_gsea_go_old_vs_young)
cnetplot(Senescence_gsea_go_old_vs_young)

#gene set enrichment in young vs old senescence cells
Senescence.markers_old_vs_young <- read.csv("Senescence_markers_old_vs_young.csv")
mm <- org.Mm.eg.db
my.symbols <- Senescence.markers_old_vs_young$X
Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                                      keytype = "SYMBOL")
df2 <- Senescence.markers_old_vs_young[Senescence.markers_old_vs_young$X %in% Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young$SYMBOL,]
df2$id <- Senescence.markers_genes_for_gseGO_and_gseKEGG_old_vs_young$ENTREZID
Senescence.markers_genes_for_gseGO_and_gseKEGG_young_vs_old <- -(df2$avg_logFC)
names(Senescence.markers_genes_for_gseGO_and_gseKEGG_young_vs_old) <- df2$id
Senescence.markers_genes_for_gseGO_and_gseKEGG_young_vs_old <- na.omit(Senescence.markers_genes_for_gseGO_and_gseKEGG_young_vs_old)
Senescence.markers_genes_for_gseGO_and_gseKEGG_young_vs_old <- sort(Senescence.markers_genes_for_gseGO_and_gseKEGG_young_vs_old, decreasing = TRUE)
Senescence_gsea_kegg_young_vs_old <- gseKEGG(Senescence.markers_genes_for_gseGO_and_gseKEGG_young_vs_old, "mmu")
Senescence_gsea_kegg_young_vs_old <- setReadable(Senescence_gsea_kegg_young_vs_old, OrgDb = mm, "ENTREZID")
dotplot(Senescence_gsea_kegg_young_vs_old)
emapplot(Senescence_gsea_kegg_young_vs_old)
cnetplot(Senescence_gsea_kegg_young_vs_old)
Senescence_gsea_go_young_vs_old <- gseGO(Senescence.markers_genes_for_gseGO_and_gseKEGG_young_vs_old, OrgDb = mm)
Senescence_gsea_go_young_vs_old <- setReadable(Senescence_gsea_go_young_vs_old, OrgDb = mm, "ENTREZID")
dotplot(Senescence_gsea_go_young_vs_old)
emapplot(Senescence_gsea_go_young_vs_old)
cnetplot(Senescence_gsea_go_young_vs_old)

#over expression analysis of senescence cells
my.symbols <- Senescence.markers$X[Senescence.markers$diffexpressed == "Senescence.pos"]
Senescence.markers_genes_senescence_pos <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                         keytype = "SYMBOL")
Senescence_over_enrich_kegg <- enrichKEGG(Senescence.markers_genes_senescence_pos$ENTREZID, "mmu")
Senescence_over_enrich_kegg <- setReadable(Senescence_over_enrich_kegg, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_kegg)
emapplot(Senescence_over_enrich_kegg)
cnetplot(Senescence_over_enrich_kegg)
Senescence_over_enrich_go <- enrichGO(Senescence.markers_genes_senescence_pos$ENTREZID, OrgDb = mm)
Senescence_over_enrich_go <- setReadable(Senescence_over_enrich_go, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_go)
emapplot(Senescence_over_enrich_go)
cnetplot(Senescence_over_enrich_go)

#over expression analysis of old senescence cells
my.symbols <- Senescence.markers_old_vs_young$X[Senescence.markers_old_vs_young$diffexpressed == "Old"]
Senescence.markers_genes_old <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                  keytype = "SYMBOL")
Senescence_over_enrich_kegg_old <- enrichKEGG(Senescence.markers_genes_old$ENTREZID, "mmu")
Senescence_over_enrich_kegg_old <- setReadable(Senescence_over_enrich_kegg_old, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_kegg_old)
emapplot(Senescence_over_enrich_kegg_old)
cnetplot(Senescence_over_enrich_kegg_old)
Senescence_over_enrich_go_old <- enrichGO(Senescence.markers_genes_old$ENTREZID, OrgDb = mm)
Senescence_over_enrich_go_old <- setReadable(Senescence_over_enrich_go_old, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_go_old)
emapplot(Senescence_over_enrich_go_old)
cnetplot(Senescence_over_enrich_go_old)

#over expression analysis of young senescence cells
my.symbols <- Senescence.markers_old_vs_young$X[Senescence.markers_old_vs_young$diffexpressed == "Young"]
Senescence.markers_genes_young <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                  keytype = "SYMBOL")
Senescence_over_enrich_kegg_young <- enrichKEGG(Senescence.markers_genes_young$ENTREZID, "mmu")
Senescence_over_enrich_kegg_young <- setReadable(Senescence_over_enrich_kegg_young, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_kegg_young)
emapplot(Senescence_over_enrich_kegg_young)
cnetplot(Senescence_over_enrich_kegg_young)
Senescence_over_enrich_go_young <- enrichGO(Senescence.markers_genes_young$ENTREZID, OrgDb = mm)
Senescence_over_enrich_go_young <- setReadable(Senescence_over_enrich_go_young, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_go_young)
emapplot(Senescence_over_enrich_go_young)
cnetplot(Senescence_over_enrich_go_young)


Senescence_old_gene_list <- df2$avg_logFC
names(Senescence_old_gene_list) <- df2$id
Senescence_old_gene_list <- na.omit(Senescence_old_gene_list)
Senescence_old_gene_list <- sort(Senescence_old_gene_list, decreasing = TRUE)
Senescence_old_gsea_kegg <- gseKEGG(Senescence_old_gene_list, "mmu", pvalueCutoff = NONE)

dotplot(Senescence_old_gsea_kegg)
emapplot(Senescence_old_gsea_kegg)
