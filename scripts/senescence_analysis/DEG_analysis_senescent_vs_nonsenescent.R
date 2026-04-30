Senis.big <- readRDS(file = "TMS_merge_TSNE_UMAP.rds")



Senescence.pos.markers_senescent_vs_nonsenescent <- FindMarkers(Senis.big, group.by = "Senescence.group",
                                                                ident.1 = "Senescence.pos", ident.2 = "Senescence.neg", 
                                                                logfc.threshold = -Inf, only.pos = FALSE)

saveRDS(Senescence.pos.markers_senescent_vs_nonsenescent, "Senescence_pos_markers_senescent_vs_nonsenescent.R")

# add a column of NAs
Senescence.pos.markers_senescent_vs_nonsenescent$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
Senescence.pos.markers_senescent_vs_nonsenescent$diffexpressed[
  Senescence.pos.markers_senescent_vs_nonsenescent$avg_logFC > 0.4 & 
    Senescence.pos.markers_senescent_vs_nonsenescent$p_val_adj < 0.05] <- "Senescence.pos"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
Senescence.pos.markers_senescent_vs_nonsenescent$diffexpressed[
  Senescence.pos.markers_senescent_vs_nonsenescent$avg_logFC < -0.4 & 
    Senescence.pos.markers_senescent_vs_nonsenescent$p_val_adj < 0.05] <- "Senescence.neg"
Senescence.pos.markers_senescent_vs_nonsenescent$delabel <- NA
Senescence.pos.markers_senescent_vs_nonsenescent$delabel[
  Senescence.pos.markers_senescent_vs_nonsenescent$diffexpressed != "NO"] <- 
  rownames(Senescence.pos.markers_senescent_vs_nonsenescent)[
    Senescence.pos.markers_senescent_vs_nonsenescent$diffexpressed != "NO"]
ggplot(Senescence.pos.markers_senescent_vs_nonsenescent, 
       aes(x = avg_logFC, y = -log10(p_val_adj), col=diffexpressed, label = delabel))+geom_point()+ 
  theme_minimal()+ geom_vline(xintercept=c(-0.4, 0.4), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")+ scale_color_manual(values=c("grey", "#F8766D", "#00BFC4")) +
  geom_text_repel() + xlab("ln(fold change)")

library(org.Mm.eg.db)

#GO and KEGG
#gene set enrichment in senescence positive cells
mm <- org.Mm.eg.db
#my.symbols <- Senescence.pos.markers_senescent_vs_nonsenescent$delabel[
#  Senescence.pos.markers_senescent_vs_nonsenescent$diffexpressed == "Senescence.pos"]

my.symbols <- Senescence.pos.markers_senescent_vs_nonsenescent$delabel

Senescence.markers_genes_for_gseGO_and_gseKEGG <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                         keytype = "SYMBOL")
df2 <- Senescence.pos.markers_senescent_vs_nonsenescent[Senescence.pos.markers_senescent_vs_nonsenescent$delabel %in% 
                            Senescence.markers_genes_for_gseGO_and_gseKEGG$SYMBOL,]
df2$id <- Senescence.markers_genes_for_gseGO_and_gseKEGG$ENTREZID
Senescence.markers_genes_list_for_gseGO_and_gseKEGG <- df2$avg_logFC
names(Senescence.markers_genes_list_for_gseGO_and_gseKEGG) <- df2$id
Senescence.markers_genes_list_for_gseGO_and_gseKEGG <- na.omit(Senescence.markers_genes_list_for_gseGO_and_gseKEGG)
Senescence.markers_genes_list_for_gseGO_and_gseKEGG <- sort(Senescence.markers_genes_list_for_gseGO_and_gseKEGG, 
                                                            decreasing = TRUE)
Senescence_gsea_kegg <- gseKEGG(Senescence.markers_genes_list_for_gseGO_and_gseKEGG, "mmu", minGSSize = 3, 
                                pAdjustMethod = "none")
Senescence_gsea_kegg <- setReadable(Senescence_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(Senescence_gsea_kegg, color = "p.adjust")
emapplot(Senescence_gsea_kegg)
cnetplot(Senescence_gsea_kegg, foldChange=Senescence.markers_genes_list_for_gseGO_and_gseKEGG)
browseKEGG(Senescence_gsea_kegg, "mmu00071")
library("pathview")
mmu00071 <- pathview(gene.data  = Senescence.markers_genes_list_for_gseGO_and_gseKEGG,
                     pathway.id = "mmu00071", species    = "mmu",
                     limit      = list(gene=max(abs(Senescence.markers_genes_list_for_gseGO_and_gseKEGG)), cpd=1))

Senescence_gsea_go <- gseGO(Senescence.markers_genes_list_for_gseGO_and_gseKEGG, OrgDb = mm, minGSSize = 3,
                            ont = "All")
Senescence_gsea_go <- setReadable(Senescence_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(Senescence_gsea_go)
emapplot(Senescence_gsea_go)
cnetplot(Senescence_gsea_go, foldChange=Senescence.markers_genes_list_for_gseGO_and_gseKEGG)
gseaplot(Senescence_gsea_go, by = "all", title = Senescence_gsea_go$Description[12], geneSetID = 1)
pmcplot(terms, 2016:2022, proportion=FALSE)


#over expression analysis of senescence cells
Senescence_over_enrich_kegg <- enrichKEGG(Senescence.markers_genes_for_gseGO_and_gseKEGG$ENTREZID, "mmu")
Senescence_over_enrich_kegg <- setReadable(Senescence_over_enrich_kegg, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_kegg)
emapplot(Senescence_over_enrich_kegg)
cnetplot(Senescence_over_enrich_kegg)
upsetplot(Senescence_over_enrich_kegg)
Senescence_over_enrich_go <- enrichGO(Senescence.markers_genes_for_gseGO_and_gseKEGG$ENTREZID, OrgDb = mm)
Senescence_over_enrich_go <- setReadable(Senescence_over_enrich_go, OrgDb = mm, "ENTREZID")
dotplot(Senescence_over_enrich_go)
emapplot(Senescence_over_enrich_go)
cnetplot(Senescence_over_enrich_go)
upsetplot(Senescence_over_enrich_go)
