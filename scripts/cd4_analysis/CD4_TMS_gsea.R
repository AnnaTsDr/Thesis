CD4_T_Old <- subset(CD4_T, subset = Age_group == "Old")

CD4_T_Old.markers_cytotoxic <- FindMarkers(CD4_T_Old, ident.1 = "Cytotoxic", logfc.threshold = -Inf)
CD4_T_Old.markers_naive <- FindMarkers(CD4_T_Old, ident.1 = "Naive", logfc.threshold = -Inf)
CD4_T_Old.markers_naive_isg15 <- FindMarkers(CD4_T_Old, ident.1 = "Naive_Isg15", logfc.threshold = -Inf)
CD4_T_Old.markers_TEM <- FindMarkers(CD4_T_Old, ident.1 = "TEM", logfc.threshold = -Inf)
CD4_T_Old.markers_exhausted <- FindMarkers(CD4_T_Old, ident.1 = "Exhausted", logfc.threshold = -Inf)
CD4_T_Old.markers_aTregs <- FindMarkers(CD4_T_Old, ident.1 = "aTregs", logfc.threshold = -Inf)
CD4_T_Old.markers_rTregs <- FindMarkers(CD4_T_Old, ident.1 = "rTregs", logfc.threshold = -Inf)
write.csv(CD4_T_Old.markers_cytotoxic, "CD4_T_Old_markers_cytotoxic.csv")
write.csv(CD4_T_Old.markers_naive, "CD4_T_Old_markers_naive.csv")
write.csv(CD4_T_Old.markers_naive_isg15, "CD4_T_Old_markers_naive_isg15.csv")
write.csv(CD4_T_Old.markers_TEM, "CD4_T_Old_markers_TEM.csv")
write.csv(CD4_T_Old.markers_exhausted, "CD4_T_Old_markers_exhausted.csv")
write.csv(CD4_T_Old.markers_aTregs, "CD4_T_Old_markers_aTregs.csv")
write.csv(CD4_T_Old.markers_rTregs, "CD4_T_Old_markers_rTregs.csv")



#GO and KEGG
#gene set enrichment in CD4 CTLs old
CD4_T_Old_markers_cytotoxic <- read.csv("CD4_T_Old_markers_cytotoxic.csv")
mm <- org.Mm.eg.db
my.symbols <- CD4_T_Old_markers_cytotoxic$X
CD4_T_Old_markers_cytotoxic_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                                         keytype = "SYMBOL")
df2 <- CD4_T_Old_markers_cytotoxic[CD4_T_Old_markers_cytotoxic$X %in% CD4_T_Old_markers_cytotoxic_genes$SYMBOL,]
df2$id <- CD4_T_Old_markers_cytotoxic_genes$ENTREZID
CD4_T_Old_markers_cytotoxic_gene_list <- df2$avg_logFC
names(CD4_T_Old_markers_cytotoxic_gene_list) <- df2$id
CD4_T_Old_markers_cytotoxic_gene_list <- na.omit(CD4_T_Old_markers_cytotoxic_gene_list)
CD4_T_Old_markers_cytotoxic_gene_list <- sort(CD4_T_Old_markers_cytotoxic_gene_list, decreasing = TRUE)
CD4_T_Old_markers_cytotoxic_gsea_kegg <- gseKEGG(CD4_T_Old_markers_cytotoxic_gene_list, "mmu")
CD4_T_Old_markers_cytotoxic_gsea_kegg <- setReadable(CD4_T_Old_markers_cytotoxic_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_cytotoxic_gsea_kegg)
emapplot(CD4_T_Old_markers_cytotoxic_gsea_kegg)
cnetplot(CD4_T_Old_markers_cytotoxic_gsea_kegg)
CD4_T_Old_markers_cytotoxic_gsea_go <- gseGO(CD4_T_Old_markers_cytotoxic_gene_list, OrgDb = mm)
CD4_T_Old_markers_cytotoxic_gsea_go <- setReadable(CD4_T_Old_markers_cytotoxic_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_cytotoxic_gsea_go)
emapplot(CD4_T_Old_markers_cytotoxic_gsea_go)
cnetplot(CD4_T_Old_markers_cytotoxic_gsea_go)

#gene set enrichment in CD4 naive old
CD4_T_Old_markers_naive <- read.csv("CD4_T_Old_markers_naive.csv")
mm <- org.Mm.eg.db
my.symbols <- CD4_T_Old_markers_naive$X
CD4_T_Old_markers_naive_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                            keytype = "SYMBOL")
df2 <- CD4_T_Old_markers_naive[CD4_T_Old_markers_naive$X %in% CD4_T_Old_markers_naive_genes$SYMBOL,]
df2$id <- CD4_T_Old_markers_naive_genes$ENTREZID
CD4_T_Old_markers_naive_gene_list <- df2$avg_logFC
names(CD4_T_Old_markers_naive_gene_list) <- df2$id
CD4_T_Old_markers_naive_gene_list <- na.omit(CD4_T_Old_markers_naive_gene_list)
CD4_T_Old_markers_naive_gene_list <- sort(CD4_T_Old_markers_naive_gene_list, decreasing = TRUE)
CD4_T_Old_markers_naive_gsea_kegg <- gseKEGG(CD4_T_Old_markers_naive_gene_list, "mmu")
CD4_T_Old_markers_naive_gsea_kegg <- setReadable(CD4_T_Old_markers_naive_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_naive_gsea_kegg)
emapplot(CD4_T_Old_markers_naive_gsea_kegg)
cnetplot(CD4_T_Old_markers_naive_gsea_kegg)
CD4_T_Old_markers_naive_gsea_go <- gseGO(CD4_T_Old_markers_naive_gene_list, OrgDb = mm)
CD4_T_Old_markers_naive_gsea_go <- setReadable(CD4_T_Old_markers_naive_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_naive_gsea_go)
emapplot(CD4_T_Old_markers_naive_gsea_go)
cnetplot(CD4_T_Old_markers_naive_gsea_go)

#gene set enrichment in CD4 naive_isg15 old
CD4_T_Old_markers_naive_isg15 <- read.csv("CD4_T_Old_markers_naive_isg15.csv")
mm <- org.Mm.eg.db
my.symbols <- CD4_T_Old_markers_naive_isg15$X
CD4_T_Old_markers_naive_isg15_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                        keytype = "SYMBOL")
df2 <- CD4_T_Old_markers_naive_isg15[CD4_T_Old_markers_naive_isg15$X %in% CD4_T_Old_markers_naive_isg15_genes$SYMBOL,]
df2$id <- CD4_T_Old_markers_naive_isg15_genes$ENTREZID
CD4_T_Old_markers_naive_isg15_gene_list <- df2$avg_logFC
names(CD4_T_Old_markers_naive_isg15_gene_list) <- df2$id
CD4_T_Old_markers_naive_isg15_gene_list <- na.omit(CD4_T_Old_markers_naive_isg15_gene_list)
CD4_T_Old_markers_naive_isg15_gene_list <- sort(CD4_T_Old_markers_naive_isg15_gene_list, decreasing = TRUE)
CD4_T_Old_markers_naive_isg15_gsea_kegg <- gseKEGG(CD4_T_Old_markers_naive_isg15_gene_list, "mmu")
CD4_T_Old_markers_naive_isg15_gsea_kegg <- setReadable(CD4_T_Old_markers_naive_isg15_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_naive_isg15_gsea_kegg)
emapplot(CD4_T_Old_markers_naive_isg15_gsea_kegg)
cnetplot(CD4_T_Old_markers_naive_isg15_gsea_kegg)
CD4_T_Old_markers_naive_isg15_gsea_go <- gseGO(CD4_T_Old_markers_naive_isg15_gene_list, OrgDb = mm)
CD4_T_Old_markers_naive_isg15_gsea_go <- setReadable(CD4_T_Old_markers_naive_isg15_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_naive_isg15_gsea_go)
emapplot(CD4_T_Old_markers_naive_isg15_gsea_go)
cnetplot(CD4_T_Old_markers_naive_isg15_gsea_go)

#gene set enrichment in CD4 TEM old
CD4_T_Old_markers_TEM <- read.csv("CD4_T_Old_markers_TEM.csv")
mm <- org.Mm.eg.db
my.symbols <- CD4_T_Old_markers_TEM$X
CD4_T_Old_markers_TEM_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                              keytype = "SYMBOL")
df2 <- CD4_T_Old_markers_TEM[CD4_T_Old_markers_TEM$X %in% CD4_T_Old_markers_TEM_genes$SYMBOL,]
df2$id <- CD4_T_Old_markers_TEM_genes$ENTREZID
CD4_T_Old_markers_TEM_gene_list <- df2$avg_logFC
names(CD4_T_Old_markers_TEM_gene_list) <- df2$id
CD4_T_Old_markers_TEM_gene_list <- na.omit(CD4_T_Old_markers_TEM_gene_list)
CD4_T_Old_markers_TEM_gene_list <- sort(CD4_T_Old_markers_TEM_gene_list, decreasing = TRUE)
CD4_T_Old_markers_TEM_gsea_kegg <- gseKEGG(CD4_T_Old_markers_TEM_gene_list, "mmu")
CD4_T_Old_markers_TEM_gsea_kegg <- setReadable(CD4_T_Old_markers_TEM_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_TEM_gsea_kegg)
emapplot(CD4_T_Old_markers_TEM_gsea_kegg)
cnetplot(CD4_T_Old_markers_TEM_gsea_kegg)
CD4_T_Old_markers_TEM_gsea_go <- gseGO(CD4_T_Old_markers_TEM_gene_list, OrgDb = mm)
CD4_T_Old_markers_TEM_gsea_go <- setReadable(CD4_T_Old_markers_TEM_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_TEM_gsea_go)
emapplot(CD4_T_Old_markers_TEM_gsea_go)
cnetplot(CD4_T_Old_markers_TEM_gsea_go)

#gene set enrichment in CD4 exhausted old
CD4_T_Old_markers_exhausted <- read.csv("CD4_T_Old_markers_exhausted.csv")
mm <- org.Mm.eg.db
my.symbols <- CD4_T_Old_markers_exhausted$X
CD4_T_Old_markers_exhausted_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                      keytype = "SYMBOL")
df2 <- CD4_T_Old_markers_exhausted[CD4_T_Old_markers_exhausted$X %in% CD4_T_Old_markers_exhausted_genes$SYMBOL,]
df2$id <- CD4_T_Old_markers_exhausted_genes$ENTREZID
CD4_T_Old_markers_exhausted_gene_list <- df2$avg_logFC
names(CD4_T_Old_markers_exhausted_gene_list) <- df2$id
CD4_T_Old_markers_exhausted_gene_list <- na.omit(CD4_T_Old_markers_exhausted_gene_list)
CD4_T_Old_markers_exhausted_gene_list <- sort(CD4_T_Old_markers_exhausted_gene_list, decreasing = TRUE)
CD4_T_Old_markers_exhausted_gsea_kegg <- gseKEGG(CD4_T_Old_markers_exhausted_gene_list, "mmu")
CD4_T_Old_markers_exhausted_gsea_kegg <- setReadable(CD4_T_Old_markers_exhausted_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_exhausted_gsea_kegg)
emapplot(CD4_T_Old_markers_exhausted_gsea_kegg)
cnetplot(CD4_T_Old_markers_exhausted_gsea_kegg)
CD4_T_Old_markers_exhausted_gsea_go <- gseGO(CD4_T_Old_markers_exhausted_gene_list, OrgDb = mm)
CD4_T_Old_markers_exhausted_gsea_go <- setReadable(CD4_T_Old_markers_exhausted_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_exhausted_gsea_go)
emapplot(CD4_T_Old_markers_exhausted_gsea_go)
cnetplot(CD4_T_Old_markers_exhausted_gsea_go)

#gene set enrichment in CD4 aTregs old
CD4_T_Old_markers_aTregs <- read.csv("CD4_T_Old_markers_aTregs.csv")
mm <- org.Mm.eg.db
my.symbols <- CD4_T_Old_markers_aTregs$X
CD4_T_Old_markers_aTregs_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                            keytype = "SYMBOL")
df2 <- CD4_T_Old_markers_aTregs[CD4_T_Old_markers_aTregs$X %in% CD4_T_Old_markers_aTregs_genes$SYMBOL,]
df2$id <- CD4_T_Old_markers_aTregs_genes$ENTREZID
CD4_T_Old_markers_aTregs_gene_list <- df2$avg_logFC
names(CD4_T_Old_markers_aTregs_gene_list) <- df2$id
CD4_T_Old_markers_aTregs_gene_list <- na.omit(CD4_T_Old_markers_aTregs_gene_list)
CD4_T_Old_markers_aTregs_gene_list <- sort(CD4_T_Old_markers_aTregs_gene_list, decreasing = TRUE)
CD4_T_Old_markers_aTregs_gsea_kegg <- gseKEGG(CD4_T_Old_markers_aTregs_gene_list, "mmu")
CD4_T_Old_markers_aTregs_gsea_kegg <- setReadable(CD4_T_Old_markers_aTregs_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_aTregs_gsea_kegg)
emapplot(CD4_T_Old_markers_aTregs_gsea_kegg)
cnetplot(CD4_T_Old_markers_aTregs_gsea_kegg)
CD4_T_Old_markers_aTregs_gsea_go <- gseGO(CD4_T_Old_markers_aTregs_gene_list, OrgDb = mm)
CD4_T_Old_markers_aTregs_gsea_go <- setReadable(CD4_T_Old_markers_aTregs_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_aTregs_gsea_go)
emapplot(CD4_T_Old_markers_aTregs_gsea_go)
cnetplot(CD4_T_Old_markers_aTregs_gsea_go)

#gene set enrichment in CD4 rTregs old
CD4_T_Old_markers_rTregs <- read.csv("CD4_T_Old_markers_rTregs.csv")
mm <- org.Mm.eg.db
my.symbols <- CD4_T_Old_markers_rTregs$X
CD4_T_Old_markers_rTregs_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                         keytype = "SYMBOL")
df2 <- CD4_T_Old_markers_rTregs[CD4_T_Old_markers_rTregs$X %in% CD4_T_Old_markers_rTregs_genes$SYMBOL,]
df2$id <- CD4_T_Old_markers_rTregs_genes$ENTREZID
CD4_T_Old_markers_rTregs_gene_list <- df2$avg_logFC
names(CD4_T_Old_markers_rTregs_gene_list) <- df2$id
CD4_T_Old_markers_rTregs_gene_list <- na.omit(CD4_T_Old_markers_rTregs_gene_list)
CD4_T_Old_markers_rTregs_gene_list <- sort(CD4_T_Old_markers_rTregs_gene_list, decreasing = TRUE)
CD4_T_Old_markers_rTregs_gsea_kegg <- gseKEGG(CD4_T_Old_markers_rTregs_gene_list, "mmu")
CD4_T_Old_markers_rTregs_gsea_kegg <- setReadable(CD4_T_Old_markers_rTregs_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_rTregs_gsea_kegg)
emapplot(CD4_T_Old_markers_rTregs_gsea_kegg)
cnetplot(CD4_T_Old_markers_rTregs_gsea_kegg)
CD4_T_Old_markers_rTregs_gsea_go <- gseGO(CD4_T_Old_markers_rTregs_gene_list, OrgDb = mm)
CD4_T_Old_markers_rTregs_gsea_go <- setReadable(CD4_T_Old_markers_rTregs_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Old_markers_rTregs_gsea_go)
emapplot(CD4_T_Old_markers_rTregs_gsea_go)
cnetplot(CD4_T_Old_markers_rTregs_gsea_go)

CD4_T_Cytotoxic_vs_TEM <- read.csv("CD4_T_Cytotoxic_vs_TEM.csv")
CD4_T_Exhausted_vs_TEM <- read.csv("CD4_T_Exhausted_vs_TEM.csv")
CD4_T_aTregs_vs_rTregs <- read.csv("CD4_T_aTregs_vs_rTregs.csv")

#gene set enrichment in CD4 Cytotoxic_vs_TEM
mm <- org.Mm.eg.db
my.symbols <- CD4_T_Cytotoxic_vs_TEM$X
CD4_T_Cytotoxic_vs_TEM_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                         keytype = "SYMBOL")
df2 <- CD4_T_Cytotoxic_vs_TEM[CD4_T_Cytotoxic_vs_TEM$X %in% CD4_T_Cytotoxic_vs_TEM_genes$SYMBOL,]
df2$id <- CD4_T_Cytotoxic_vs_TEM_genes$ENTREZID
CD4_T_Cytotoxic_vs_TEM_gene_list <- df2$avg_logFC
names(CD4_T_Cytotoxic_vs_TEM_gene_list) <- df2$id
CD4_T_Cytotoxic_vs_TEM_gene_list <- na.omit(CD4_T_Cytotoxic_vs_TEM_gene_list)
CD4_T_Cytotoxic_vs_TEM_gene_list <- sort(CD4_T_Cytotoxic_vs_TEM_gene_list, decreasing = TRUE)
CD4_T_Cytotoxic_vs_TEM_gsea_kegg <- gseKEGG(CD4_T_Cytotoxic_vs_TEM_gene_list, "mmu")
CD4_T_Cytotoxic_vs_TEM_gsea_kegg <- setReadable(CD4_T_Cytotoxic_vs_TEM_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Cytotoxic_vs_TEM_gsea_kegg)
emapplot(CD4_T_Cytotoxic_vs_TEM_gsea_kegg)
cnetplot(CD4_T_Cytotoxic_vs_TEM_gsea_kegg)
CD4_T_Cytotoxic_vs_TEM_gsea_go <- gseGO(CD4_T_Cytotoxic_vs_TEM_gene_list, OrgDb = mm)
CD4_T_Cytotoxic_vs_TEM_gsea_go <- setReadable(CD4_T_Cytotoxic_vs_TEM_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Cytotoxic_vs_TEM_gsea_go)
emapplot(CD4_T_Cytotoxic_vs_TEM_gsea_go)
cnetplot(CD4_T_Cytotoxic_vs_TEM_gsea_go)

#gene set enrichment in CD4 Exhausted_vs_TEM
mm <- org.Mm.eg.db
my.symbols <- CD4_T_Exhausted_vs_TEM$X
CD4_T_Exhausted_vs_TEM_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                       keytype = "SYMBOL")
df2 <- CD4_T_Exhausted_vs_TEM[CD4_T_Exhausted_vs_TEM$X %in% CD4_T_Exhausted_vs_TEM_genes$SYMBOL,]
df2$id <- CD4_T_Exhausted_vs_TEM_genes$ENTREZID
CD4_T_Exhausted_vs_TEM_gene_list <- df2$avg_logFC
names(CD4_T_Exhausted_vs_TEM_gene_list) <- df2$id
CD4_T_Exhausted_vs_TEM_gene_list <- na.omit(CD4_T_Exhausted_vs_TEM_gene_list)
CD4_T_Exhausted_vs_TEM_gene_list <- sort(CD4_T_Exhausted_vs_TEM_gene_list, decreasing = TRUE)
CD4_T_Exhausted_vs_TEM_gsea_kegg <- gseKEGG(CD4_T_Exhausted_vs_TEM_gene_list, "mmu")
CD4_T_Exhausted_vs_TEM_gsea_kegg <- setReadable(CD4_T_Exhausted_vs_TEM_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Exhausted_vs_TEM_gsea_kegg)
emapplot(CD4_T_Exhausted_vs_TEM_gsea_kegg)
cnetplot(CD4_T_Exhausted_vs_TEM_gsea_kegg)
CD4_T_Exhausted_vs_TEM_gsea_go <- gseGO(CD4_T_Exhausted_vs_TEM_gene_list, OrgDb = mm)
CD4_T_Exhausted_vs_TEM_gsea_go <- setReadable(CD4_T_Exhausted_vs_TEM_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_Exhausted_vs_TEM_gsea_go)
emapplot(CD4_T_Exhausted_vs_TEM_gsea_go)
cnetplot(CD4_T_Exhausted_vs_TEM_gsea_go)

#gene set enrichment in CD4 aTregs_vs_rTregs
mm <- org.Mm.eg.db
my.symbols <- CD4_T_aTregs_vs_rTregs$X
CD4_T_aTregs_vs_rTregs_genes <- select(mm, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), 
                                       keytype = "SYMBOL")
df2 <- CD4_T_aTregs_vs_rTregs[CD4_T_aTregs_vs_rTregs$X %in% CD4_T_aTregs_vs_rTregs_genes$SYMBOL,]
df2$id <- CD4_T_aTregs_vs_rTregs_genes$ENTREZID
CD4_T_aTregs_vs_rTregs_gene_list <- df2$avg_logFC
names(CD4_T_aTregs_vs_rTregs_gene_list) <- df2$id
CD4_T_aTregs_vs_rTregs_gene_list <- na.omit(CD4_T_aTregs_vs_rTregs_gene_list)
CD4_T_aTregs_vs_rTregs_gene_list <- sort(CD4_T_aTregs_vs_rTregs_gene_list, decreasing = TRUE)
CD4_T_aTregs_vs_rTregs_gsea_kegg <- gseKEGG(CD4_T_aTregs_vs_rTregs_gene_list, "mmu")
CD4_T_aTregs_vs_rTregs_gsea_kegg <- setReadable(CD4_T_aTregs_vs_rTregs_gsea_kegg, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_aTregs_vs_rTregs_gsea_kegg)
emapplot(CD4_T_aTregs_vs_rTregs_gsea_kegg)
cnetplot(CD4_T_aTregs_vs_rTregs_gsea_kegg)
CD4_T_aTregs_vs_rTregs_gsea_go <- gseGO(CD4_T_aTregs_vs_rTregs_gene_list, OrgDb = mm)
CD4_T_aTregs_vs_rTregs_gsea_go <- setReadable(CD4_T_aTregs_vs_rTregs_gsea_go, OrgDb = mm, "ENTREZID")
dotplot(CD4_T_aTregs_vs_rTregs_gsea_go)
emapplot(CD4_T_aTregs_vs_rTregs_gsea_go)
cnetplot(CD4_T_aTregs_vs_rTregs_gsea_go)
