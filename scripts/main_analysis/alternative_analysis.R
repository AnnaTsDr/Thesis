Senis.big <- readRDS("TMS_merge.rds")

#Highly variable expressed genes
Senis.big <- FindVariableFeatures(Senis.big, selection.method = "mean.var.plot", dispersion.cutoff = c(0, Inf), 
                                  mean.cutoff = c(0.0125, Inf))
top10 <- head(VariableFeatures(Senis.big),15)
plot1 <- VariableFeaturePlot(Senis.big)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#Save variable genes in csv
write.csv(VariableFeatures(Senis.big), file = "Highly_Variable_Expessed_Genes.csv")

saveRDS(Senis.big, file = "VariableGenes.rds")
Senis.big <- readRDS(file = "VariableGenes.rds")

t.genes <- c("Ptprc", "Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", "Cd8b1", 
             "Cd40lg", "Cd28", "Cd27", "Cd44", "Il7r", "Gzmb", 
             "Itga1", "Ifng", "Ccl3", "Prf1", "Il17a", "Il17f", "Cxcl10", "Stat5a", "Stat5b", "Tfrc", 
             "Slc3a2", "Ebi3", "H2-Ea", "H2-Eb1", "Il27", "Itgae", "Tgfb1", "Il12a", "Smad3", "Il12b",
             "Gata3", "Rorc", "Tbx21", "Crtam", "Il4", "Tnf", "Il6", "Il13", "Il9", "Il21", "Il22",
             "Il25", "Il10", "Cd26", "Fas", "Cd38", "Cd69", "Mki67", "B3gat1", "Fcrl6", "Gzma", "Il23a",
             "Ccr1", "Ccr3", "Ccr4", "Ccr6", "Ccr8", "Ccr10", "Cxcr4", "Cxcr5", "Entpd1", "Il17rb",
             "Pdcd1", "Cd94", "Glb1", "Klrd1", "Slamf1", "Cd160", "Slamf6", "Itga2", "Ncam1", "Ncam2",
             "Il15ra", "Zbtb16", "Pdgfra", "Il1a", "Il1b", "Tgfb2", "Tgfb3", "Ahr", "Icos", "Rora", "Stat3", 
             "Runx1", "Batf", "Irf4", "Il24", "Ccl20", "Il6ra", "Il13ra1", "Il21r", "Il23r", "Il18", "Tnfsf11",
             "Havcr2", "Fasl", "Runx3", "Stat4", "Il12rb2", "Il18ra", "Il27ra", "Notch3", "Lta", "Il31", "Il33",
             "Tslp", "Stat5a", "Stat5b", "Stat6", "Gfi1", "Cd1d1",
             "Ikzf2", "Cd74", "Tnfrsf4", "Cd81", "Tnfrsf9", "Tnfrsf18", "Cst7", "Foxp3", "AW112010", "Tpt1", "Ms4a4b",
             "Cd83", "Ctla4", "Maf", "Ifit1", "Ifit3", "Isg15", "Rtp4", "Igfbp4", "Dapl1", "Lef1", "Satb1", "Ccr7",
             "Npc2", "Sell", "Gm8730", "Gm9493", "Gnb2l1", "Eef1a1", "B2m", "H2-D1", "H2-K1", "Malat1", "Prr13",
             "Cd82", "S100a10", "Ifi27l2a", "Srgn", "Izumo1r", "S100a4", "S100a11", "S100a6", "Cxcr3", "Il2rb", 
             "Tnfsf8", "Sostdc1", "Gpm6b", "Lag3", "Tbc1d4", "Ly6a", "Eea1", "Cd200", "Itgb1", "Angptl2", "Ccl5",
             "Nkg7", "Gzmk", "Ccl4", "Ctla2a", "Eomes", "Xcl1", "Ier3", "Ccr5", "Il2ra")

sen.genes <- c("Cdkn2a", "Cdkn2b", "Cdkn1a", "Mtor", "E2f2", "Lmnb1", "Tnf", "Itgax", "Krt15", "Krt18", "Sfn", "Lgals3", 
               "Igfbp2", "Ly6d", "Tnfrsf12a", "Glb1", "Cdkn2c")

t_hvg <- intersect(t.genes, VariableFeatures(Senis.big))

Senis.big <- ScaleData(Senis.big, fearutres = t.genes)
Senis.big <- ScaleData(Senis.big, fearutres = sen.genes)

Senis.big <- RunPCA(Senis.big, fearutres = t.genes)
ElbowPlot(Senis.big)

Senis.big <- RunTSNE(Senis.big, dims = 1:50, perplexity = 50)
DimPlot(Senis.big, reduction = "tsne")

Senis.big <- FindNeighbors(Senis.big, dims = 1:50)

#louvian
Senis.big <- FindClusters(Senis.big, resolution = seq(0, 1.2, by = 0.1))

clustree(Senis.big, prefix = "RNA_snn_res.")
clustree(Senis.big, node_colour = "sc3_stability") + scale_colour_viridis_c(option = 'plasma', begin = 0.3)
clustree(Senis.big, prefix = "RNA_snn_res.", node_colour_aggr = "median", node_colour = "Ptprc", exprs = 'scale.data') + 
  scale_colour_viridis_c(option = 'plasma', begin = 0.3)
