Senescence_cells_list <- as.vector(WhichCells(Senis.big, idents = "Senescence.pos"))
Senescence_reactome_cells_list <- as.vector(WhichCells(Senis.big, idents = "REACTOME_CELLULAR_SENESCENCE"))
Senescence_list <- list(p21p16 = Senescence_cells_list, reactome = Senescence_reactome_cells_list)
venn.diagram(Senescence_list, "Senescence_p21_p16_vs_reactome.jpeg")
