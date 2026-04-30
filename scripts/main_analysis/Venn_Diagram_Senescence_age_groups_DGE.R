Senescence.pos.markers_o_vs_y_young_genes <- as.vector(Senescence.pos.markers_o_vs_y$delabel[
  Senescence.pos.markers_o_vs_y$diffexpressed == "Young"])
Senescence.pos.markers_o_vs_s_Supercentenarian_genes <- as.vector(Senescence.pos.markers_o_vs_s$delabel[
  Senescence.pos.markers_o_vs_s$diffexpressed == "Supercentenarian"])
Senescence.pos.markers_o_vs_y_Old_genes <- as.vector(Senescence.pos.markers_o_vs_y$delabel[
  Senescence.pos.markers_o_vs_y$diffexpressed == "Old"])
Senescence.pos.markers_y_vs_s_Supercentenarian_genes <- as.vector(Senescence.pos.markers_y_vs_s$delabel[
  Senescence.pos.markers_y_vs_s$diffexpressed == "Supercentenarian"])

Venn_Diagram_young_supercentenarian_list <- list(Young = Senescence.pos.markers_o_vs_y_young_genes, 
                                                 Supercentenarian = Senescence.pos.markers_o_vs_s_Supercentenarian_genes)
Venn_Diagram_old_supercentenarian_list <- list(Old = Senescence.pos.markers_o_vs_y_Old_genes, 
                                                 Supercentenarian = Senescence.pos.markers_y_vs_s_Supercentenarian_genes)
venn.diagram(Venn_Diagram_young_supercentenarian_list, "Venn_Diagram_young_supercentenarian_Senescence_pos.jpeg", 
             fill = c("red", "green"), cex = c(3,1,3), cat.pos = c(180,180), ext.percent = 0.01)
venn.diagram(Venn_Diagram_old_supercentenarian_list, "Venn_Diagram_old_supercentenarian_Senescence_pos.jpeg",
             fill = c("red", "green"), cex = c(3,1,3), cat.pos = c(180,180), ext.percent = 0.001)
