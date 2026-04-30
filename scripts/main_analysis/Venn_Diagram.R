library(VennDiagram)

Senescence.pos.markers_o_vs_y <- read.csv(Senescence.pos.markers_o_vs_y, "Senescence_pos_markers_o_vs_y.csv")

senescence_list_old_and_young <- list()
senescence_list_old_and_young$A <- rownames(Old_senescence_vs_non)[Old_senescence_vs_non$diffexpressed == "Senescence.pos"]
senescence_list_old_and_young$B <- rownames(Young_senescence_vs_non)[Young_senescence_vs_non$diffexpressed == "Senescence.pos"]
senescence_list_old_and_young$C <- rownames(Senescence.pos.markers_o_vs_y)[Senescence.pos.markers_o_vs_y$diffexpressed != "NO"]
v0 <- venn.diagram(senescence_list_old_and_young, NULL, 
             category.names = c("Old_sen", "Young_sen", "sen"), fill = c("red", "blue", "yellow"))#, cex=c(1,1,1,)
grid.draw(v0)

overlaps <- calculate.overlap(senescence_list_old_and_young)
overlaps <- rev(overlaps)
indx <- as.numeric(substr(names(overlaps),2,2))
for (i in 1:length(overlaps)){
  v0[[i+4]]$label <- paste(overlaps[[i]], collapse = "\n") 
}
grid.newpage()
grid.draw(v0)
write.csv(CD4_Cytotoxic, "CD4_Cytotoxic_TMS_DEG.csv")


agingdata.markers_Cytotoxic <- read.csv("agingdata_markers_Cytotoxic.csv")
CD4_Cytotoxic <- read.csv("CD4_Cytotoxic_TMS_DEG.csv")

Cytotoxic_TMS_vs_aging_mice <- list()
Cytotoxic_TMS_vs_aging_mice$A <- rownames(agingdata.markers_Cytotoxic)[agingdata.markers_Cytotoxic$diffexpressed == "Cytotoxic"]
Cytotoxic_TMS_vs_aging_mice$B <- rownames(CD4_T_Cytotoxic_vs_other )[CD4_T_Cytotoxic_vs_other $diffexpressed == "Cytotoxic"]
v0 <- venn.diagram(Cytotoxic_TMS_vs_aging_mice, NULL, 
                   category.names = c("aging mice cytotoxic", "TMS cytotoxic"), fill = c("red", "blue"))#, cex=c(1,1,1,)
grid.draw(v0)

overlaps <- calculate.overlap(Cytotoxic_TMS_vs_aging_mice)
overlaps <- rev(overlaps)
indx <- as.numeric(substr(names(overlaps),2,2))
for (i in 1:length(overlaps)){
  v0[[i+4]]$label <- paste(overlaps[[i]], collapse = "\n") 
}
grid.newpage()
grid.draw(v0)
write.csv(as.data.frame(overlaps$a3), "cytotoxic_genes_overlaps_TMS_aging_mice.csv")

venn.diagram(Cytotoxic_TMS_vs_aging_mice, "Venn_Diagram_TMS_aging_cytotoxic.jpeg",
             fill = c("red", "green"), cex = c(3,1,3), cat.pos = c(180,180), ext.percent = 0.001)

Cytotoxic_TMS_vs_aging_mice_Cytotoxic_TEM <- list()
Cytotoxic_TMS_vs_aging_mice_Cytotoxic_TEM$Aging_mice <- rownames(agingdata.markers_Cytotoxic_vs_TEM)
Cytotoxic_TMS_vs_aging_mice_Cytotoxic_TEM$TMS <- rownames(CD4_T_Cytotoxic_vs_TEM)
venn.diagram(Cytotoxic_TMS_vs_aging_mice_Cytotoxic_TEM, "Venn_Diagram_TMS_aging_cytotoxic_vs_TEM.jpeg",
             fill = c("red", "green"), cex = c(2,3,3), cat.pos = c(225,135), ext.percent = 0.001, cat.cex = 1)

Cytotoxic_TMS_vs_aging_mice_Exhausted_TEM <- list()
Cytotoxic_TMS_vs_aging_mice_Exhausted_TEM$Aging_mice <- rownames(agingdata.markers_Exhausted_vs_TEM)
Cytotoxic_TMS_vs_aging_mice_Exhausted_TEM$TMS <- rownames(CD4_T_Exhausted_vs_TEM)
venn.diagram(Cytotoxic_TMS_vs_aging_mice_Exhausted_TEM, "Venn_Diagram_TMS_aging_Exhausted_vs_TEM.jpeg",
             fill = c("red", "green"), cex = c(2,3,3), cat.pos = c(225,135), ext.percent = 0.001, cat.cex = 1)

Cytotoxic_TMS_vs_aging_mice_aTregs_rTregs <- list()
Cytotoxic_TMS_vs_aging_mice_aTregs_rTregs$Aging_mice <- rownames(agingdata.markers_aTregs_vs_rTregs)
Cytotoxic_TMS_vs_aging_mice_aTregs_rTregs$TMS <- rownames(CD4_T_aTregs_vs_rTregs)
venn.diagram(Cytotoxic_TMS_vs_aging_mice_aTregs_rTregs, "Venn_Diagram_TMS_aging_aTregs_rTregs.jpeg",
             fill = c("red", "green"), cex = c(2,3,3), cat.pos = c(225,135), ext.percent = 0.001, cat.cex = 1)

CD4_CTL_joind_markers <- list()
CD4_CTL_joind_markers$Humah_Aging_vs_Bladder_Cancer <- rownames(CD4_CTL.markers_Super_vs_tomor)
CD4_CTL_joind_markers$Humah_Aging_vs_Mouse_Aging <- rownames(CD4_CTL.markers_Super_vs_aging)
CD4_CTL_joind_markers$Human_vs_Mouse <- rownames(CD4_CTL.markers_mouse_vs_human)
CD4_CTL_joind_markers$Old_vs_Young <- rownames(CD4_CTL.markers_aging_vs_yuong)
CD4_CTL_joind_markers$Mouse_Aging_vs_Bladder_Cancer <- rownames(CD4_CTL.markers_aging_vs_tomor)
venn.diagram(CD4_CTL_joind_markers, "Venn_Diagram_CD4_CTL_joind_markers.jpeg", ext.percent = 0.001, cat.cex = 0.5, 
             fill = c("red", "yellow", "green", "blue", "purple"), cat.pos = c(0, 340, 200, 145, 0))
CD4_CTL_joind_markers_list <- calculate.overlap(CD4_CTL_joind_markers)
write.csv(CD4_CTL_joind_markers_list$a31, "CD4_CTL_joind_markers_list_a31.csv")
saveRDS(CD4_CTL_joind_markers_list, "CD4_CTL_joind_markers_list.rds")
