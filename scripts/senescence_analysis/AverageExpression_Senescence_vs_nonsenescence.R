Idents(Senis.big) <- Senis.big$Senescence.group

Senescence.average <- AverageExpression(Senis.big, return.seurat = TRUE, add.ident = "Age_group")
CellScatter(Senescence.average, cell1 = "Senescence.pos", cell2 = "Senescence.neg")
CellScatter(Senescence.average, cell1 = "Senescence.pos_Old", cell2 = "Senescence.neg_Old")
m <- Senescence.pos.markers_senescent_vs_nonsenescent %>% group_by(diffexpressed) %>% top_n(n = 30, wt = avg_logFC) %>% 
  top_n(n = 30, wt = p_val_adj)
DoHeatmap(Senescence.average, features = m$delabel, draw.lines = FALSE) + scale_fill_viridis()
Senescence.average <- AverageExpression(Senis.big, return.seurat = TRUE)
Idents(Senis.big) <- Senis.big$Age_group
levels(Senis.big) <- c("Young", "Old", "Supercentenarian")
Age_group.average <- AverageExpression(Senis.big, return.seurat = TRUE, add.ident = "Senescence.group")
DoHeatmap(Age_group.average, features = m$delabel, draw.lines = FALSE, 
          group.colors = c("#65EFFF", "#996035", "#331900")) + scale_fill_viridis()
Age_group.average <- AverageExpression(Senis.big, return.seurat = TRUE)
DoHeatmap(Age_group.average, features = m$delabel, draw.lines = FALSE, 
          group.colors = c("#65EFFF", "#996035", "#331900")) + scale_fill_viridis()

Idents(Senis.big) <- Senis.big$Age_group_2
levels(Senis.big) <- c("Young", "Old")
Age_group.average <- AverageExpression(Senis.big, return.seurat = TRUE, add.ident = "Senescence.group")
DoHeatmap(Age_group.average, features = m$delabel, draw.lines = FALSE, 
          group.colors = c("#65EFFF", "#996035")) + scale_fill_viridis()
Age_group.average <- AverageExpression(Senis.big, return.seurat = TRUE)
DoHeatmap(Age_group.average, features = m$delabel, draw.lines = FALSE, 
          group.colors = c("#65EFFF", "#996035")) + scale_fill_viridis()