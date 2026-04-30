Tongue <- readRDS(file = "Clustred_Tongue.rds")
Heart_and_Aorta <- readRDS(file = "Clustred_Heart_and_Aorta.rds")
Marrow <- readRDS(file = "Clustred_Marrow.rds")
Mammary_Gland <- readRDS(file = "Clustred_Mammary_Gland.rds")
Fat <- readRDS(file = "Clustred_Fat.rds")
Kidney <- readRDS(file = "Clustred_Kidney.rds")
Liver <- readRDS(file = "Clustred_Liver.rds")
Lung <- readRDS(file = "Clustred_Lung.rds")
Limb_Muscle <- readRDS(file = "Clustred_Limb_Muscle.rds")
Pancreas <- readRDS(file = "Clustred_Pancreas.rds")
Spleen <- readRDS(file = "Clustred_Spleen.rds")
Thymus <- readRDS(file = "Clustred_Thymus.rds")
Bladder <- readRDS(file = "Clustred_Bladder.rds")
Skin <- readRDS(file = "Clustred_Skin.rds")
Large_Intestine <- readRDS(file = "Clustred_Large_Intestine.rds")
Trachea <- readRDS(file = "Clustred_Trachea.rds")

library(msigdbr)
library(SCINA)

CP_REACTOME_gene_sets = msigdbr(species = "mouse", category = "C2", subcategory = "CP:REACTOME")
msigdbr_list = split(x = CP_REACTOME_gene_sets$gene_symbol, f = CP_REACTOME_gene_sets$gs_name)
Senescence_genes_list <- list()
Senescence_genes_list$"REACTOME_CELLULAR_SENESCENCE" <- msigdbr_list[["REACTOME_CELLULAR_SENESCENCE"]]

Tongue_ <- SCINA(Tongue@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, convergence_rate = 0.9, 
                 sensitivity_cutoff = 1)
Heart_and_Aorta_ <- SCINA(Heart_and_Aorta@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, 
                          convergence_rate = 0.9, sensitivity_cutoff = 1)
Marrow_ <- SCINA(Marrow@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, convergence_rate = 0.9,
                 sensitivity_cutoff = 1)
Mammary_Gland_ <- SCINA(Mammary_Gland@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, 
                        convergence_rate = 0.9, sensitivity_cutoff = 1)
Fat_ <- SCINA(Fat@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, convergence_rate = 0.9, 
              sensitivity_cutoff = 1)
Kidney_ <- SCINA(Kidney@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, convergence_rate = 0.9,
                 sensitivity_cutoff = 1)
Liver_ <- SCINA(Liver@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, convergence_rate = 0.9,
                sensitivity_cutoff = 1)
Lung_ <- SCINA(Lung@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, convergence_rate = 0.9,
               sensitivity_cutoff = 1)
Limb_Muscle_ <- SCINA(Limb_Muscle@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, 
                      convergence_rate = 0.9, sensitivity_cutoff = 1)
Pancreas_ <- SCINA(Pancreas@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, 
                   convergence_rate = 0.9, sensitivity_cutoff = 1)
Spleen_ <- SCINA(Spleen@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, convergence_rate = 0.9,
                 sensitivity_cutoff = 1)
Thymus_ <- SCINA(Thymus@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, convergence_rate = 0.9,
                 sensitivity_cutoff = 1)
Bladder_ <- SCINA(Bladder@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, 
                  convergence_rate = 0.9, sensitivity_cutoff = 1)
Skin_ <- SCINA(Skin@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, convergence_rate = 0.9, 
               sensitivity_cutoff = 1)
Large_Intestine_ <- SCINA(Large_Intestine@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, 
                          convergence_rate = 0.9, sensitivity_cutoff = 1)
Trachea_ <- SCINA(Trachea@assays$RNA@data, Senescence_genes_list, max_iter = 100, rm_overlap = 0, 
                  convergence_rate = 0.9, sensitivity_cutoff = 1)


saveRDS(Tongue_, "Tongue_senescence.rds")
saveRDS(Heart_and_Aorta_, "Heart_and_Aorta_senescence.rds")
saveRDS(Marrow_, "Marrow_senescence.rds")
saveRDS(Mammary_Gland_, "Mammary_Gland_senescence.rds")
saveRDS(Fat_, "Fat_senescence.rds")
saveRDS(Kidney_, "Kidney_senescence.rds")
saveRDS(Liver_, "Liver_senescence.rds")
saveRDS(Lung_, "Lung_senescence.rds")
saveRDS(Limb_Muscle_, "Limb_Muscle_senescence.rds")
saveRDS(Pancreas_, "Pancreas_senescence.rds")
saveRDS(Spleen_, "Spleen_senescence.rds")
saveRDS(Thymus_, "Thymus_senescence.rds")
saveRDS(Bladder_, "Bladder_senescence.rds")
saveRDS(Skin_, "Skin_senescence.rds")
saveRDS(Large_Intestine_, "Large_Intestine_senescence.rds")
saveRDS(Trachea_, "Trachea_senescence.rds")

Tongue$Senescence_reactome <- Tongue_$cell_labels
Heart_and_Aorta$Senescence_reactome <- Heart_and_Aorta_$cell_labels
Marrow$Senescence_reactome <- Marrow_$cell_labels
Mammary_Gland$Senescence_reactome <- Mammary_Gland_$cell_labels
Fat$Senescence_reactome <- Fat_$cell_labels
Kidney$Senescence_reactome <- Kidney_$cell_labels
Liver$Senescence_reactome <- Liver_$cell_labels
Lung$Senescence_reactome <- Lung_$cell_labels
Limb_Muscle$Senescence_reactome <- Limb_Muscle_$cell_labels
Pancreas$Senescence_reactome <- Pancreas_$cell_labels
Spleen$Senescence_reactome <- Spleen_$cell_labels
Thymus$Senescence_reactome <- Thymus_$cell_labels
Bladder$Senescence_reactome <- Bladder_$cell_labels
Skin$Senescence_reactome <- Skin_$cell_labels
Large_Intestine$Senescence_reactome <- Large_Intestine_$cell_labels
Trachea$Senescence_reactome <- Trachea_$cell_labels

saveRDS(Tongue, file = "Clustred_Senescence_Annotated_Tongue.rds")
saveRDS(Heart_and_Aorta, file = "Clustred_Senescence_Annotated_Heart_and_Aorta.rds")
saveRDS(Marrow, file = "Clustred_Senescence_Annotated_Marrow.rds")
saveRDS(Mammary_Gland, file = "Clustred_Senescence_Annotated_Mammary_Gland.rds")
saveRDS(Fat, file = "Clustred_Senescence_Annotated_Fat.rds")
saveRDS(Kidney, file = "Clustred_Senescence_Annotated_Kidney.rds")
saveRDS(Liver, file = "Clustred_Senescence_Annotated_Liver.rds")
saveRDS(Lung, file = "Clustred_Senescence_Annotated_Lung.rds")
saveRDS(Limb_Muscle, file = "Clustred_Senescence_Annotated_Limb_Muscle.rds")
saveRDS(Pancreas, file = "Clustred_Senescence_Annotated_Pancreas.rds")
saveRDS(Spleen, file = "Clustred_Senescence_Annotated_Spleen.rds")
saveRDS(Thymus, file = "Clustred_Senescence_Annotated_Thymus.rds")
saveRDS(Bladder, file = "Clustred_Senescence_Annotated_Bladder.rds")
saveRDS(Skin, file = "Clustred_Senescence_Annotated_Skin.rds")
saveRDS(Large_Intestine, file = "Clustred_Senescence_Annotated_Large_Intestine.rds")
saveRDS(Trachea, file = "Clustred_Senescence_Annotated_Trachea.rds")
