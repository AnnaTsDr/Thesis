agingdata <- readRDS("AgingData&tSNE&UMAP.rds")

ref_aging <- seurat_ref(agingdata, cluster_col = "Subset_Idan")
saveRDS(ref_aging, "ref_aging.rds")
ref_aging <- readRDS("ref_aging.rds")
ref_CD4_mouse_aging_Hezi <- seurat_ref(agingdata, cluster_col = "RNA_snn_res.0")
saveRDS(ref_CD4_mouse_aging_Hezi, "ref_CD4_mouse_aging_Hezi.rds")
