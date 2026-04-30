library(biomaRt)

agingdata <- readRDS("AgingDataClustered.rds")

genes <- rownames(agingdata)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "jul2016.archive.ensembl.org")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "jul2016.archive.ensembl.org")
genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes , mart = mouse, 
                 attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

humanx <- unique(genesV2[,1])
genesV2_2 <- genesV2[which(humanx %in% genesV2[,1]),]
agingdata_2 <- subset(agingdata, features = humanx)
genesV2_2 <- genesV2_2[,2]
rownames(agingdata_2@assays$RNA@counts) <- genesV2_2
rownames(agingdata_2@assays$RNA@data) <- genesV2_2

agingdata_2 <- subset(agingdata_2, ident = "Cytotoxic")

saveRDS(agingdata_2, "AgingDataCytotoxic.rds")

mouse_to_human <- genesV2
mouse_to_human <- mouse_to_human[!duplicated(mouse_to_human$MGI.symbol),]

g <- rownames(agingdata)
common_genes <- intersect(g, mouse_to_human[,1])
g <- match(common_genes, mouse_to_human)
agingdata <- agingdata[common_genes,]
mouse_to_human <- arrange(mouse_to_human, common_genes)
row.names(mouse_to_human) <- mouse_to_human$MGI.symbol
mouse_to_human <- mouse_to_human[match(common_genes,mouse_to_human$MGI.symbol),]
mouse_to_human <- merge(genesV2, humanx)
row.names(agingdata@assays$RNA) <- row.names(mouse_to_human)

saveRDS(mouse_to_human, "mouse_vs_human_genes_for_analysis.rds")
saveRDS(genesV2, "mouse_vs_human_genes.rds")
write.csv(genesV2, "mouse_vs_human_genes.csv")
genesV2 <- readRDS("mouse_vs_human_genes.rds")

#genesV2_sup <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id', 'gene_biotype'), filters = "ensembl_gene_id", 
#                     values = genes_sup , mart = human, uniqueRows=T)

write.csv(genesV2_sup, "biomart_genes_id_to_name.csv")

convertMouseGeneList <- function(x){
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "mar2016.archive.ensembl.org")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "mar2016.archive.ensembl.org")
  genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, 
                   attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

ensembl84 <- useEnsembl(biomart = 'genes', 
                        dataset = 'hsapiens_gene_ensembl',
                        version = 84)


Supercentenarians.lognorm <- read.table(file = "Centenarians/02.UMI.lognorm.txt.gz")

genesV4 <- read.delim("HOM_MouseHumanSequence.txt")
genesV4 <- genesV4[,1:4]
g <- data.table()
g$ID <- "0"
g$mouse_symbol <- "0"
g$human_symbol <- "0"
for (n in 1:length(genesV4[,1])) {
  for (m in n:length(genesV4[,1])) {
    if (genesV4[n,1] == genesV4[m,1]){ 
    g[n,1] <- genesV4[n,1]
    g[n,2] <- genesV4[n,4]
    g[n,3] <- genesV4[m,4]
    }
  }
}
