# list of genes are obtained from https://www.nature.com/articles/nbt.4231#Sec32
# Paired-cell sequencing enables spatial gene expression mapping of liver endothelial cells Nat Biotech 2018
#
landmark_genes_cv <- read.csv('/home/samba/pihome/2021/Shintaku/genes_cv.csv',header = F)
landmark_genes_pn <- read.csv('/home/samba/pihome/2021/Shintaku/genes_pn.csv',header = F)
#
#

vivo <- subset(x=pbmc, subset=dish =="E00")

landmark_genes_cv <- subset(landmark_genes_cv, (V1 %in% ms_ref$gene_short_name))
landmark_genes_cv$type <- "pericentral"
landmark_genes_pn <- subset(landmark_genes_pn, (V1 %in% ms_ref$gene_short_name))
landmark_genes_pn$type <- "periportal"
landmark_genes_combine <- rbind(landmark_genes_cv,landmark_genes_pn)

landmark_genes_exp <- vivo[["RNA"]]@data[landmark_genes_combine$V1,]
rownames(landmark_genes_combine) <- landmark_genes_combine$V1
landmark_genes_type <- data.frame(landmark_genes_combine$type)
rownames(landmark_genes_type) <- landmark_genes_combine$V1

library(pheatmap)

pheatmap(landmark_genes_exp,
         annotation_row = landmark_genes_type,
         cluster_cols = TRUE,cluster_rows = FALSE,scale="row")
