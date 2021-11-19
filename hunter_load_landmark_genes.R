library(R.matlab)
library(tools)
#https://www.nature.com/articles/nbt.4231#MOESM96
#hepa <- subset(pbmc,subset=plate=="p02",invert=TRUE)
#hepa <- subset(hepa,subset=plate=="P15",invert=TRUE)

pbmc_monocle <- as.CellDataSet(pbmc)
pbmc_monocle <- estimateSizeFactors(pbmc_monocle)
pbmc_monocle <- estimateDispersions(pbmc_monocle)
disp_table <- dispersionTable(pbmc_monocle)
rownames(disp_table)<-disp_table$gene_id

genes_zonation <- readMat("/home/samba/public/shintaku/matlab_code/Zonation_params.mat")
genes_cv <- data.frame(unlist(genes_zonation$genes.cv))
colnames(genes_cv)<-"gene_id"
rownames(genes_cv)<- genes_cv$gene_id
genes_pn <- data.frame(unlist(genes_zonation$genes.pn))
colnames(genes_pn)<-"gene_id"
rownames(genes_pn)<-genes_pn$gene_id

ordering_genes <- rbind(genes_pn,genes_cv)