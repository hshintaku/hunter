require(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(R.utils) 
library(Matrix)
library(openxlsx)
library(dplyr)
library(Seurat)
library(SingleCellSignalR)
library(seqinr)

# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/storage0/shintaku/20210216HiSeqX002/"
wdir <- "/home/samba/storage0/shintaku/HUNTER/"
rdir <- "/home/samba/storage0/shintaku/github/hunter"

barcode <- read.table(file.path("/home/samba/storage0/shintaku/HUNTER/RTbarcodes.txt"))
# load functions for barcode decoding
source(file.path(rdir,"whitelist_encode.R"))
# laod whitelist and check the batch effect
source(file.path(rdir,'hunter_preprocess_whitelist.R'))

# preprocess the count data and load reference
source(file.path(rdir,'hunter_preprocess_data.R'))

# download reference data from ensembl with biomaRt
gene_list <- unique(data.frame(str_replace(allData$gene,"_intron","")))
colnames(gene_list) <- "gene"
source(file.path(rdir,'hunter_biomart_ref.R'))


# save count data with 10x format
source(file.path(rdir, 'hunter_preprocess_save_10x_format.R'))


#
# you can restart from here
# load data from 10x formatted files
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
gene_list <- data.frame(rownames(pbmc))
colnames(gene_list) <- "gene"
source(file.path(rdir,'hunter_biomart_ref.R'))

# technical check, pca and umap clustering for cell typing
source(file.path(rdir,'hunter_Seurat_technicalcheck.R'))

# load FCS data
indexdir =paste0(wdir,"index/")
source(file.path(rdir,'hunter_Seurat_load_adt_data.R'))

# subset analysis: clustering subset and checking the mCherry expression
source(file.path(rdir,"hunter_Seurat_subset_analysis.R"))
source(file.path(rdir, "hunter_Seurat_subset_scatter.R"))

# WGCNA for gene network module
source(file.path(rdir,"hunter_WGCNA_main.R"))

#clusterProfiler for GO analysis
source(file.path(rdir,"hunter_clusterProfiler.R"))

#SingleCellSignale.R
source(file.path(rdir,"hunter_SingleCellSignal.R"))

library(VennDiagram)

corr_module_normGFP <- names(datExpr)[moduleColors=="purple"]
acorr_module_normGFP <- names(datExpr)[moduleColors=="black"]
normGFP_module <- c(corr_module_normGFP,acorr_module_normGFP)

corr_moudle_mCherry <- names(datExpr)[moduleColors=="yellow"]
acorr_module_mCherry <- names(datExpr)[moduleColors=="green"]
mCherry_module <- c(corr_moudle_mCherry,acorr_module_mCherry)

PC_1_gene <- PCASigGenes(object=AML,pcs.use=1,pval.cut=0.1)
PC_2_gene <- PCASigGenes(object=AML,pcs.use=2,pval.cut=0.1)

gene_list<- list(normGFP=corr_module_normGFP, mCherry=mCherry_module,PC_1=PC_1_gene, PC_2=PC_2_gene)
venn.diagram(gene_list,filename = file.path(wdir,"gene.jpg"), fill=c(2,3,4,5), alpha=0.4, lty=3)


#VlnPlot(pbmc, features = c("pg-GAPDH", "pg-S100A2"), slot = "counts", log = TRUE)
#FeaturePlot(pbmc, features = c("hs-MT-ND4", "pg-GAPDH"))




