require(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(R.utils)
library(RCurl)
library(Matrix)
library(openxlsx)
library(dplyr)
library(Seurat)
library(SingleCellSignalR)
library(seqinr)
library(stringr)

# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/storage0/Shiomi/20210427MiSeq017Ana"
wdir <- "/home/samba/storage0/shintaku/20210427MiSeq017/"
rdir <- "/home/samba/storage0/shintaku/github/hunter"

barcode <- read.table(file.path("/home/samba/storage0/shintaku/github/hunter/cell_id_list.txt"))
barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))

source(file.path(rdir,"hunter_first_data_process.R"))

#
# you can restart from here
# load data from 10x formatted files
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
#
# create reference table with gene_short_name
# source(file.path(rdir,'hunter_biomart_ref.R'))
# gene_list <- data.frame(rownames(pbmc))
# colnames(gene_list) <- "gene"
# #hs_ref <- func.biomart.ref(hs_mart,gene_list,"hgnc_symbol")
# filter="mgi_symbol"
# symbol="mgi_symbol"
# ms_ref <- unique(func.biomart.ref(ms_mart,gene_list,filter,symbol))
# missing_ref <- subset(gene_list,!(gene %in% ms_ref$gene_short_name))
# adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
# colnames(adding_ref) <- colnames(ms_ref)
# rownames(adding_ref) <- adding_ref$ensembl_gene_id
# ms_ref <- rbind(adding_ref,ms_ref)


# technical check, pca and umap clustering for cell typing
source(file.path(rdir,'hunter_Seurat_technicalcheck.R'))

# load FCS data
#indexdir =paste0(wdir,"index/")
indexdir <- "/home/samba/storage0/Shiomi/hunterindex"
channel <- c("Events","FSC","SSC","Venus","mCherry")
#c("Events","FSC","SSC","Venus","APC","mCherry")

source(file.path(rdir,'hunter_Seurat_load_adt_data.R'))
#
# load cite-seq-count=FLD data
#
source(file.path(rdir,'hunter_Seurat_load_fld_data.R'))
#
# 
#
source(file.path(rdir,"shiomi_Seurat_cellcycle_dependence.R"))






# subset analysis: clustering subset and checking the mCherry expression
source(file.path(rdir,"hunter_Seurat_subset_analysis.R"))
source(file.path(rdir, "hunter_Seurat_subset_scatter.R"))

# WGCNA for gene network module
#source(file.path(rdir,"hunter_WGCNA_main.R"))

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




