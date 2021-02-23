require(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(R.utils) 
library(Matrix)
library(openxlsx)
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/watson/sanger/shintaku/20210216HiSeqX002/"
wdir <- "/home/watson/public/shintaku/HUNTER/"

barcode <- read.table("/home/watson/sanger/shintaku/HUNTER/RTbarcodes.txt")
# load functions for barcode decoding
source("/home/watson/sanger/shintaku/HUNTER/whitelist_encode.R")
# laod whitelist and check the batch effect
source('/home/watson/public/shintaku/HUNTER/hunter_preprocess_whitelist.R')

# preprocess the count data and load reference
source('/home/watson/public/shintaku/HUNTER/hunter_preprocess_data.R')

# save count data with 10x format
source('/home/watson/public/shintaku/HUNTER/hunter_preprocess_save_10x_format.R.R')

# load data from 10x formatted files
source("/home/watson/public/shintaku/HUNTER/hunter_Seurat_load_dataset.R")

# technical check, pca and umap clustering for cell typing
source('/home/watson/public/shintaku/HUNTER/hunter_Seurat_technicalcheck.R')

# load FCS data
indexdir ="/home/watson/sanger/shintaku/HUNTER/index/"
source('/home/watson/public/shintaku/HUNTER/hunter_Seurat_load_adt_data.R')

# subset analysis: clustering subset and checking the mCherry expression
source("/home/watson/public/shintaku/HUNTER/hunter_Seurat_subset_analysis.R")
source("/home/watson/public/shintaku/HUNTER/hunter_Seurat_subset_scatter.R")

# WGCNA for gene network module
source("/home/watson/public/shintaku/HUNTER/hunter_WGCNA_main.R")

#clusterProfiler for GO analysis
source("/home/watson/public/shintaku/HUNTER/hunter_clusterProfiler.R")

#SingleCellSignale.R
source("/home/watson/public/shintaku/HUNTER/hunter_SingleCellSignal.R")

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


# Ligand/Receptor analysis using SingleCellSignalR
#signal = cell_signaling(data=data,genes=all.genes,cluster=cluster)

# Visualization
#visualize(signal)
#intra = intra_network("S1PR1",data,all.genes,cluster,"cluster 1",signal = signal)


