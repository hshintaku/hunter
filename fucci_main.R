library(stringr)
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
source('./geom_split_violin.R')
rdir <- "/home/samba/storage0/shintaku/github/hunter"



# run following commands for the first
source("/home/samba/storage0/shintaku/github/hunter/fucci_preprocess_RNAseq_data.R")
# download reference data from ensembl with biomaRt

gene_list <- unique(data.frame(str_replace(allData$gene,"_intron","")))
colnames(gene_list) <- "gene"
source(file.path(rdir,'hunter_biomart_ref.R'))
# save count data with 10x format
source(file.path(rdir, 'hunter_preprocess_save_10x_format.R'))


# load data from 10x formatted files
# Sawano
# Sawano CAGE-RNA-seq
datadir <- "/home/samba/storage0/shintaku/Fucci3.2 CAGE＿unpublished/STAR/"
wdir <- "/home/samba/storage0/shintaku/Fucci3.2 CAGE＿unpublished/"
source(file.path(rdir,"fucci_cage_Seurat_load_dataset.R"))
# Kaneko
datadir <- "/home/samba/storage0/shintaku/20210324MiSeq016_Kaneko/"
wdir <- "/home/samba/storage0/shintaku/20210324MiSeq016_Kaneko/"
source(file.path(rdir,"./fucci_Kaneko_Seurat_load_dataset.R"))


source(file.path(rdir,'hunter_biomart_ref.R'))
gene_list <- data.frame(unique(c(rownames(fucci_cage),rownames(fucci))))
colnames(gene_list) <- "gene"
hs_ref <- unique(func.biomart.ref(hs_mart,"hgnc_symbol",gene_list,"hgnc_symbol"))
all_ref <-hs_ref


source('fucci_Seurat_cell_cycle.R')