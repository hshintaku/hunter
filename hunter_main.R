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
library(VennDiagram)
library(stringdist)
library(biomaRt)

# load data from 10x formatted files

# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/"
wdir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/"
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
hepa2<- pbmc
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/"
wdir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/"
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
hepa1 <- pbmc
allcell<- merge(hepa1, y = hepa2,  project = "hunter")

allcell <- subset(allcell,subset = gate ==c("g1"),invert=TRUE)
allcell <- subset(allcell,subset = gate ==c("g4"),invert=TRUE)
allcell <- subset(allcell,subset = plate ==c("p04"),invert=TRUE)
allcell <- subset(allcell,subset = plate ==c("p05"),invert=TRUE)


hepa <- subset(allcell,subset = plate ==c("p02"),invert=TRUE)
hepa <- subset(hepa,subset = plate ==c("P15"),invert=TRUE)



pbmc<-hepa
pbmc <- liver#allcell

cellids <- colnames(pbmc)



pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 500)

# technical check, pca and umap clustering for cell typing
source(file.path(rdir,'hunter_Seurat_technicalcheck.R'))
source(file.path(rdir,'hunter_Seurat_clustering.R'))
# load FCS data
#indexdir =paste0(wdir,"index/")
indexdir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/index/"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
#c("Events","FSC","SSC","Venus","APC","mCherry")
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
# check cell cycle dependence 
source(file.path(rdir,"shiomi_Seurat_cellcycle_dependence.R"))
# compute pseudotime and order cells along the gene expression
source(file.path(rdir,"shiomi_Seurat_monocle_pseudotime.R"))

#
# http://yulab-smu.top/clusterProfiler-book/index.html
#
# GO analysis
source(file.path(rdir,"hunter_clusterProfiler_GO.R"))
# pathway analysis
source(file.path(ridir,"hunter_clusterProfiler_GSEA.R"))






