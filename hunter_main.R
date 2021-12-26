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
rdir <- "/home/samba/public/shintaku/github/hunter2/"
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/"
wdir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/"
indexdir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/index/"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
hepa2<- pbmc
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/"
wdir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/"
indexdir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/index/"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
hepa1 <- pbmc
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20211124HiSeqX006_hunter/"
wdir <- "/home/samba/public/shintaku/20211124HiSeqX006_hunter/"
indexdir <- "/home/samba/public/shintaku/20211124HiSeqX006_hunter/index/"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
hepa3 <- pbmc


hepa.list <-list(hepa1,hepa2)
anchors <- FindIntegrationAnchors(object.list = hepa.list)
integrated <- IntegrateData(anchorset = anchors)



# technical check, pca and umap clustering for cell typing
source(file.path(rdir,'hunter_Seurat_technicalcheck.R'))
source(file.path(rdir,'hunter_Seurat_clustering.R'))
# load FCS data
#indexdir =paste0(wdir,"index/")
#c("Events","FSC","SSC","Venus","APC","mCherry")
# check cell cycle dependence 
source(file.path(rdir,"shiomi_Seurat_cellcycle_dependence.R"))
#
# zonation
#
#source(file.path(rdir,"hunter_load_landmark_genes.R"))
# compute pseudotime and order cells along the gene expression
#source(file.path(rdir,"hunter_Seurat_pseudotime.R"))
# compute zonation via diffusion map 
source(file.path(rdir,"hunter_Seurat_diffusionmap.R"))
#
# http://yulab-smu.top/clusterProfiler-book/index.html
#
# GO analysis
source(file.path(rdir,"hunter_clusterProfiler_GO.R"))
# pathway analysis
source(file.path(ridir,"hunter_clusterProfiler_GSEA.R"))






