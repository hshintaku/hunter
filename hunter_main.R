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

# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20211124HiSeqX006_Islet/"
wdir <- "/home/samba/public/shintaku/20211124HiSeqX006_Islet/"
rdir <- "/home/samba/public/shintaku/github/hunter2"

barcode <- read.table(file.path(rdir,"cell_id_list.txt"))
barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))

#symbol="mgi_symbol"
#symbol="hgnc_symbol"
source(file.path(rdir,"hunter_first_data_process.R"))

#
#
# you can restart from here
# load data from 10x formatted files
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
huvec <- subset(pbmc, subset=species=="human")
islet <- subset(pbmc,subset=species=="rattus")

pbmc<-merge(huvec,y=islet)
pbmc<-islet
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000)
# technical check, pca and umap clustering for cell typing
source(file.path(rdir,'nashi_Seurat_clustering.R'))
# DESeq2
source(file.path(rdir,"bulk_DESeq2.R"))





# compute pseudotime and order cells along the gene expression
source(file.path(rdir,"shiomi_Seurat_monocle_pseudotime.R"))
#
# http://yulab-smu.top/clusterProfiler-book/index.html
#
# GO analysis
source(file.path(rdir,"hunter_clusterProfiler_GO.R"))
# pathway analysis
source(file.path(ridir,"hunter_clusterProfiler_GSEA.R"))






