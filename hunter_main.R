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
#library(VennDiagram)

# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/sanger/shintaku/20210728HiSeqX004_10x_TIG/"
#wdir <- "/home/samba/sanger/shintaku/20210323MiSeq015Ana10X/"
wdir <- "/home/samba/public/shintaku/20210728HiSeqX004_10x_TIG/"
rdir <- "/home/samba/public/shintaku/hunter/"

#
# cell_id_list2.txt contains all barcodes
# cell_id_list.txt contains selected barcodes by GC percent.
#barcode <- read.table(file.path(rdir,"cell_id_list2.txt"))
barcode <- read.table(file.path("/home/samba/public/shintaku/cellranger-6.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"))
rownames(barcode)<-barcode$V1
barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))

# 10x for big data, facs for small data
source(file.path(rdir,"10x_first_data_process.R"))
# preprocess FLD data
source(file.path(rdir,"preprocess/preprocess_FLD_data.R"))

#
#
# you can restart from here
# load data from 10x formatted files
source(file.path(rdir,"/io/hunter_Seurat_load_dataset.R"))
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 500)
# annotate cells
source(file.path(rdir,"shiomi_Seurat_annotate_cells.R"))

# technical check, pca and umap clustering for cell typing
source(file.path(rdir,'shiomi_Seurat_10x_technical.R'))

# load FCS data
#indexdir =paste0(wdir,"index/")
indexdir <- "/home/samba/sanger/Shiomi/ELASTomicsindex"
channel <- c("Events","FSC","SSC","Venus","mCherry")
#c("Events","FSC","SSC","Venus","APC","mCherry")
source(file.path(rdir,'/io/hunter_Seurat_load_adt_data.R'))


#
# load cite-seq-count=FLD data
#
source(file.path(rdir,'io/hunter_Seurat_load_fld_data.R'))
#
source(file.path(rdir,"shiomi_fld_external_control_analysis.R"))
# first overview
# analyze data with PCA and UMAP
# find clusters and marker genes
#source(file.path(rdir,"shiomi_Seurat_technicalcheck.R"))

# check cell cycle dependence 
source(file.path(rdir,"shiomi_Seurat_cellcycle_dependence.R"))
# find marker genes
source(file.path(rdir,"shiomi_Seurat_Marker.genes.R"))
# compute pseudotime and order cells along the gene expression
source(file.path(rdir,"shiomi_Seurat_monocle_pseudotime.R"))

#
# http://yulab-smu.top/clusterProfiler-book/index.html
#
# GO analysis
source(file.path(rdir,"hunter_clusterProfiler_GO.R"))
# pathway analysis
source(file.path(ridir,"hunter_clusterProfiler_GSEA.R"))






