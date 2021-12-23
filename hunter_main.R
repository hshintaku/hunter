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
datadir <- "/home/samba/sanger/shintaku/20211124HiSeqX006_TIG/"
#wdir <- "/home/samba/sanger/shintaku/20210728HiSeqX004_10x_cellranger/C01TIG-H12/outs/raw_feature_bc_matrix/"
wdir <- "/home/samba/public/shintaku/20211124HiSeqX006_TIG/"
rdir <- "/home/samba/public/shintaku/hunter/"

#
# cell_id_list2.txt contains all barcodes
# cell_id_list.txt contains selected barcodes by GC percent.
barcode <- read.table(file.path(rdir,"cell_id_list2.txt"))
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
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# annotate cells
source(file.path(rdir,"shiomi_Seurat_annotate_cells.R"))

# technical check, pca and umap clustering for cell typing
source(file.path(rdir,'shiomi_Seurat_10x_technical.R'))
# clustering
source(file.path(rdir,"shiomi_Seurat_clustering.R"))

# load FCS data
#indexdir =paste0(wdir,"index/")
indexdir <- "/home/samba/sanger/Shiomi/ELASTomicsindex"
channel <- c("Events","FSC","SSC","Venus","mCherry")
#c("Events","FSC","SSC","Venus","APC","mCherry")
source(file.path(rdir,'/io/hunter_Seurat_load_adt_data.R'))


#
# load cite-seq-count=FLD data as pbmc.tag
#
source(file.path(rdir,'io/hunter_Seurat_load_fld_data.R'))
# annotate cells with HTO and create FLD total
source(file.path(rdir,'20210816HiSeqX004_annotate_condition.R'))
#
pbmc[["FLD"]] <- CreateAssayObject(counts=new_pbmc.tag[c("FLD004","FLD010","FLD070","FLD150","FLD500","FLDtotal"),])
# check FLD by visualizing results
source(file.path(rdir,'20210816HiSeqX004_visualize_condition.R'))
# add HTO to Seurat object
pbmc[["HTO"]] <- CreateAssayObject(counts=pbmc.tag[c("T20CTL","T50CTL","TAZCTL"),])
#pbmc <- NormalizeData(pbmc, normalization.method = "CLR", scale.factor = 1e5,assay = "HTO")

DimPlot(pbmc,reduction="pca")
#p1<-DimPlot(pbmc,reduction="umap")
FeaturePlot(pbmc,features=c("FLD500","T20CTL","T50CTL","TAZCTL"),reduction = "pca")

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






