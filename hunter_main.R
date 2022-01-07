require(readr)
library(plyr)
library(dplyr)
conflict_prefer("combine", "dplyr")
conflict_prefer("arrange", "dplyr")
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
datadir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/downsample/"
wdir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/downsample/"
indexdir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/index/"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
hepa2<- pbmc
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/downsample/"
wdir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/downsample/"
indexdir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/index/"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
hepa1 <- pbmc
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20211124HiSeqX006_hunter/downsample/"
wdir <- "/home/samba/public/shintaku/20211124HiSeqX006_hunter/downsample/"
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

active_E0771_vivo <- subset(regulonTargetsInfo,subset=TF==c("Rela","Klf2","Irf8"))
receptor_E0771_vivo <- signal$`GFP+-E0771vivo`

active_E0771_vivo_gene <- unique(active_E0771_vivo$gene)
receptor_E0771_vivo_gene <- unique(receptor_E0771_vivo$E0771vivo)
data=list(scenic=active_E0771_vivo_gene,signalR=receptor_E0771_vivo_gene)
venn.diagram(data, filename="senic_signalR.svg",
             imagetype="svg", height=5, width=5, fill=c(4,7), lty=2,
             scaled=F, cex=c(2,2,2), cat.pos=c(330,30), cat.dist=c(0.05,0.05), cat.cex=c(1.2,1.2))
active_E0771_vivo_gene[data$scenic %in% data$signalR]
