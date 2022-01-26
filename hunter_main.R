require(readr)
library(plyr)
library(dplyr)
library(conflicted)
conflict_prefer("combine", "dplyr")
conflict_prefer("arrange", "dplyr")
library(tidyr)
library(ggplot2)
library(tidyverse)
library(R.utils)
library(RCurl)
library(Matrix)
library(openxlsx)
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
datadir <- "/home/samba/sanger/shintaku/20211026HiSeqX005_hunter/downsample/"
wdir <- "/home/samba/sanger/shintaku/20211026HiSeqX005_hunter/downsample/"
indexdir <- "/home/samba/sanger/shintaku/20211026HiSeqX005_hunter/index/"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
hepa2<- pbmc
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/sanger/shintaku/20210216HiSeqX002_HUNTER/downsample/"
wdir <- "/home/samba/sanger/shintaku/20210216HiSeqX002_HUNTER/downsample/"
indexdir <- "/home/samba/sanger/shintaku/20210216HiSeqX002_HUNTER/index/"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
hepa1 <- pbmc
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/sanger/shintaku/20211124HiSeqX006_hunter/downsample/"
wdir <- "/home/samba/sanger/shintaku/20211124HiSeqX006_hunter/downsample/"
indexdir <- "/home/samba/sanger/shintaku/20211124HiSeqX006_hunter/index/"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
hepa3 <- pbmc
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/sanger/shintaku/20220109HiSeqX008_hunter/"
wdir <- "/home/samba/sanger/shintaku/20220109HiSeqX008_hunter/"
indexdir <- "/home/samba/sanger/shintaku/20220109HiSeqX008_hunter/index"
channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
source(file.path(rdir,'io/hunter_Seurat_load_adt_data.R'))
hepa4 <- pbmc


datadir <- "/home/samba/public/shintaku/20220109HiSeqX008_hunter_10x/TK10x-02-S/"
wdir <- "/home/samba/public/shintaku/20220109HiSeqX008_hunter_10x/TK10x-02-S/"
#indexdir <- "/home/samba/sanger/shintaku/20211124HiSeqX006_hunter/index/"
#channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
hepa10x02S<-pbmc
datadir <- "/home/samba/public/shintaku/20220109HiSeqX008_hunter_10x/"
wdir <- "/home/samba/public/shintaku/20220109HiSeqX008_hunter_10x/TK10x-01-P/"
#indexdir <- "/home/samba/sanger/shintaku/20211124HiSeqX006_hunter/index/"
#channel <- c("Events","FSC","SSC","Venus","Azrite","mCherry")
source(file.path(rdir,"hunter_Seurat_load_dataset.R"))
hepa10x01P<-pbmc
#hepa.list <-list(hepa1,hepa2)
#anchors <- FindIntegrationAnchors(object.list = hepa.list)
#integrated <- IntegrateData(anchorset = anchors)
hepa10x01P[["batch"]]<-"01P"
hepa10x02S[["batch"]]<-"02S"
hepa10x01P[["percent.mt"]] <- PercentageFeatureSet(hepa10x01P, pattern = "^mt-")
hepa10x02S[["percent.mt"]] <- PercentageFeatureSet(hepa10x02S, pattern = "^mt-")
hepa10x <- merge(hepa10x01P,y=hepa10x02S)
FeatureScatter(hepa10x,feature1 = "nCount_RNA","nFeature_RNA",group.by = "batch")
umi_count <- hepa10x02S[["nFeature_RNA"]]
#umi_count <- umi_count[order(umi_count$nFeature_RNA),]
index <- order(umi_count$nFeature_RNA,decreasing = TRUE)
umi_count <- data.frame(umi_count[index,])
colnames(umi_count) <- "UMI"
umi_count$index <- as.numeric(rownames(umi_count))
umi_count$mtpercent <- as.numeric(unlist(hepa10x02S[["percent.mt"]]))

ggplot(umi_count,aes(x=index,y=UMI,fill=mtpercent))+geom_point()+scale_x_log10()+scale_y_log10()

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# annotate cells
source(file.path(rdir,"shiomi_Seurat_annotate_cells.R"))

source(file.path(rdir,'hunter_Seurat_technicalcheck.R'))
source(file.path(rdir,'hunter_Seurat_clustering.R'))
# load FCS data
#indexdir =paste0(wdir,"index/")
#c("Events","FSC","SSC","Venus","APC","mCherry")
# check cell cycle dependence 
source(file.path(rdir,"hunter_Seurat_cellcycle.R"))
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
