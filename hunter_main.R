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

# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/sanger/shintaku/20210708HiSeqX003Piezo/"
#wdir <- "/home/samba/sanger/shintaku/20210323MiSeq015Ana10X/"
wdir <- "/home/samba/public/shintaku/Piezo/"
rdir <- "/home/samba/public/shintaku/hunter"

barcode <- read.table(file.path(rdir,"cell_id_list.txt"))
barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))

source(file.path(rdir,"hunter_first_data_process.R"))
#
#
piezo.data <- Read10X(data.dir = wdir)
piezo <- CreateSeuratObject(counts = piezo.data, project = "pbmc3k", min.cells = 1, min.features = 100)
piezo_count <- data.frame(piezo[["RNA"]]@data[ms_ref$gene_short_name,])

