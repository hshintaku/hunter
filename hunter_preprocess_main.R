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
datadir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/"
wdir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/"

# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/"
wdir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/"

rdir <- "/home/samba/public/shintaku/github/hunter2/"

barcode <- read.table(file.path(rdir,"cell_id_list.txt"))
barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))

#
symbol="mgi_symbol"
#symbol="hgnc_symbol"
source(file.path(rdir,"hunter_first_data_process.R"))
