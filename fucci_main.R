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

datadir <- "/home/samba/storage0/shintaku/Fucci3.2 CAGE＿unpublished/STAR/"
wdir <- "/home/samba/storage0/shintaku/Fucci3.2 CAGE＿unpublished/"
source("/home/samba/storage0/shintaku/github/hunter/fucci_preprocess_data.R")
# download reference data from ensembl with biomaRt
source(file.path(rdir,'hunter_biomart_ref.R'))
# save count data with 10x format
source(file.path(rdir, 'hunter_preprocess_save_10x_format.R'))
# load data from 10x formatted files
source(file.path(rdir,"fucci_Seurat_load_dataset.R"))



