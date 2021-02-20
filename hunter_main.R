require(readr)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(R.utils) 
library(Matrix)
library(openxlsx)
# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/watson/sanger/shintaku/20210216HiSeqX002/"
wdir <- "/home/watson/public/shintaku/HUNTER/"

barcode <- read.table("/home/watson/sanger/shintaku/HUNTER/RTbarcodes.txt")
# load whitelist and check the batch effect
source("/home/watson/sanger/shintaku/HUNTER/whitelist_encode.R")
source('/home/watson/public/shintaku/HUNTER/hunter_whitelist.R')

# preprocess the count data and load reference
source('/home/watson/public/shintaku/HUNTER/hunter_preprocess_data.R')

# save count data with 10x format
source('/home/watson/public/shintaku/HUNTER/hunter_10x_data.R')

# load data from 10x formatted files
source("/home/watson/public/shintaku/HUNTER/hunter_Seurat_load_dataset.R")

# technical check, pca and umap clustering for cell typing
source('/home/watson/public/shintaku/HUNTER/hunter_Seurat_technicalcheck.R')

# load FCS data
indexdir ="/home/watson/sanger/shintaku/HUNTER/index/"
source('/home/watson/public/shintaku/HUNTER/hunter_Seurat_adt_data_io.R')

# subset analysis: clustering subset and checking the mCherry expression
source("/home/watson/public/shintaku/HUNTER/hunter_Seurat_subset_analysis.R")
source("/home/watson/public/shintaku/HUNTER/hunter_Seurat_subset_scatter.R")

