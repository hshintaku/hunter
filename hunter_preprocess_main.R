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
datadir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/downsample/"
wdir <- "/home/samba/public/shintaku/20210216HiSeqX002_HUNTER/downsample/"

# decode the single cell data from whitelist of UMI-tools output
datadir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/downsample/"
wdir <- "/home/samba/public/shintaku/20211026HiSeqX005_hunter/downsample/"

datadir <- "/home/samba/public/shintaku/20211124HiSeqX006_hunter/downsample/"
wdir <- "/home/samba/public/shintaku/20211124HiSeqX006_hunter/downsample/"

datadir <- "/home/samba/public/shintaku/20220109HiSeqX008_hunter/"
wdir <- "/home/samba/public/shintaku/20220109HiSeqX008_hunter/"

rdir <- "/home/samba/public/shintaku/github/hunter2/"

barcode <- read.table(file.path("/home/samba/public/Program/cellranger-6.1.2/lib/python/cellranger/barcodes/3M-february-2018.txt"))
barcode <- read.table(file.path(rdir,"cell_id_list.txt"))

barcode$GC <- as.numeric(lapply(lapply(as.character(barcode$V1),s2c),GC))

#
symbol="mgi_symbol"
symbol="hgnc_symbol"
symbol="rgd_symbol"

filter="ensembl_gene_id"
#filter="mgi_symbol"

# load functions for barcode decoding
source(file.path(rdir,"util/whitelist_encode.R"))
# laod whitelist and check the batch effect
source(file.path(rdir,'preprocess/preprocess_whitelist.R'))
active_barcode <- barcode[sort(unique(allencoded$first_index)),]
encode_barcode=TRUE
# preprocess the count data and load reference
source(file.path(rdir,'preprocess/preprocess_RNAseq_data.R'))
#
# download reference data from ensembl with biomaRt
gene_list <- unique(data.frame(str_replace(allData$gene,"_intron","")))
colnames(gene_list) <- "gene"
#
source(file.path(rdir,'util/hunter_biomart_ref.R'))
#hs_ref <- func.biomart.ref(hs_mart,gene_list,"hgnc_symbol")

if (symbol=="mgi_symbol"){
  ms_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(ms_mart,gene_list,filter,symbol))
}else if(symbol=="hgnc_symbol"){
  hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(hs_mart,gene_list,filter,symbol))
}else{
  rgd_mart <- useMart(biomart="ensembl", dataset="rnorvegicus_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(rgd_mart,gene_list,filter,symbol))
}

#
#pig_mart <- useMart(biomart="ensembl", dataset="sscrofa_gene_ensembl")


missing_ref <- subset(gene_list,!(gene %in% ms_ref$ensembl_gene_id))
adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
colnames(adding_ref) <- colnames(ms_ref)
rownames(adding_ref) <- adding_ref$ensembl_gene_id
ms_ref <- rbind(adding_ref,ms_ref)
rm(missing_ref,adding_ref)

# save count data with 10x format
source(file.path(rdir, 'preprocess/preprocess_save_10x_format.R'))
