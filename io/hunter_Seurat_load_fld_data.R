library(Seurat)
library(stringr)
library(dplyr)

source(file.path(rdir,"preprocess/preprocess_FLD_data.R"))
FLDmapALL <- load.fld(datadir,"umi_count",barcode,FALSE)
tagnames <- rownames(FLDmapALL)
tagnames <- data.frame(strsplit(tagnames,"-"))
rownames(FLDmapALL) <- tagnames[1,]
rownames(FLDmapALL) <- toupper(rownames(FLDmapALL))

# extract cellids shared with cDNA
cellids_fld <- colnames(FLDmapALL)
selcellids <- intersect(cellids,cellids_fld)
seladt.FLD.csv <- FLDmapALL[,selcellids] # extract cells detected in RNA-seq
rm(FLDmapALL)
#create empty data frame with matchig rows and cols of RNA-seq
pbmc.tag<-data.frame((matrix(0,nrow=nrow(seladt.FLD.csv),ncol=length(cellids))))
rownames(pbmc.tag)<-rownames(seladt.FLD.csv)
colnames(pbmc.tag)<-cellids
#replace the created data frame with FLD data
pbmc.tag[,selcellids]<-seladt.FLD.csv[,selcellids]
rm(cellids_fld,selcellids)
#seladt.FLD.csv[is.na(seladt.FLD.csv$romin),]<-0

#pbmc.adt.FLD <-as.sparse(t(seladt.FLD.csv)) # convert the format to sparse matrix
source(file.path(rdir,'20210816HiSeqX004_annotate_condition.R'))
pbmc[["FLD"]] <- CreateAssayObject(counts=new_pbmc.tag[c("FLD004","FLD010","FLD070","FLD150","FLD500","FLDtotal"),])
source(file.path(rdir,'20210816HiSeqX004_visualize_condition.R'))


#pbmc <- NormalizeData(pbmc, normalization.method = "CLR", margin=2,assay = "FLD")
pbmc[["HTO"]] <- CreateAssayObject(counts=pbmc.tag[c("T20CTL","T50CTL","TAZCTL"),])
pbmc <- NormalizeData(pbmc, normalization.method = "CLR", scale.factor = 1e5,assay = "HTO")

DimPlot(pbmc,reduction="pca")
#p1<-DimPlot(pbmc,reduction="umap")
FeaturePlot(pbmc,features=c("FLD500","T20CTL","T50CTL","TAZCTL"),reduction = "pca")


rm(seladt.FLD.csv,pbmc.adt.FLD,FLDcomb,FLDmelt)

