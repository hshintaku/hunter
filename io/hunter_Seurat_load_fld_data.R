library(Seurat)
library(stringr)
library(dplyr)

#source(file.path(rdir,"shiomi_preprocess_FLD_data.R"))
FLDmapALL <- load.fld(datadir,"umi_count",barcode)
colnames(FLDmapALL) <- c("romin","FLD004","FLD010","FLD070","FLD500","FLDcon","Unmapped","GC")
rownames(FLDmapALL) <- toupper(rownames(FLDmapALL))

# adding FLD data
seladt.FLD.csv <- t(FLDmapALL[cellids,]) # extract cells detected in RNA-seq
pbmc.adt.FLD <-as.sparse(seladt.FLD.csv) # convert the format to sparse matrix
pbmc[["FLD"]] <- CreateAssayObject(counts=pbmc.adt.FLD)

#pbmc <- NormalizeData(pbmc,assay="ADT",normalization.method = "CLR",margin=2)
#pbmc <- NormalizeData(pbmc,assay="FLD",normalization.method = "CLR",margin=2)

FLDcomb <- data.frame(t(rbind(pbmc[["FLD"]]@data['FLD004',],
                              pbmc[["FLD"]]@data['FLD010',],
                              pbmc[["FLD"]]@data['FLD070',],
                              pbmc[["FLD"]]@data['FLD500',],
                              pbmc[["ADT"]]@data["Venus",])))
colnames(FLDcomb) <- c("FLD004","FLD010","FLD070","FLD500","Venus")
#FLDcomb$GC <- FLDsub[,]$GC
#ggplot(FLDcomb,aes(x=Venus,y=FLD070,color=GC))+geom_point()+scale_x_log10(limits=c(1,1e5))+scale_y_log10(limits=c(1,1e5))

cor(FLDcomb, method="spearman")

rm(seladt.FLD.csv,pbmc.adt.FLD,FLDcomb)

