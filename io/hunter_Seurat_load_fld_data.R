library(Seurat)
library(stringr)
library(dplyr)

#source(file.path(rdir,"shiomi_preprocess_FLD_data.R"))
FLDmapALL <- load.fld(datadir,"read_count",barcode)
colnames(FLDmapALL) <- c("romin","FLD004","FLD010","FLD070","FLD500","FLDcon","Unmapped","GC")
rownames(FLDmapALL) <- toupper(rownames(FLDmapALL))

# adding FLD data
seladt.FLD.csv <- FLDmapALL[cellids,] # extract cells detected in RNA-seq
rownames(seladt.FLD.csv) <- cellids
seladt.FLD.csv[is.na(seladt.FLD.csv$romin),]<-0
pbmc.adt.FLD <-as.sparse(t(seladt.FLD.csv)) # convert the format to sparse matrix
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

FeatureScatter(pbmc,feature1="Venus",feature2="FLD500")+scale_x_log10()+scale_y_log10()

FLDmelt <- melt(FLDcomb,id.vars = "Venus")
ggplot(FLDmelt,aes(x=Venus+1,y=value+1,color=variable))+geom_point()+scale_x_log10()+scale_y_log10()
cor(FLDcomb, method="spearman")

rm(seladt.FLD.csv,pbmc.adt.FLD,FLDcomb,FLDmelt)

