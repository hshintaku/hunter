library(Seurat)
library(stringr)
library(dplyr)

#source(file.path(rdir,"io/preprocess_FLD_data.R"))
FLDmapALL <- load.fld(datadir,"read_count",barcode)
colnames(FLDmapALL) <- c("romin","FLD004","FLD010","FLD070","FLD500","FLDcon","Unmapped","GC")
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



FLDcomb <- data.frame(t(rbind(pbmc[["FLD"]]@data['FLD004',],
                              pbmc[["FLD"]]@data['FLD010',],
                              pbmc[["FLD"]]@data['FLD070',],
                              pbmc[["FLD"]]@data['FLD500',],
                              pbmc[["FLD"]]@data['GC',],
                              pbmc[["ADT"]]@data["Venus",])))
colnames(FLDcomb) <- c("FLD004","FLD010","FLD070","FLD500","GC","Venus")
FLDcomb <- subset(FLDcomb,Venus!=0) # remoce no cell
#FLDcomb$GC <- FLDsub[,]$GC
#ggplot(FLDcomb,aes(x=Venus,y=FLD070,color=GC))+geom_point()+scale_x_log10(limits=c(1,1e5))+scale_y_log10(limits=c(1,1e5))

FeatureScatter(pbmc,feature1="Venus",feature2="FLD004")+scale_x_log10()+scale_y_log10(limits=c(1,1e6))

FLD_high_GC <- subset(FLDcomb,GC>0.3)
FLDmelt <- melt(FLD_high_GC,id.vars = c("Venus","GC"))
ggplot(FLDmelt,aes(x=Venus+1,y=value+1,color=variable))+geom_point()+scale_y_log10()+scale_x_log10()
cor(FLD_high_GC, method="spearman")

rm(seladt.FLD.csv,pbmc.adt.FLD,FLDcomb,FLDmelt)

