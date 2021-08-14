library(Seurat)
library(stringr)
library(dplyr)

source(file.path(rdir,"preprocess/preprocess_FLD_data.R"))
FLDmapALL <- load.fld(datadir,"umi_count",barcode,FALSE)
tagnames <- rownames(FLDmapALL)
tagnames <- data.frame(strsplit(tagnames,"-"))
rownames(FLDmapALL) <- tagnames[1,]
rownames(FLDmapALL) <- toupper(rownames(FLDmapALL))
tFLDmapALL <- t(FLDmapALL)

# extract cellids shared with cDNA
cellids_fld <- colnames(FLDmapALL)
selcellids <- intersect(cellids,cellids_fld)
seladt.FLD.csv <- FLDmapALL[,selcellids] # extract cells detected in RNA-seq

#create empty data frame
pbmc.fld<-data.frame((matrix(0,nrow=nrow(seladt.FLD.csv),ncol=length(cellids))))
rownames(pbmc.fld)<-rownames(seladt.FLD.csv)
colnames(pbmc.fld)<-cellids

#replace the created data frame with FLD data
pbmc.fld[,selcellids]<-seladt.FLD.csv[,selcellids]

#seladt.FLD.csv[is.na(seladt.FLD.csv$romin),]<-0

#pbmc.adt.FLD <-as.sparse(t(seladt.FLD.csv)) # convert the format to sparse matrix

pbmc[["FLD"]] <- CreateAssayObject(counts=pbmc.fld)

p0<-DimPlot(pbmc,reduction="pca")
#p1<-DimPlot(pbmc,reduction="umap")
p2<-FeaturePlot(pbmc,features=c("FLD004","T20CTL","T50CTL","TAZCTL"),reduction = "pca")
p0+p2

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

