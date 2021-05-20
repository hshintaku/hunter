library(Seurat)
library(stringr)
library(dplyr)
library(seqinr)
library(ggplot2)

wFdir <- "/home/samba/storage0/Shiomi/20210427MiSeq017Ana/"
sample_text <- "/home/samba/storage0/Shiomi/20210427MiSeq017Ana/citesample_Mix.txt"
sampleID <- read.table(sample_text) 
FLDGC <- matrix(nrow=10, ncol=0)
for (i in 1:19){
  FLDmap <- matrix(nrow=6, ncol=0)
  samfol = paste(file.path(wFdir), "CITE-seq/", sampleID[i,1], "/umi_count/",sep="")
  FLD.data <- Read10X(data.dir = samfol, gene.column=1)
  iName = substr(sampleID[i,1], 1, 10)
  romin <- whitelist.umi_tools.encode(colnames(FLD.data),barcode$V1)
  FLData <- as.data.frame(t(as.matrix(FLD.data)))
  nameF <- colnames(FLData)
  FLData <- FLData %>% 
    mutate(count = select(., one_of(nameF)) %>% 
             rowSums(na.rm = T))
  FLData <- as.data.frame(cbind(FLData, romin))
  FLData$GC <- as.numeric(lapply(lapply(as.character(rownames(FLData)),s2c),GC))
  rownames(FLData) <- paste (iName, rownames(FLData), sep = "_")
  FLData$Lot <- iName
  FLDGC <-rbind(FLDGC, FLData)
}

FLDGC$GC <- paste (FLDGC$GC, "/", sep = "")
FLDGC1 <- filter(FLDGC, str_detect(Lot, "rg"))
ggplot(FLDGC1,aes(x=GC,y=count))+ geom_violin() + scale_y_log10()

