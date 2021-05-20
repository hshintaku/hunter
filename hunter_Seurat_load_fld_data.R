library(Seurat)
library(stringr)
library(dplyr)

wFdir <- "/home/samba/storage0/Shiomi/20210427MiSeq017Ana/"
sample_text <- "/home/samba/storage0/Shiomi/20210427MiSeq017Ana/citesample.txt"
sampleID <- read.table(sample_text) 
barcode <- read.table(file.path("/home/samba/storage0/Shiomi/20210323MiSeq015Ana10X/RTbarcodes.csv"))

FLDmapALL <- matrix(nrow=7, ncol=0)
for (i in 1:19){
  FLDmap <- matrix(nrow=7, ncol=0)
  samfol = paste(file.path(wFdir), "CITE-seq/", sampleID[i,1], "/umi_count/",sep="")
  FLD.data <- Read10X(data.dir = samfol, gene.column=1)
  iName = substr(sampleID[i,1], 1, 10)
  romin <- whitelist.umi_tools.encode(colnames(FLD.data),barcode$V1)
  FLData <- as.data.frame(as.matrix(FLD.data))
  romin1 <- romin$index * ifelse(romin$value>2, 0, 1)
  FLData <- as.data.frame(cbind(t(FLData), romin1))
  FLData <- group_by(FLData, romin1)
  FLData <- summarise_each(FLData, dplyr::funs(sum))
  FLDmap <- as.data.frame(FLDmap)
  FLDempty <- FLData[1,]
  for (t in 1:16){
    if (is.na(match(t, FLData$romin1))){
      FLDempty[1,] <- list(t,0,0,0,0,0,0)
      FLDmap <-rbind(FLDmap, FLDempty)
    } else{
      FLDmap <-rbind(FLDmap, filter(FLData, romin1 == t))
    }
  }
  rownames(FLDmap) <- paste(iName, 1:16, sep="-")
  FLDmapALL <-rbind(FLDmapALL, FLDmap)
}

FLDmapE <- as.data.frame(t(FLDmapALL[, c(1, 2, 3, 4, 5)]))
rownames(FLDmapE) <- substr(rownames(FLDmapE),1,6)
colnames(FLDmapE) <- toupper(rownames(FLDmapALL))

# adding FLD data
seladt.FLD.csv <- FLDmapE[,cellids] # extract cells detected in RNA-seq
pbmc.adt.FLD <-as.sparse(seladt.FLD.csv) # convert the format to sparse matrix
pbmc[["ADT"]] <- CreateAssayObject(counts=pbmc.adt.FLD)


# adding FLD+adt data
Total <- rbind(adt.csv[,cellids], FLDmapE[,cellids])
pbmc.adt.Total <-as.sparse(Total) # convert the format to sparse matrix
pbmc[["ADT"]] <- CreateAssayObject(counts=pbmc.adt.Total)

#pbmc <- NormalizeData(pbmc,assay="ADT",normalization.method = "LogNormalize",scale.factor = 1e5)
pbmc <- NormalizeData(pbmc,assay="ADT",normalization.method = "LogNormalize",margin=2,scale.factor = 1e5)

#NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)
pbmc <- ScaleData(pbmc,assay="ADT")

# correlation in PC_1 vs GFP
TOTAL <- subset(x=pbmc)
p1 <- FeatureScatter(object=pbmc, feature1='Venus',feature2='FLD004',group.by = 'Sort')
p2 <- FeatureScatter(object=TOTAL, feature1='Venus',feature2='FLD010',group.by = 'Sort')
p3 <- FeatureScatter(object=TOTAL, feature1='Venus',feature2='FLD070',group.by = 'Sort')
p4 <- FeatureScatter(object=TOTAL, feature1='Venus',feature2='FLD500',group.by = 'Sort')
p5 <- FeatureScatter(object=TOTAL, feature1='Venus',feature2='FSC',group.by = 'Sort')
p6 <- FeatureScatter(object=TOTAL, feature1='Venus',feature2='SSC',group.by = 'Sort')
grid.arrange(p1, p2, p3, p4,p5,p6, nrow = 2)



Total[1,] <- substr(colnames(Total),4, 8)
tTotal <- as.data.frame(t(Total))
for (t in 2:11){
tTotal[, t] <- as.numeric(tTotal[, t])
}
GC05 <- c("1", "2", "5", "7", "8", "9", "11", "12", "13", "14", "15")
GC02 <- c("4", "6", "16")
GC03 <- c("3", "10")
Live <- c("NN", "RG")
tTotal1 <- subset(tTotal, !(is.na(tTotal$normGFP)))
tTotal1 <- filter(tTotal1, str_detect(romin1, paste(GC05, collapse = "|")))
tTotal1 <- filter(tTotal1, str_detect(Events, paste(Live, collapse = "|")))


p_0 <- ggplot(data = tTotal1, mapping = aes(x = FLD004,
                                           y =  mCherry, 
                                           color = Events))
p_0 <- p_0 + scale_x_log10()
p_0 <- p_0 + scale_y_log10()
p_0 <- p_0 +geom_point()
p_0
ggsave(file = "GFP_NoemalizeFLD.pdf", plot = p_0, dpi = 100, width = 5, height = 5)


