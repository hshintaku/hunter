load.adt <- function (indexfiles,batch) {
  adt.xlsx <- read.xlsx(indexfiles[1],sheet=1, rows=16:111,cols=1:7,colNames = FALSE,rowNames = TRUE)
  colnames(adt.xlsx)<- c("Events","FSC","SSC","Venus","APC","mCherry")
  
  # create normalized GFP
  adt.xlsx$normGFP <- (adt.xlsx$Venus/adt.xlsx$mCherry)*1e3
  
  
  wellid <- rownames(adt.xlsx)
  for (ij in 1:length(wellid)){
    icnt=as.numeric(which(LETTERS==substr(wellid[ij],1,1)))
    jcnt=as.numeric(substr(wellid[ij],2,3))
    pindex= (jcnt+1)%/%2
    rtindex <- icnt+((jcnt-1)%%2)*8
    adt.xlsx$cell[ij] <- paste0(batch,"-",pindex,"-",rtindex)
  }
  rownames(adt.xlsx)<-adt.xlsx$cell
  adt.xlsx <- t(adt.xlsx)
  adt.xlsx <- adt.xlsx[-nrow(adt.xlsx),]
  return(adt.xlsx)
}

#load FACS index data file name is 8 char long.
batch <- unique(substr(files$name,1,8))
indexfiles <- file.path(indexdir,paste0(batch,".xlsx"))

#indexlist <- rbind(indexfiles,batch)

adt.xlsx <- mapply(load.adt,indexfiles,batch)
adt.csv <- do.call("cbind",adt.xlsx)

# adding adt data
seladt.csv <- adt.csv[,cellids] # extract cells detected in RNA-seq
pbmc.adt <-as.sparse(seladt.csv) # convert the format to sparse matrix

pbmc[["ADT"]] <- CreateAssayObject(counts=pbmc.adt)

#pbmc <- NormalizeData(pbmc,assay="ADT",normalization.method = "LogNormalize",scale.factor = 1e5)
pbmc <- NormalizeData(pbmc,assay="ADT",normalization.method = "LogNormalize",margin=1)
#NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)

pbmc <- ScaleData(pbmc,assay="ADT")

