
files <- data.frame(list.files(datadir,pattern="ReadsPerGene.out.tab", recursive=T))
colnames(files)<-"name"
for (icnt in 1:nrow(files)){
  #icnt=1
  myData <- read.table(file.path(datadir,files[icnt,]), skip=4, header = FALSE)
  myData <- myData[,c(1,2)]
  #white<-subset(allencoded,allencoded$batch==str_sub(files[icnt,],1,10))
  #encoded <- whitelist.umi_tools.encode(myData$cell,barcode$V1)
  myData$cell <- str_sub(files[icnt,],1,14)#paste0(str_sub(files[icnt,],1,10),"-",encoded$index)
  myData <- myData[myData$V2>0,]
  colnames(myData) <- c("gene","count","cell")
  if (icnt>1){
    allData <- rbind(allData,myData)
  }
  else{
    allData <- myData
  }
}
rm(myData)
