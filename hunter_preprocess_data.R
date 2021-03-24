
# load count data
#datadir <- dir #"/home/watson/sanger/shintaku/20210216HiSeqX002/count/"
files <- data.frame(list.files(file.path(datadir,"count"),pattern="counts.tsv.gz"))
colnames(files)<-"name"

for (icnt in 1:nrow(files)){
  #icnt=1
  myData <- read.table(file.path(datadir,"count",files[icnt,]), header = TRUE)
  white<-subset(allencoded,allencoded$batch==str_sub(files[icnt,],1,10))
  data.frame(substr(whitelist$represent, 1,8))
  encoded1 <- whitelist.umi_tools.encode(substr(myData$cell,1,8),barcode$V1)
  encoded2 <- whitelist.umi_tools.encode(substr(myData$cell,9,16),barcode$V1)
  encoded3 <- whitelist.umi_tools.encode(substr(myData$cell,17,24),barcode$V1)
  myData$cell <- paste0(str_sub(files[icnt,],1,10),"-",encoded1$index,"-",encoded2$index,"-",encoded3$index)
  
  myData <- myData[encoded1$value<1,]
  
  if (icnt>1){
    allData <- rbind(allData,myData)
  }
  else{
    allData <- myData
  }
}
rm(myData)

