
# load count data
#datadir <- dir #"/home/watson/sanger/shintaku/20210216HiSeqX002/count/"
files <- data.frame(list.files(datadir,pattern="counts.tsv.gz"))
colnames(files)<-"name"

for (icnt in 1:nrow(files)){
  #icnt=1
  myData <- read.table(file.path(datadir,files[icnt,]), header = TRUE)
  white<-subset(allencoded,allencoded$batch==str_sub(files[icnt,],1,10))
  encoded <- whitelist.umi_tools.encode(myData$cell,barcode$V1)
  myData$cell <- paste0(str_sub(files[icnt,],1,10),"-",encoded$index)
  
  myData <- myData[encoded$value<1,]
  
  if (icnt>1){
    allData <- rbind(allData,myData)
  }
  else{
    allData <- myData
  }
}
rm(myData)

# download reference data from ensembl with biomaRt
source('/home/watson/public/shintaku/HUNTER/hunter_biomart_ref.R')
missing_ref <- subset(gene_list,!(gene %in% ms_ref$ensembl_gene_id))
adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
colnames(adding_ref) <- colnames(ms_ref)
rownames(adding_ref) <- adding_ref$ensembl_gene_id
ms_ref <- rbind(adding_ref,ms_ref)