library(DESeq2)

piezo[["genotype"]] <- c("ko4ohm","dmso","ko4ohu","dmso","","test","test","test","test2","test2","test2","test2","test2","test2","test2","test2")

pbmc12 <- subset(piezo,subset=genotype=="test",invert=TRUE)
test1 <- subset(piezo,subset=genotype=="test1")
test2 <- subset(piezo,subset=genotype=="test2")
all <- merge(test1,y=test2, project = "hirano")



pbmc12.data <- pbmc12[["RNA"]]@counts
#CellScatter(pbmc,cell1="Ida-Lib-01-7",cell2="Ida-Lib-02-9")+scale_x_log10()+scale_y_log10()

pbmc.mtx <- as.matrix(pbmc12.data)
group <- data.frame(con = pbmc12[["genotype"]])

#colnames(pbmc.mtx) <- group$batch

dds <- DESeqDataSetFromMatrix(countData = pbmc.mtx, colData = group, design = ~ genotype)
dds <- DESeq(dds)
res <- results(dds)
View(data.frame(res))
head(res)
plotMA(res, alpha = 1)