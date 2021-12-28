library(DESeq2)
onchip <- subset(pbmc,subset=condition=="before",invert=TRUE)
onchip <- subset(onchip,subset=condition=="control",invert=TRUE)
pbmc12.data <- onchip[["RNA"]]@counts
#CellScatter(pbmc,cell1="Ida-Lib-01-7",cell2="Ida-Lib-02-9")+scale_x_log10()+scale_y_log10()

pbmc.mtx <- as.matrix(pbmc12.data)
group <- data.frame(con = onchip[["condition"]])

#colnames(pbmc.mtx) <- group$batch

dds <- DESeqDataSetFromMatrix(countData = pbmc.mtx, colData = group, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
View(data.frame(res))
head(res)
plotMA(res, alpha = 0.05)
