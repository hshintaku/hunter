library(DESeq2)
#onchip <- subset(islet,subset=condition=="minus",invert=TRUE)
#onchip <- subset(onchip,subset=condition=="control",invert=TRUE)
pbmc12.data <- islet[["RNA"]]@counts

pbmc.mtx <- as.matrix(pbmc12.data)
group <- data.frame(con = islet[["pool"]])

#colnames(pbmc.mtx) <- group$batch

dds <- DESeqDataSetFromMatrix(countData = pbmc.mtx, colData = group, design = ~ pool)
dds <- DESeq(dds)
res <- results(dds)
View(data.frame(res))
head(res)
plotMA(res, alpha = 0.01)

result <- data.frame(res)
de_genes <- subset(result, subset=padj<0.001)

library(pheatmap)
conflict_prefer("pheatmap", "pheatmap")
pheatmap(pbmc[["RNA"]]@data[rownames(pbmc) %in% rownames(de_genes),])
