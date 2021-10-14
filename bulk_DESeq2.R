
library(DESeq2)
pbmc.data <- Read10X(data.dir = wdir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 1, min.features = 100)
pbmc12 <- subset(pbmc, subset = batch =="03", invert = TRUE)

pbmc12.data <- data.frame(pbmc12[["RNA"]]@counts)
#CellScatter(pbmc,cell1="Ida-Lib-01-7",cell2="Ida-Lib-02-9")+scale_x_log10()+scale_y_log10()

pbmc.mtx <- as.matrix(pbmc12.data)
group <- data.frame(con = pbmc12[["batch"]])
#colnames(pbmc.mtx) <- group$batch

dds <- DESeqDataSetFromMatrix(countData = pbmc.mtx, colData = group, design = ~ batch)
dds <- DESeq(dds)
res <- results(dds)
head(res)
