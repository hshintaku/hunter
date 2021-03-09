library(dplyr)
library(Seurat)
#library(patchwork)
#library(SingleCellSignalR)

# Load the PBMC dataset
#wdir <- "/home/watson/sanger/shintaku/HUNTER/"
pbmc.data <- Read10X(data.dir = wdir)


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 600, min.features = 1000)

cellids <- colnames(pbmc)
pbmc[['plates']] <- substr(cellids,1,3)
pbmc[['dish']] <- substr(cellids,4,6)
pbmc[['gate']] <- substr(cellids,7,8)
pbmc[['pool']] <- substr(cellids,10,10)
pbmc[['rtid']] <- substr(cellids,12,13)


head(x = FetchData(object = pbmc, vars = c('ident')),n=1000)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 200)

