
#library(patchwork)

# Load the PBMC dataset
#wdir <- "/home/watson/sanger/shintaku/HUNTER/"
fucci.data <- Read10X(data.dir = wdir)


# Initialize the Seurat object with the raw (non-normalized data).
fucci <- CreateSeuratObject(counts = fucci.data, project = "helafucci", min.cells = 1, min.features = 1)

cellids <- colnames(fucci)
fucci[['exp']] <- substr(cellids,1,9)
fucci[['batch']] <- substr(cellids,11,14)
fucci[['time']] <- substr(cellids,9,9)
#pbmc[['batch']] <- substr(cellids,6,7)
#pbmc[['pool']] <- substr(cellids,10,10)
#pbmc[['rtid']] <- substr(cellids,12,13)


head(x = FetchData(object = fucci, vars = c('ident')),n=1000)

fucci <- NormalizeData(fucci, normalization.method = "LogNormalize", scale.factor = 1e5)
fucci <- FindVariableFeatures(fucci, selection.method = "vst", nfeatures = 200)

