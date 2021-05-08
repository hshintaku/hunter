
#library(patchwork)

# Load the PBMC dataset
#wdir <- "/home/watson/sanger/shintaku/HUNTER/"
fucci.data <- Read10X(data.dir = wdir)


# Initialize the Seurat object with the raw (non-normalized data).
fucci <- CreateSeuratObject(counts = fucci.data, project = "elp", min.cells = 1, min.features = 1)

cellids <- colnames(fucci)
fucci[['cell']] <- substr(cellids,1,3)
fucci[['exp']] <- substr(cellids,4,5)
fucci[['batch']] <- substr(cellids,6,7)
#pbmc[['pool']] <- substr(cellids,10,10)
#pbmc[['rtid']] <- substr(cellids,12,13)


head(x = FetchData(object = fucci, vars = c('ident')),n=1000)

fucci <- NormalizeData(fucci, normalization.method = "LogNormalize", scale.factor = 1e5)
fucci <- FindVariableFeatures(fucci, selection.method = "vst", nfeatures = 200)

