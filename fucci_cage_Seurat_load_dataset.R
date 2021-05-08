
#library(patchwork)

# Load the PBMC dataset
#wdir <- "/home/watson/sanger/shintaku/HUNTER/"
fucci.data <- Read10X(data.dir = wdir)


# Initialize the Seurat object with the raw (non-normalized data).
fucci_cage <- CreateSeuratObject(counts = fucci.data, project = "helafucci", min.cells = 1, min.features = 1)

cellids <- colnames(fucci_cage)
fucci_cage[['exp']] <- substr(cellids,1,9)
fucci_cage[['batch']] <- substr(cellids,11,14)
fucci_cage[['time']] <- substr(cellids,9,9)
#pbmc[['batch']] <- substr(cellids,6,7)
#pbmc[['pool']] <- substr(cellids,10,10)
#pbmc[['rtid']] <- substr(cellids,12,13)


head(x = FetchData(object = fucci_cage, vars = c('ident')),n=1000)

fucci_cage <- NormalizeData(fucci_cage, normalization.method = "LogNormalize", scale.factor = 1e5)
fucci_cage <- FindVariableFeatures(fucci_cage, selection.method = "vst", nfeatures = 200)

