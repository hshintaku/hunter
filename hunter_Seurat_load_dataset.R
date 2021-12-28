
#library(patchwork)

# Load the PBMC dataset
#wdir <- "/home/watson/sanger/shintaku/HUNTER/"
pbmc.data <- Read10X(data.dir = wdir)


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 1, min.features = 100)

cellids <- colnames(pbmc)
#pbmc[['lib']] <- substr(cellids,1,8)
#pbmc[['batch']] <- substr(cellids,9,10)
#pbmc[['rtid']] <- substr(cellids,12,13)
#pbmc[['pool']] <- substr(cellids,10,10)
#pbmc[['rtid']] <- substr(cellids,12,13)

dim(pbmc)
#head(x = FetchData(object = pbmc, vars = c('ident')),n=1000)
rm(pbmc.data)

