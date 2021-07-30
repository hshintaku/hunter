
#library(patchwork)

# Load the PBMC dataset
#wdir <- "/home/watson/sanger/shintaku/HUNTER/"
pbmc.data <- Read10X(data.dir = wdir)


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 1, min.features = 1)

cellids <- colnames(pbmc)
pbmc[['plates']] <- substr(cellids,1,3)
pbmc[['treat']] <- substr(cellids,4,5)
pbmc[['cell']] <- substr(cellids,4,6)
pbmc[['gate']] <- substr(cellids,7,8)
pbmc[['pool']] <- substr(cellids,10,10)
pbmc[['rtid']] <- substr(cellids,12,13)

dim(pbmc)
#head(x = FetchData(object = pbmc, vars = c('ident')),n=1000)
rm(pbmc.data)

