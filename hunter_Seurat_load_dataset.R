library(dplyr)
library(Seurat)
#library(patchwork)
#library(SingleCellSignalR)

# Load the PBMC dataset
#wdir <- "/home/watson/sanger/shintaku/HUNTER/"
pbmc.data <- Read10X(data.dir = wdir)


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 100, min.features = 50)

cellids <- colnames(pbmc)
pbmc[['plates']] <- substr(cellids,1,3)
pbmc[['dish']] <- substr(cellids,4,6)
pbmc[['gate']] <- substr(cellids,7,8)
pbmc[['pool']] <- substr(cellids,10,10)
pbmc[['rtid']] <- substr(cellids,12,13)


head(x = FetchData(object = pbmc, vars = c('ident')),n=1000)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 200)



#VlnPlot(pbmc, features = c("pg-GAPDH", "pg-S100A2"), slot = "counts", log = TRUE)
#FeaturePlot(pbmc, features = c("hs-MT-ND4", "pg-GAPDH"))


# Ligand/Receptor analysis using SingleCellSignalR
#signal = cell_signaling(data=data,genes=all.genes,cluster=cluster)

# Visualization
#visualize(signal)
#intra = intra_network("S1PR1",data,all.genes,cluster,"cluster 1",signal = signal)


