library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellSignalR)

# Load the PBMC dataset
datadir <- "/home/watson/sanger/shintaku/HUNTER/"
pbmc.data <- Read10X(data.dir = datadir)


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 100, min.features = 50)

#seurat_object = CreateSeuratObject(counts = B$`protein_coding`)

#pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "-MT-")
#pbmc[['percent.hs']] <- PercentageFeatureSet(pbmc, pattern = "^hs-")
#pbmc[['percent.ms']] <- PercentageFeatureSet(pbmc, pattern = "^ms-")
#pbmc[['percent.pg']] <- PercentageFeatureSet(pbmc, pattern = "^pg-")
pbmc[['nCount_hs']] <- pbmc[['percent.hs']]$percent.hs*pbmc[['nCount_RNA']]$nCount_RNA
pbmc[['nCount_ms']] <- pbmc[['percent.ms']]$percent.ms*pbmc[['nCount_RNA']]$nCount_RNA
cellids <- colnames(pbmc)
pbmc[['plates']] <- substr(cellids,1,3)
pbmc[['dish']] <- substr(cellids,4,6)
pbmc[['gate']] <- substr(cellids,7,8)
pbmc[['pool']] <- substr(cellids,10,10)
pbmc[['rtid']] <- substr(cellids,12,13)

head(x = FetchData(object = pbmc, vars = c('ident')),n=1000)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 200)






#VlnPlot(pbmc, features = c("percent.hs", "percent.ms"), ncol = 2)


#pbmc_subset <- SubsetData(object = pbmc, cells=pca_topcells)




pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:7)
JackStrawPlot(pbmc, dims = 1:7)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:7)
pbmc <- FindClusters(pbmc, resolution = 0.1)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(pbmc))
data = data.frame(pbmc[["RNA"]]@data)

pbmc <- RunUMAP(pbmc, dims = 1:4)
DimPlot(pbmc, reduction = "umap")

# cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, min.pct = 0)
# head(cluster0.markers, n = 5)
# cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0)
# head(cluster1.markers, n = 5)
# cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0)
# head(cluster2.markers, n = 5)
# cluster3.markers <- FindMarkers(pbmc, ident.1 = 3, min.pct = 0)
# head(cluster3.markers, n = 5)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.1)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

VlnPlot(pbmc, features = c("pg-GAPDH", "pg-S100A2"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("hs-MT-ND4", "pg-GAPDH"))


# Ligand/Receptor analysis using SingleCellSignalR
signal = cell_signaling(data=data,genes=all.genes,cluster=cluster)

# Visualization
visualize(signal)
intra = intra_network("S1PR1",data,all.genes,cluster,"cluster 1",signal = signal)


