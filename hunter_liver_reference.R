#https://tabula-muris.ds.czbiohub.org/
library(Seurat)
liver.data <- Read10X(data.dir = "/home/samba/public/tabula_muris/droplet/Liver-10X_P7_0/")
liver <- CreateSeuratObject(counts = liver.data, project = "liver", min.cells = 10, min.features = 1000)
liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = 1e5)
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 500)
VlnPlot(liver, features = c("nCount_RNA","nFeature_RNA"))
FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
top10 <- head(VariableFeatures(liver), 200)
# plot variable features with labels
plot1 <- VariableFeaturePlot(liver)
plot1 <- LabelPoints(plot = plot1, points = top10)
plot1

all.genes <- rownames(liver)
liver <- ScaleData(liver, features = all.genes)
liver <- RunPCA(liver, npcs=20, features = VariableFeatures(object = liver))

# plot plates/dishes/gates/pools/rtid
DimPlot(liver, reduction = "pca")

liver <- JackStraw(liver, num.replicate = 100)
liver <- ScoreJackStraw(liver, dims = 1:20)
JackStrawPlot(liver, dims = 1:20)
ElbowPlot(liver)

liver <- FindNeighbors(liver, dims = 1:9)
liver <- FindClusters(liver, resolution = 0.4)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(liver))
liver <- RunUMAP(liver, dims = 1:9)
p1 <- DimPlot(liver, reduction = "pca")
p2 <- DimPlot(liver, reduction = "umap")
p1+p2

liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.3)
liver.markers %>% group_by(cluster) %>% top_n(n = 2)

cellannotation <-read.csv2("/home/samba/public/tabula_muris/annotations_droplets.csv")
