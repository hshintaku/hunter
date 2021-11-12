#https://tabula-muris.ds.czbiohub.org/
library(Seurat)
liver.data <- Read10X(data.dir = "/home/samba/public/tabula_muris/droplet/")
liver <- CreateSeuratObject(counts = liver.data, project = "liver", min.cells = 10, min.features = 1000)
lung.data <- Read10X(data.dir = "/home/samba/public/tabula_muris/droplet/Lung-10X_P7_8/")
lung <- CreateSeuratObject(counts = lung.data, project = "liver", min.cells = 10, min.features = 1000)
liver <- merge(liver, y = lung, add.cell.ids = c("liver", "lung"), project = "hunter")
kidney.data <- Read10X(data.dir = "/home/samba/public/tabula_muris/droplet/Kidney-10X_P4_5/")
kidney <- CreateSeuratObject(counts = kidney.data, project = "kideny", min.cells = 10, min.features = 1000)
liver<-merge(liver, y = kidney, add.cell.ids = c("liver", "kideny"), project = "hunter")

liver <- NormalizeData(liver, normalization.method = "LogNormalize", scale.factor = 1e5)
liver <- FindVariableFeatures(liver, selection.method = "vst", nfeatures = 100)
VlnPlot(liver, features = c("nCount_RNA","nFeature_RNA"))
FeatureScatter(liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
top10 <- head(VariableFeatures(liver), 20)
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

liver <- FindNeighbors(liver, dims = 1:2)
liver <- FindClusters(liver, resolution = 0.1)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(liver))
liver <- RunUMAP(liver, dims = 1:9)
p1 <- DimPlot(liver, reduction = "pca")

p1

liver.markers <- FindAllMarkers(liver, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.4)
liver.markers <- liver.markers[liver.markers$p_val_adj <0.001 & liver.markers$cluster==2 & liver.markers$avg_log2FC>5,] %>% group_by(cluster)

cellannotation <-read.csv2("/home/samba/public/tabula_muris/annotations_droplets.csv")
p0<-DimPlot(liver,reduction="pca")
p1<-FeaturePlot(liver,features="Cyp2c29",reduction="pca")
p2<-FeaturePlot(liver,features="Azgp1",reduction="pca")
p0+p1+p2
