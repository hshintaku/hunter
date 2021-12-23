
# show number of counts
pbmc <- subset(pbmc, subset = percent.mt <5)
VlnPlot(pbmc, features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol = 3)
pbmc <- subset(pbmc, subset = percent.mt < 5)
tenx <- FetchData(pbmc,vars="nCount_RNA")
median(tenx$nCount_RNA)
tenx <- FetchData(pbmc,vars="nFeature_RNA")
median(tenx$nFeature_RNA)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 40)
# plot variable features with labels
p1 <- VariableFeaturePlot(pbmc)
p1 <- LabelPoints(plot = plot1, points = top10)
p1
p2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" )
p2

