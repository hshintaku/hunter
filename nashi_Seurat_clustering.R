
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with labels
plot1 <- VariableFeaturePlot(pbmc)
plot1 <- LabelPoints(plot = plot1, points = top10)
plot1

# running PCA npcs in the RunPCA function must be less than the number of samples
# default is 50

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, npcs=10, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:2, nfeatures = 50)
DimPlot(pbmc)
DimPlot(pbmc,group.by = "species")+DimPlot(pbmc,group.by = "condition")

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:10)
JackStrawPlot(pbmc, dims = 1:10)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:2)
pbmc <- FindClusters(pbmc, resolution = 0.3)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:2)
p1 <- DimPlot(pbmc, reduction = "pca")
p2 <- DimPlot(pbmc, reduction = "umap")
p1+p2

#find marker genes in each cluster
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0, logfc.threshold = 0.1)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 5)

rm(p1,p2,p3,plot1,plot2,count_summary,count_summary_mean,pca_topcells,top10,all.genes,tenx)


