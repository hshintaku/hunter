
# show number of counts
VlnPlot(pbmc, features = c("nCount_RNA","nFeature_RNA"),
        ncol = 2,group.by = "plate")
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "gate" )
p1<-VlnPlot(pbmc, features = c("percent.mt"),group.by = "plate" )
p2<-FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "plate" )
p1+p2
pbmc <- subset(pbmc, subset= percent.mt<7)

p1<-VlnPlot(pbmc, features = c("percent.mt"),group.by = "plate" )
p2<-FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "plate" )
p1+p2
count_summary <- pbmc[[c("nCount_RNA","nFeature_RNA","gate","cell")]]
count_summary_mean <- count_summary %>%
  dplyr::group_by(gate,cell) %>%
  dplyr::summarise(nCount.mean=mean(nCount_RNA), nFeature.mean=mean(nFeature_RNA))


# number of UMI counts under various conditions
ggplot(count_summary_mean,aes(x=factor(gate),y=factor(cell),
                              color=nCount.mean,size=nCount.mean))+
  geom_point()+
  scale_size(range = c(10, 20))

tenx <- FetchData(pbmc,vars="nCount_RNA")
median(tenx$nCount_RNA)
tenx <- FetchData(pbmc,vars="nFeature_RNA")
median(tenx$nFeature_RNA)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 20)
# plot variable features with labels
plot1 <- VariableFeaturePlot(pbmc)
plot1 <- LabelPoints(plot = plot1, points = top10)
plot1
#plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" )
#plot2

# running PCA npcs in the RunPCA function must be less than the number of samples
# default is 50

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, npcs=20, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:2, nfeatures = 50)

# plot plates/dishes/gates/pools/rtid
p1 <- DimPlot(pbmc, reduction = "pca",group.by = "plate")
p2 <- DimPlot(pbmc, reduction = "pca",group.by = "cell")
p3 <- DimPlot(pbmc, reduction = "pca",group.by = "gate")
p1+p2+p3

DimPlot(pbmc, reduction = "pca",group.by = "pool")
DimPlot(pbmc, reduction = "pca",group.by = "rtid")

# gene expression scatter
pca_topcells <- TopCells(object = pbmc[['pca']], balanced = FALSE)
CellScatter(object = pbmc, cell1 = pca_topcells[1], cell2 = pca_topcells[2])


pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.8)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:5)
p1 <- DimPlot(pbmc, reduction = "pca",group.by = "plate")
p2 <- DimPlot(pbmc, reduction = "umap",group.by = "plate")
p3<-DimPlot(pbmc)
p1+p2+p3

#find marker genes in each cluster
#DimPlot(pbmc, reduction = "umap")
pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.1, logfc.threshold =0.25 )
#pbmc.markers %>% group_by(cluster) %>% top_n(n = 2)

#FeaturePlot(pbmc,features="Serpinala")

p1 <- DimPlot(pbmc, reduction = "umap",group.by = "plate")
p2<-FeaturePlot(pbmc,features="Saa1")
p3<-DimPlot(pbmc)

cluster <-as.numeric(Idents(pbmc))
cluster <- data.frame(cluster)
rownames(cluster) <- colnames(pbmc)
cluster$gate <- pbmc[['gate']]
cluster$plate <- pbmc[['plate']]


cluster_density <- cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::count(plate, name = 'count')

p4<-ggplot(cluster_density,aes(y=cluster,x=plate$plate,fill=count))+ 
  geom_tile()+ theme_classic()+
  scale_fill_gradientn(colours = c("white", "red"))+geom_text(aes(label = count)) 

p1+p2+p3+p4

rm(p1,p2,p3,plot1,plot2,count_summary,count_summary_mean,pca_topcells,top10,all.genes,tenx)

