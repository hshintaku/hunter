
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


