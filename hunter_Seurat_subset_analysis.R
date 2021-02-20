
#create AML object
AML <- subset(x=pbmc, subset=gate==c("g1","g4"))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(AML), 40)
plot1 <-VariableFeaturePlot(AML)
plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2 <- FeatureScatter(AML, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
grid.arrange(plot2,plot1, nrow = 1)

#RUN PCA
AML <- RunPCA(AML, npcs=10, features = VariableFeatures(object = AML))
print(AML[["pca"]], dims = 1:2, nfeatures = 50)
DimHeatmap(AML, dims = 1:10, cells = 500, balanced = TRUE)

AML <- JackStraw(AML, num.replicate = 10)
AML <- ScoreJackStraw(AML, dims = 0:2)
JackStrawPlot(AML, dims = 0:2)
ElbowPlot(AML)

AML <- FindNeighbors(AML, dims = 1:5)
AML <- FindClusters(AML, resolution = 0.1)
AML <- RunUMAP(AML, dims = 1:5)

p1 <- DimPlot(AML, reduction = "pca",group.by = "dish")
p2 <- DimPlot(AML, reduction = "umap",group.by = "dish")
p1+p2





FeatureScatter(object=AML,feature1="Venus",feature2="normGFP")

p1 <- FeatureScatter(object=AML,feature2='Venus',feature1='CS2nGFPT2AH2BmCherG01')
p2 <- FeatureScatter(object=AML,feature2='mCherry',feature1='CS2nGFPT2AH2BmCherG01')
p3 <- FeatureScatter(object=AML,feature2='normGFP',feature1='CS2nGFPT2AH2BmCherG01')

p4 <- FeatureScatter(object=AML,feature2='APC',feature1='CS2nGFPT2AH2BmCherG01')
p5 <- FeatureScatter(object=AML,feature2='FSC',feature1='CS2nGFPT2AH2BmCherG01')
p6 <- FeatureScatter(object=AML,feature2='SSC',feature1='CS2nGFPT2AH2BmCherG01')


grid.arrange(p1, p2, p3, p4,p5,p6, nrow = 2)
