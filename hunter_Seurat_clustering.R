
# running PCA npcs in the RunPCA function must be less than the number of samples
# default is 50

all.genes <- rownames(pbmc)
#all.genes <- ordering_genes_disp$gene_id
#pbmc<-hepa
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, npcs=20, features = VariableFeatures(object = pbmc))

print(pbmc[["pca"]], dims = 1:2, nfeatures = 50)

# plot plates/dishes/gates/pools/rtid
p1 <- DimPlot(pbmc, reduction = "pca",group.by = "plate")
p2 <- DimPlot(pbmc, reduction = "pca",group.by = "cell")
p3 <- DimPlot(pbmc, reduction = "pca",group.by = "gate")
p4<-DimPlot(pbmc, reduction = "pca",group.by = "pool")
p5<-DimPlot(pbmc, reduction = "pca",group.by = "rtid")
p1+p2+p3+p4+p5

# gene expression scatter
pca_topcells <- TopCells(object = pbmc[['pca']], balanced = FALSE)
CellScatter(object = pbmc, cell1 = pca_topcells[1], cell2 = pca_topcells[2])


pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:6)
pbmc <- FindClusters(pbmc, resolution = 0.8)


# Retreiving the results of the preprocessing from the Seurat object
cluster = as.numeric(Idents(pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:5)
p1 <- DimPlot(pbmc, reduction = "pca",group.by = "plate")
p2 <- DimPlot(pbmc, reduction = "umap",group.by = "plate")
p3<-DimPlot(pbmc)
p1+p2+p3

#find marker genes in each cluster
pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE, min.pct = 0.1, logfc.threshold =0.25 )

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
