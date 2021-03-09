#
# hunter_SingleCellSignal
# http://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html
# browseVignettes("SingleCellSignalR")


# Ligand/Receptor analysis using SingleCellSignalR
#signal = cell_signaling(data=data,genes=all.genes,cluster=cluster)

# Visualization
#visualize(signal)
#intra = intra_network("S1PR1",data,all.genes,cluster,"cluster 1",signal = signal)

library(SingleCellSignalR)

#scdata <- data.frame(GetAssayData(object=pbmc[["RNA"]]))
#scdata <- scdata-min(min(scdata))
# Data clustering
#all.genes <- rownames(pbmc)
pbmc <- NormalizeData(pbmc,scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.1)

cluster = as.numeric(Idents(pbmc))
scdata = data.frame(pbmc[["RNA"]]@data)

all.genes <- row.names(scdata)
clust <- clustering(data=scdata, n.cluster=4, n=10,method="simlr",write=TRUE,pdf=FALSE)

signal = cell_signaling(data=scdata,genes=all.genes,cluster=clust$cluster,species ="mus musculus",
                        logFC=log2(1.0001),s.score=0.0001,int.type = "paracrine",write=TRUE)


#clust.analysis <- cluster_analysis(data=scdata,genes=rownames(scdata),cluster=clust$cluster,write=FALSE)


signal <- cell_signaling(data = data, genes = all.genes, 
                         cluster = cluster, species ="mus musculus",write = FALSE,
                         logFC=log2(1.0001),s.score=0.001)


sc_names <- data.frame(rownames(scdata)) # list gene short names
scnames_ref <- ms_ref[match(sc_names$rownames.scdata.,ms_ref$gene_short_name),] # create reference from ms_ref regardless of entrez annotation

  
#sc_annot_index <- data.frame(which(!is.na(datExpr_ref$entrez_annotation))) #
#ms_ref_subset <- datExpr_ref[which(!is.na(datExpr_ref$entrez_annotation)),]#ms_ref[which(!is.na(ms_ref$entrez_annotation)),]


t_sne <- data.frame(clust$`t-SNE`)
t_sne$dish <- pbmc[['dish']]
colnames(t_sne) <- c("tsne1","tsne2","dish")
ggplot(t_sne,aes(x=tsne1,y=tsne2,shape=dish))+geom_point()

