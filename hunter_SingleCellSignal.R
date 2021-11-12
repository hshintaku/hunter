#
# hunter_SingleCellSignal
# http://www.bioconductor.org/packages/release/bioc/html/SingleCellSignalR.html
# browseVignettes("SingleCellSignalR")


# Ligand/Receptor analysis using SingleCellSignalR
#signal = cell_signaling(data=data,genes=all.genes,cluster=cluster)
# user guide
# https://rdrr.io/bioc/SingleCellSignalR/f/vignettes/UsersGuide.Rmd

# Visualization
#visualize(signal)
#intra = intra_network("S1PR1",data,all.genes,cluster,"cluster 1",signal = signal)

library(SingleCellSignalR)

#scdata <- data.frame(GetAssayData(object=pbmc[["RNA"]]))
#scdata <- scdata-min(min(scdata))
# Data clustering
scdata = data.frame(pbmc[["RNA"]]@data)

all.genes <- row.names(scdata)
#clust <- clustering(data=scdata, n.cluster=4, n=10,method="simlr",write=TRUE,pdf=FALSE)

#clust <-cluster
signal = cell_signaling(data=scdata,genes=all.genes,cluster=cluster,species ="mus musculus",
                        logFC=log2(1.5),s.score=0.6,int.type = "paracrine",write=TRUE,c.names=c("E0771","E0771vivo","GFP+","GFP-","mCherry+semihelth","healthy"))

inter.net <- inter_network(data = scdata, signal = signal, genes = all.genes, cluster = cluster, write = FALSE)
visualize_interactions(signal = signal)

