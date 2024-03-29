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

library(Seurat)
library(SingleCellSignalR)

#scdata <- data.frame(GetAssayData(object=pbmc[["RNA"]]))
#scdata <- scdata-min(min(scdata))
# Data clustering
pbmc<-hepa
p1<- DimPlot(pbmc)
p2<- DimPlot(pbmc,group.by = "plate")
p1+p2
#pbmc <- allcell

scdata = data.frame(pbmc[["RNA"]]@data)

all.genes <- row.names(scdata)
#clust <- clustering(data=scdata, n.cluster=4, n=10,method="simlr",write=TRUE,pdf=FALSE)

#clust <-cluster
signal = cell_signaling(data=scdata,genes=all.genes,cluster=cluster$cluster,species ="mus musculus",
                        logFC=log2(2),s.score=0.4,int.type = "paracrine",write=TRUE)
#                        c.names=levels(cell_all))

inter.net <- inter_network(data = scdata, signal = signal, genes = all.genes, cluster = cluster$cluster, write = FALSE)
visualize_interactions(signal = signal)


#cellinker <-read.table("/home/samba/public/genome/GRCm38-mouse/cellinker_mus_musculus.txt",
#                       header = TRUE, sep="\t",fill=TRUE)
#hunter.interaction <- signal$`mCherry+-E0771vivo`
#mCherry.ligand <- cellinker[cellinker$Ligand_symbol %in% hunter.interaction$`mCherry+`, ]
#mCherry.receptor <- cellinker[cellinker$Ligand_symbol %in% hunter.interaction$`mCherry+`, ]
