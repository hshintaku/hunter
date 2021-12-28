
#library(patchwork)

# Load the PBMC dataset
#wdir <- "/home/watson/sanger/shintaku/HUNTER/"
pbmc.data <- Read10X(data.dir = wdir)


# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 1, min.features = 100)

cellids <- colnames(pbmc)

pbmc[['rtid']] <- substr(cellids,17,18)
pbmc[['pool']] <- substr(cellids,14,14)
cellmetadata <- read.xlsx("/home/samba/public/shintaku/20211124HiSeqX006_Islet/cellid_list.xlsx")
rownames(cellmetadata)<-cellmetadata$cellid
pbmc[['condition']]<-cellmetadata[cellids,]$cell
pbmc[['species']]<-cellmetadata[cellids,]$species

dim(pbmc)
rm(pbmc.data)

