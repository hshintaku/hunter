#library(remotes)
#install_github("theislab/destiny")
#https://github.com/theislab/destiny/blob/master/vignettes/Diffusion-Maps.ipynb
library(destiny)
library(ggplot2)
library(conflicted)
suppressPackageStartupMessages(library(scran))
library(purrr)
library(wordspace)
library(pheatmap)
library(R.matlab)
library(tools)
genes_zonation <- readMat("/home/samba/public/shintaku/matlab_code/Zonation_params.mat")
genes_cv <- data.frame(unlist(genes_zonation$genes.cv))
colnames(genes_cv)<-"gene_id"
rownames(genes_cv)<- genes_cv$gene_id
genes_pn <- data.frame(unlist(genes_zonation$genes.pn))
colnames(genes_pn)<-"gene_id"
rownames(genes_pn)<-genes_pn$gene_id

ordering_genes <- rbind(genes_pn,genes_cv)

#hepa.data <- liver[["RNA"]]@data
#hepa.data <- hepa[["RNA"]]@data
tabula <- Read10X(data.dir = "/home/samba/public/shintaku/tabula_muris/droplet/Liver-10X_P4_2/")
tabulamuris <- CreateSeuratObject(counts = tabula, project = "pbmc3k", min.cells = 1, min.features = 1000)
#tabula <- hepa10x01P
tabulamuris <- NormalizeData(tabulamuris, normalization.method = "LogNormalize", scale.factor = 1e5)
hepa.data <- data.frame(tabulamuris[["RNA"]]@data)

hepa.data.zone <- data.frame(t(hepa.data[rownames(hepa.data) %in% capitalize(ordering_genes$gene_id),]))
cellids<- rownames(hepa.data.zone)
#hepa.data.zone$plate <- substr(cellids,1,3)
dm <- DiffusionMap(hepa.data.zone)
p1<-plot(dm,1:2,
     col_by = 'Cyp2e1',
     legend_main = 'Cyp2e1')
p2<-plot(dm,1:2,
         col_by = 'Cyp2f2',
         legend_main = 'Cyp2f2')
# p3<-plot(dm,1:2,
#          col_by = 'plate',
#          legend_main = 'plate')
p1+p2#+p3

# qplot(y = eigenvalues(dm)) + theme_minimal() +
#   labs(x = 'Diffusion component (DC)', y = 'Eigenvalue')
# 
# set.seed(1)
# dms <- c('euclidean', 'cosine', 'rankcor') %>% #, 'l2'
#   set_names() %>%
#   map(~ DiffusionMap(hepa.data.zone, distance = ., knn_params = list(method = 'covertree')))
# options(repr.plot.width = 14, repr.plot.height = 4)
# dms %>%
#   imap(function(dm, dist) plot(dm, 1:2,col_by='Cyp2f2') + ggtitle(dist)) %>%
#   cowplot::plot_grid(plotlist = ., nrow = 1)

hepa.data.zone$norm <- dm$DC1#rowNorms(cbind(dm$DC1,dm$DC2), method = "euclidean", p = 2)
ggplot(hepa.data.zone,aes(x=norm,y=Cyp2e1))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Cyp2c37))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Cyp1a2))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Cyp2f2))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Hsd17b13))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Pck1))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Arg1))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Ass1))+geom_point()
ggplot(hepa.data.zone,aes(x=plate,y=norm,fill=plate))+geom_violin()+
  geom_jitter(shape = 16, position = position_jitter(0.07))
hepa.heat <- hepa.data.zone[order(hepa.data.zone$norm,decreasing = FALSE),colnames(hepa.data.zone)!=c("plate","norm")]
hepa.pn <-data.frame(hepa.heat[,colnames(hepa.heat) %in% toTitleCase(genes_pn$gene_id)])
hepa.cv <-hepa.heat[,colnames(hepa.heat) %in% toTitleCase(genes_cv$gene_id)]
cbind(hepa.pn,hepa.cv)

pheatmap(hepa.cv,
                cluster_cols = TRUE,cluster_rows = FALSE,scale="column",clustering_distance_cols = "canberra")
pheatmap(hepa.pn,
             cluster_cols = TRUE,cluster_rows = FALSE,scale="column",clustering_distance_cols = "canberra")

