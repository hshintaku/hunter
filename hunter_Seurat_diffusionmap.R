#library(remotes)
#install_github("theislab/destiny")
#https://github.com/theislab/destiny/blob/master/vignettes/Diffusion-Maps.ipynb
library(destiny)
library(ggplot2)
library(conflicted)
suppressPackageStartupMessages(library(scran))
library(purrr)
library(wordspace)

#hepa.data <- liver[["RNA"]]@data
hepa.data <- hepa[["RNA"]]@data
#hepa.data <- pbmc[["RNA"]]@data

hepa.data.zone <- data.frame(t(hepa.data[ordering_genes_disp$gene_id,]))
cellids<- rownames(hepa.data.zone)
hepa.data.zone$plate <- substr(cellids,1,3)
dm <- DiffusionMap(hepa.data.zone)
p1<-plot(dm,1:2,
     col_by = 'Cyp2e1',
     legend_main = 'Cyp2e1')
p2<-plot(dm,1:2,
         col_by = 'Cyp2f2',
         legend_main = 'Cyp2f2')
p3<-plot(dm,1:2,
         col_by = 'plate',
         legend_main = 'plate')
p1+p2+p3

qplot(y = eigenvalues(dm)) + theme_minimal() +
  labs(x = 'Diffusion component (DC)', y = 'Eigenvalue')

set.seed(1)
dms <- c('euclidean', 'cosine', 'rankcor') %>% #, 'l2'
  set_names() %>%
  map(~ DiffusionMap(hepa.data.zone, distance = ., knn_params = list(method = 'covertree')))
options(repr.plot.width = 14, repr.plot.height = 4)
dms %>%
  imap(function(dm, dist) plot(dm, 1:2,col_by='Cyp2f2') + ggtitle(dist)) %>%
  cowplot::plot_grid(plotlist = ., nrow = 1)

hepa.data.zone$norm <- dm$DC1#rowNorms(cbind(dm$DC1,dm$DC2), method = "euclidean", p = 2)
ggplot(hepa.data.zone,aes(x=norm,y=Cyp2e1,colour=plate))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Cyp2c37,colour=plate))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Cyp1a2,colour=plate))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Cyp2f2,colour=plate))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Hsd17b13,colour=plate))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Pck1,colour=plate))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Arg1,colour=plate))+geom_point()+
  ggplot(hepa.data.zone,aes(x=norm,y=Ass1,colour=plate))+geom_point()
ggplot(hepa.data.zone,aes(x=plate,y=norm,fill=plate))+geom_violin()+
  geom_jitter(shape = 16, position = position_jitter(0.07))
library(pheatmap)
hepa.heat <- hepa.data.zone[order(hepa.data.zone$norm,decreasing = FALSE),colnames(hepa.data.zone)!=c("plate","norm")]
hepa.pn <-data.frame(hepa.heat[,colnames(hepa.heat) %in% toTitleCase(genes_pn$gene_id)])
hepa.cv <-hepa.heat[,colnames(hepa.heat) %in% toTitleCase(genes_cv$gene_id)]
cbind(hepa.pn,hepa.cv)
pheatmap(hepa.cv,
                cluster_cols = TRUE,cluster_rows = FALSE,scale="column",clustering_distance_cols = "canberra")
  pheatmap(hepa.pn,
             cluster_cols = TRUE,cluster_rows = FALSE,scale="column",clustering_distance_cols = "canberra")

ggplot(hepa.pn,aes())
