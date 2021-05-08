library(RCurl)
library(pheatmap)
library(Seurat)
library(monocle3)
library(ggplot2)
library(reshape2)
library(dplyr)






#bulk <- subset(x=pbmc, subset=exp==c("EL"))
cellcycle_analysis <- function(bulk, list_of_times){
  #load cell cycle genes
  source(file.path(rdir,"fucci_cellcycle_genes.R"))
  
  sub_ref <- all_ref %>%
    dplyr::filter(gene_short_name %in% rownames(bulk))
  genes <- fucci_cellcycle_genes(sub_ref)
  cell_cycle_markers<-genes[[1]]
  s_genes <- genes[[2]]
  g2m_genes <- genes[[3]]
  
  bulk <- CellCycleScoring(bulk,g2m.features = g2m_genes,s.features = s_genes)
  bulk <- ScaleData(bulk, vars.to.regress = c("S.Score", "G2M.Score"), features = c(s_genes, g2m_genes))

cell_cycle_result <- cbind(bulk[["Phase"]],bulk[["S.Score"]],bulk[["G2M.Score"]],bulk[["exp"]],bulk[["batch"]],
                           bulk[["nCount_RNA"]],bulk[["nFeature_RNA"]],bulk[["RNA"]]@data["ACTB",],bulk[["RNA"]]@data["GAPDH",],
                           bulk[["RNA"]]@data["SKP2",],
                           bulk[["RNA"]]@data["FZR1",],
                           bulk[["RNA"]]@data["CDH1",],
                           bulk[["RNA"]]@data["CCNA2",],
                           bulk[["RNA"]]@data["CCNB1",])


cellids <- colnames(bulk)
cell_cycle_result <- data.frame(cell_cycle_result[sort(cellids,index=T)$ix,])
colnames(cell_cycle_result) <- c("Phase","S.Score","G2M.Score","exp","batch",
                                 "nCount_RNA","nFeature_RNA","ACTB","GAPDH",
                                 "SKP2","FZR1","CDH1","CCNA2","CCNB1")


#both
cycle_gene_expression <- t(as.data.frame(bulk[["RNA"]]@data[c(cell_cycle_markers$gene_short_name),]))
cycle_gene_expression <- cycle_gene_expression[sort(cellids,index=T)$ix,]
cell_cycle_result_EL <- cbind(cell_cycle_result,cycle_gene_expression)
cell_cycle_result_EL$hours <- list_of_times
cell_cycle_result_EL <- subset(cell_cycle_result_EL, exp !="EL")

phased_genes <- data.frame(cell_cycle_markers[,1])
rownames(phased_genes) <- cell_cycle_markers[,5]

phase_result <- data.frame(cell_cycle_result_EL[,c(1,2,3,10)])
rownames(phase_result) <- rownames(cell_cycle_result_EL)
colnames(phase_result) <- c("phase","S.Score","G2M.Score","time")

out <- pheatmap(cell_cycle_result_EL[,c("ACTB","GAPDH", c(cell_cycle_markers$gene_short_name))],
         annotation_row=phase_result,annotation_col = phased_genes,
         cluster_cols = TRUE,cluster_rows = FALSE,scale="column")
#both
# cycle_gene_expression <- t(as.data.frame(bulk[["RNA"]]@data[c(cell_cycle_markers$gene_short_name),]))
# cycle_gene_expression <- cycle_gene_expression[sort(cellids,index=T)$ix,]
# cell_cycle_result_EL <- cbind(cell_cycle_result,cycle_gene_expression)
# cell_cycle_result_EL$hours <- list_of_times
# cell_cycle_result_EL <- subset(cell_cycle_result_EL, exp =="EL")
# 
# phased_genes <- data.frame(cell_cycle_markers[,1])
# rownames(phased_genes) <- cell_cycle_markers[,5]
# 
# phase_result <- data.frame(cell_cycle_result_EL[,c(1,2,3,10)])
# rownames(phase_result) <- rownames(cell_cycle_result_EL)
# colnames(phase_result) <- c("phase","S.Score","G2M.Score","time")
# cell_cycle_result_ordered <- cell_cycle_result_EL[,c("ACTB","GAPDH", c(cell_cycle_markers$gene_short_name))]
# cell_cycle_result_ordered <- cell_cycle_result_ordered[,out$tree_col[["order"]]]
# pheatmap(cell_cycle_result_ordered,
#          annotation_row=phase_result,annotation_col = phased_genes,
#          cluster_cols = FALSE,cluster_rows = FALSE,scale="column")
# clear pattern with scale="row" using Fucci3.2 cage data.
#cell_cycle_result_EL_exp <- melt(t(cell_cycle_result_EL[,c("ACTB","GAPDH", g2m_genes)]))
#bulk[["hour"]] <- c(-1,0,4,10,16,20,24,-2,-3,-4,0,0,4,4,17,17,21,21,25,25,0,0,4,4,17,17,21,21,25,25)
cycle_gene_expression <- t(as.data.frame(bulk[["RNA"]]@data[c(cell_cycle_markers$gene_short_name),]))
cycle_gene_expression <- cycle_gene_expression[sort(cellids,index=T)$ix,]
cell_cycle_result_EL <- cbind(cell_cycle_result,cycle_gene_expression)
cell_cycle_result_EL$hour <- list_of_times
genelist <- c("ACTB","GAPDH", c(cell_cycle_markers$gene_short_name))
genelist <- genelist[out$tree_col[["order"]]]
return(cell_cycle_result_EL)


}




list_of_times_cage<- c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7)
cell_cycle_result <- cellcycle_analysis(fucci_cage,list_of_times_cage)

list_of_times<- c(-1,0,4,10,16,20,24,-2,-3,-4,0,0,4,4,17,17,21,21,25,25,0,0,4,4,17,17,21,21,25,25)
cell_cycle_result <- cellcycle_analysis(fucci,list_of_times)

fucci_elp <- subset(x=fucci, subset=exp==c("EL"))
list_of_times_elp<- c(0,0,4,4,17,17,21,21,25,25,0,0,4,4,17,17,21,21,25,25)
cell_cycle_result <- cellcycle_analysis(fucci_elp,list_of_times_elp)

fucci_elp1 <- subset(x=fucci_elp, subset=batch==c("01"))
list_of_times_elp1<- c(0,0,4,4,17,17,21,21,25,25)

list_of_times <- c(list_of_times_cage,list_of_times_elp)
fucci.combined <- merge(fucci_cage, y = fucci_elp, add.cell.ids = c("cage", "el"), project = "cellcycle")
cell_cycle_result <- cellcycle_analysis(fucci.combined,list_of_times)
#cell_cycle_result <- unlist(cell_cycle_result)
#ordered_genes <- cell_cycle_result[[1]]
#cell_cycle_result <- cell_cycle_result[[1]]

#Kaneko
cell_cycle_result_EL <- subset(cell_cycle_result, exp=="EL")

rownames(cell_cycle_result_EL) <- paste0("S", seq(from = 1, to = 20, by = 1),"-",
                                         cell_cycle_result_EL$exp,cell_cycle_result_EL$batch,"-H",cell_cycle_result_EL$hours,"-",
                                         cell_cycle_result_EL$Phase)

cell_cycle_result[cell_cycle_result$exp!="EL",]$exp <- "cage"
cell_cycle_result[cell_cycle_result$exp=="cage" & cell_cycle_result$hour==1,]$hour=5.076029142-8
  cell_cycle_result[cell_cycle_result$exp=="cage" & cell_cycle_result$hour==2,]$hour=3.195500002-8
  cell_cycle_result[cell_cycle_result$exp=="cage" & cell_cycle_result$hour==3,]$hour=5.429744736-8
  cell_cycle_result[cell_cycle_result$exp=="cage" & cell_cycle_result$hour==4,]$hour=7.350208855-8
  cell_cycle_result[cell_cycle_result$exp=="cage" & cell_cycle_result$hour==5,]$hour=9.905043657-8
  cell_cycle_result[cell_cycle_result$exp=="cage" & cell_cycle_result$hour==6,]$hour=12.7410301-8
  cell_cycle_result[cell_cycle_result$exp=="cage" & cell_cycle_result$hour==7,]$hour=18.20921783-8

cell_cycle_result_EL$sample <-rownames(cell_cycle_result_EL)

# number of detected genes/UMI
cell_cycle_result_melt <-
  melt(cell_cycle_result_EL,id.vars=c("sample","hour","batch","exp"),
       measure.vars = c("nFeature_RNA"))
P1<- ggplot(cell_cycle_result_melt,aes(x=factor(hour),y=value,color=exp))+
  geom_boxplot()+scale_y_log10()+
  xlab("time (hour)")+ylab("Number of detected genes")
cell_cycle_result_melt <-
  melt(cell_cycle_result_EL,id.vars=c("sample","hour","batch","exp"),
       measure.vars = c("nCount_RNA"))
P2 <- ggplot(cell_cycle_result_melt,aes(x=factor(hour),y=value,color=exp))+
  geom_boxplot()+scale_y_log10()+
  xlab("time (hour)")+ylab("UMI count")
P1+P2



#cell_cycle_result_EL <- cell_cycle_result

cell_cycle_result_melt <-
  melt(cell_cycle_result_EL,id.vars=c("hour","exp"),
       measure.vars = c("ACTB"))
ggplot(cell_cycle_result_melt,aes(x=factor(hour),y=value,color=exp))+
  geom_boxplot()+ylim(0,7)
  xlab("time (hour)")+ylab("ACTB expression")

  
  cell_cycle_result_melt <-
    melt(cell_cycle_result_EL,id.vars=c("hour","exp"),
         measure.vars = c("G2M.Score"))
  P1<- ggplot(cell_cycle_result_melt,aes(x=factor(hour),y=value,fill="variable"))+
    geom_boxplot()+
    xlab("time (hour)")+ylab("G2M.Score")+
    scale_fill_manual(values=c("variable"="green"))
  
  cell_cycle_result_melt <-
    melt(cell_cycle_result_EL,id.vars=c("hour","exp"),
         measure.vars = c("TOP2A","ACTB"))
  P2<- ggplot(cell_cycle_result_melt,aes(x=factor(hour),y=value,fill=variable))+
    geom_boxplot()+ylim(0,7)+
    xlab("time (hour)")+ylab("TOP2A expression")+
    scale_fill_manual(values=c(variable="green","gray"))
  
  
cell_cycle_result_melt <-
  melt(cell_cycle_result_EL,id.vars=c("hour","exp"),
       measure.vars = c("S.Score"))
P4 <- ggplot(cell_cycle_result_melt,aes(x=factor(hour),y=value,fill="variable"))+
  geom_boxplot()+
  xlab("time (hour)")+ylab("S.Score")+
  scale_fill_manual(values=c("variable"="orange"))

cell_cycle_result_melt <-
  melt(cell_cycle_result_EL,id.vars=c("hour","exp"),
       measure.vars = c("PCNA","ACTB"))
P5 <- ggplot(cell_cycle_result_melt,aes(x=factor(hour),y=value,fill=variable))+
  geom_boxplot()+ylim(0,7)+
  xlab("time (hour)")+ylab("PCNA expression")+
  scale_fill_manual(values=c(variable="orange","gray"))


gridExtra::grid.arrange(P5,P4,P2,P1, nrow = 2)
#P1+P2+P3+P4


#ggplot(cell_cycle_result_EL,aes(x=S.Score,y=G2M.Score, color=Phase))+geom_point()




FeatureScatter(bulk,feature1="hour",feature2="TIPIN")

#house keeping gene expressions
P1 <- VlnPlot(bulk,features = "GAPDH",y.max=6)+ scale_y_continuous(limits = c(0,6))
P2 <- VlnPlot(bulk,features = "ACTB")+ scale_y_continuous(limits = c(0,6))
gridExtra::grid.arrange(P2,P1, nrow = 1)

bulk <- RunPCA(bulk, npcs=10, features = c(s_genes, g2m_genes))
#DimPlot(bulk)
DimHeatmap(bulk, dims = 1:2, cells = 10, balanced = TRUE)
DimPlot(bulk, group.by = "Phase", reduction = "pca")
DoHeatmap(bulk, features = c(s_genes,g2m_genes),slot="scale.data",group.by="Phase", label=TRUE)


bulk <- RunUMAP(bulk, dims=1:5, reduction.name = "UMAP")


celldataset <- as.cell_data_set(bulk)
celldataset <- cluster_cells(cds = celldataset, reduction_method = "UMAP")
celldataset <- learn_graph(celldataset, use_partition = TRUE)
celldataset <- order_cells(celldataset,reduction_method = "UMAP")
plot_cells(cds = celldataset,
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE)
bulk <- AddMetaData(object=bulk,
                    metadata = celldataset@principal_graph_aux@listData$UMAP$pseudotime,
                    col.name='pseudotime')


cellids= WhichCells (pbmc)

ggplot(cell_cycle_result,aes(x=hour,y=pseudotime,color=Phase))+geom_point()
ggplot(cell_cycle_result_EL,aes(x=hour,y=pseudotime,color=Phase))+geom_point()



#+xlim(c(0,25))

#par(new=TRUE)
#ggplot(cell_cycle_result_EL,aes(x=hour,color=exp,y=G2M.Score))+geom_point()

#CellScatter(object = pbmc, cell1 = cellids[9] , cell2 = cellids[10])
