library(monocle)
library(R.matlab)
library(tools)

hepa <- subset(pbmc,subset=plate=="p02",invert=TRUE)
hepa <- subset(hepa,subset=plate=="P15",invert=TRUE)

pbmc_monocle <- as.CellDataSet(liver)
pbmc_monocle <- estimateSizeFactors(pbmc_monocle)
pbmc_monocle <- estimateDispersions(pbmc_monocle)
disp_table <- dispersionTable(pbmc_monocle)
rownames(disp_table)<-disp_table$gene_id

genes_zonation <- readMat("/home/samba/public/shintaku/matlab_code/Zonation_params.mat")
genes_cv <- data.frame(unlist(genes_zonation$genes.cv))
colnames(genes_cv)<-"gene_id"
rownames(genes_cv)<- genes_cv$gene_id
genes_pn <- data.frame(unlist(genes_zonation$genes.pn))
colnames(genes_pn)<-"gene_id"
rownames(genes_pn)<-genes_pn$gene_id

ordering_genes <- rbind(genes_pn,genes_cv)
ordering_genes$gene_id <- toTitleCase(ordering_genes$gene_id)
ordering_genes_disp <- disp_table[disp_table$gene_id %in% ordering_genes$gene_id,]
ordering_genes <- ordering_genes[ordering_genes$gene_id %in% disp_table$gene_id,]

pbmc_monocle <- setOrderingFilter(pbmc_monocle, ordering_genes_disp)
pbmc_monocle <- reduceDimension(pbmc_monocle)

pbmc_monocle <- orderCells(pbmc_monocle)
plot_cell_trajectory(pbmc_monocle, markers = c("Cyp2e1","Cyp2f2"), use_color_gradient = TRUE)
#p2<-plot_cell_trajectory(pbmc_monocle, color_by = )
#p1+p2

plot_genes_in_pseudotime(pbmc_monocle, color_by = "cell")


Pseudotime <- data.frame(pbmc_monocle$Pseudotime)
rownames(Pseudotime) <- colnames(pbmc_monocle)



hepa[["Pseudotime"]] <- Pseudotime

FeatureScatter(hepa,feature1 = "Pseudotime",feature2 = "Arg1",group.by = "plate")
#hepa$Pseudotime
#colnames(Pseudotime)<-"Pseudotime"
#cellnames <- data.frame(rownames(Pseudotime))
#cellnames <- cellnames[order(Pseudotime),]
#p<-DoHeatmap(hepa,features = ordering_genes$gene_id)
#p$data$cell <- factor(p$data$Cell, levels=cellnames)
#p



hepa.data <- pbmc[["RNA"]]@data
hepa.data.zone <- hepa.data[ordering_genes[1:30],]
library(pheatmap)
pheatmap(hepa.data.zone,
         cluster_rows = FALSE,cluster_cols = TRUE,
         scale = "row")
