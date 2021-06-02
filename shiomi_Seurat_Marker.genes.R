#
# 5Aza cells/RG cells
#
RG.markers <- FindMarkers(pbmc,logfc.threshold = 0.1,ident.1=colnames(subset(pbmc,subset=gate=="RG")))
perturbed_gene_RG <- subset(RG.markers, p_val_adj<0.2)
Aza.markers <- FindMarkers(pbmc,logfc.threshold = 0.1,ident.1=colnames(subset(pbmc,subset=cell=="HEA")))
perturbed_gene_HEA <- subset(Aza.markers, p_val_adj<0.2)

perturbed_gene <- list(RG=rownames(perturbed_gene_RG),HEA=rownames(perturbed_gene_HEA))
p0 <- venn.diagram(perturbed_gene,filename =NULL, fill=c(2,3), alpha=0.4, lty=3)
grid::grid.newpage()
grid::grid.draw(p0)
#
# intersection of 5Aza and RG
#
shared_genes_RG_HEA<- intersect(rownames(perturbed_gene_HEA),rownames(perturbed_gene_RG))
shared_genes_RG_HEA <- ms_ref[ms_ref$gene_short_name %in% shared_genes_RG_HEA,]
perturbed_expression <- pbmc[["RNA"]]@data[shared_genes_RG_HEA$gene_short_name,]
pheatmap(perturbed_expression,
         annotation_col = pbmc[[c("gate","cell")]],
         cluster_cols = TRUE,cluster_rows = TRUE,scale='row')

#
# control cells
#
cntl_HeLa <- subset(x=pbmc,subset=dish!=c("HEA") )
elp.markers <- FindMarkers(cntl_HeLa,logfc.threshold = 0.1,ident.1=colnames(subset(cntl_HeLa,subset=gate=="RG")))
perturbed_gene_elp <- subset(elp.markers, p_val_adj<0.05)
perturbed_gene_elp
write.csv(RG_marker_regardless_of_5Aza,"/home/samba/pihome/2021/Shintaku/Shiomi/RG_marker_regardless_of_5Aza_out_of_15genes.csv")


#
# intersection of 5Aza and RG HeLa
#
RG_marker_regardless_of_5Aza <-intersect(shared_genes_RG_HEA$gene_short_name,rownames(perturbed_gene_elp))
RG_marker_regardless_of_5Aza <-ms_ref[ms_ref$gene_short_name %in% RG_marker_regardless_of_5Aza,]

perturbed_expression <- pbmc[["RNA"]]@data[RG_marker_regardless_of_5Aza$gene_short_name,]
pheatmap(perturbed_expression,
         annotation_col = pbmc[[c("gate","cell")]],
         cluster_cols = TRUE,cluster_rows = TRUE,scale='row')

FeatureScatter(subset(x=cntl_HeLa,subset=cell=="HEE"),feature1="S.Score",feature2 = "G2M.Score",group.by = "gate")

#
#
#
write.csv(perturbed_gene_HEA,file="/home/samba/pihome/2021/Shintaku/Shiomi/marker.genes_HEA.csv")
write.csv(perturbed_gene_RG,file="/home/samba/pihome/2021/Shintaku/Shiomi/marker.genes_RG.csv")
write.csv(shared_genes_RG_HEA,file="/home/samba/pihome/2021/Shintaku/Shiomi/shared_marker.genes_RG_HEA.csv")
