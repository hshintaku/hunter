cluster_marker_entrez <- function(pbmc.markers,ms_ref,cluster.num,p.thres,log2FC.thres){
  cluster.marker.n <- pbmc.markers[pbmc.markers$cluster==cluster.num &
                                     pbmc.markers$p_val_adj<p.thres &
                                     pbmc.markers$avg_log2FC > log2FC.thres,]
  perturbed_gene.n <- ms_ref[ms_ref$gene_short_name %in% cluster.marker.n$gene,]
  perturbed_gene.n <- perturbed_gene.n[!is.na(perturbed_gene.n$entrez_annotation),]
  return(perturbed_gene.n)
}
#
# GOenrichmentAnalysis (experimental)
#BiocManager::install("org.Mm.eg.db")
#library("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
#install_github("https://github.com/YuLab-SMU/clusterProfiler")
library("org.Hs.eg.db")
library(clusterProfiler)
#
ensmusg <- data.frame(unlist(as.list(org.Hs.egENSEMBL2EG)))
# entrez annotation
ms_ref$entrez_annotation <- ensmusg[ms_ref$ensembl_gene_id,]

cluster0 <- cluster_marker_entrez(pbmc.markers,ms_ref,0,0.001,1)
cluster1 <- cluster_marker_entrez(pbmc.markers,ms_ref,1,0.001,1)
cluster2 <- cluster_marker_entrez(pbmc.markers,ms_ref,2,0.01,1)
cluster3 <- cluster_marker_entrez(pbmc.markers,ms_ref,3,0.001,1)
cluster4 <- cluster_marker_entrez(pbmc.markers,ms_ref,4,0.001,1)

cluster5 <- cluster_marker_entrez(pbmc.markers,ms_ref,5,0.05,0)
#
# https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
#
# clusterProfiler
# 
#gene_module_go <- allLLIDs[moduleColors_subset=="pink"]
#
# ont = CC, BP, MF
ego_result <- enrichGO(gene          = cluster5$entrez_annotation, 
                       OrgDb         = org.Hs.eg.db,
                       ont           = "CC",
                         pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)

barplot(ego_result, drop=TRUE, showCategory=30)
head(as.data.frame(ego_result))
ego_result.simple<-simplify(ego_result)
head(as.data.frame(ego_result.simple))
#https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/compareCluster
perturbed_gene <- list(Cluster0=cluster0$entrez_annotation,
                            Cluster1=cluster1$entrez_annotation,
                            Cluster2=cluster2$entrez_annotation,
                            Cluster3=cluster3$entrez_annotation,
                            Cluster4=cluster4$entrez_annotation,
                            Cluster5=cluster5$entrez_annotation)
goterm <- compareCluster(perturbed_gene, fun="enrichGO",
                            OrgDb         = org.Hs.eg.db)
#p1 <- DimPlot(pbmc, reduction = "umap",group.by = "plate")
#p2<-FeaturePlot(pbmc,features="Saa1")
#p3<-DimPlot(pbmc)
#p4<-ggplot(cluster_density,aes(y=cluster,x=plate$plate,fill=count))+ 
#  geom_tile()+ theme_classic()+
#  scale_fill_gradientn(colours = c("white", "red"))+geom_text(aes(label = count)) 

dotplot(goterm)
#gridExtra::grid.arrange(p5,p6,nrow=1)
#goplot(ego_result.simple)

p1 <- DimPlot(pbmc, reduction = "pca")
p2 <- DimPlot(pbmc, reduction = "umap")
p3<-FeaturePlot(pbmc,features="FLDtotal")
p1+p2+p3
