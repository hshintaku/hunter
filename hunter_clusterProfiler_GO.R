cluster_marer_entrez <- function(pbmc.markers,ms_ref,cluster.num,p.thres,log2FC.thres){
  cluster.marker.n <- pbmc.markers[pbmc.markers$cluster==cluster.num &
                                     pbmc.markers$p_val_adj<p.thres &
                                     pbmc.markers$avg_log2FC>log2FC.thres,]
  perturbed_gene.n <- ms_ref[ms_ref$gene_short_name %in% cluster.marker.n$gene,]
  perturbed_gene.n <- perturbed_gene.n[!is.na(perturbed_gene.n$entrez_annotation),]
  return(perturbed_gene.n)
}
#
# GOenrichmentAnalysis (experimental)
BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
#library("org.Hs.eg.db")
library(clusterProfiler)
#
#ensmusg <- data.frame(unlist(as.list(org.Hs.egENSEMBL2EG)))
ensmusg <- data.frame(unlist(as.list(org.Mm.egENSEMBL2EG)))
# entrez annotation
ms_ref$entrez_annotation <- ensmusg[ms_ref$ensembl_gene_id,]

cluster.marker <- cluster_marer_entrez(pbmc.markers,ms_ref,2,0.01,1)
#
# https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
#
# clusterProfiler
# 
#gene_module_go <- allLLIDs[moduleColors_subset=="pink"]

ego_result <- enrichGO(gene          = perturbed_gene$entrez_annotation, 
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                         pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE, showCategory=30)
#head(as.data.frame(ego_result))
ego_result.simple<-simplify(ego_result)
goplot(ego_result.simple)
#head(as.data.frame(ego_result.simple))

cluster0 <- cluster_marer_entrez(pbmc.markers,ms_ref,0,0.01,2)
cluster1 <- cluster_marer_entrez(pbmc.markers,ms_ref,1,0.01,2)
cluster2 <- cluster_marer_entrez(pbmc.markers,ms_ref,2,0.01,1)
cluster3 <- cluster_marer_entrez(pbmc.markers,ms_ref,3,0.01,2)
cluster4 <- cluster_marer_entrez(pbmc.markers,ms_ref,4,0.01,2)
cluster5 <- cluster_marer_entrez(pbmc.markers,ms_ref,5,0.01,2)
#https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/compareCluster
perturbed_gene_list <- list(Cluster0=cluster0$entrez_annotation,
                            Cluster1=cluster1$entrez_annotation,
                            Cluster2=cluster2$entrez_annotation,
                            Cluster3=cluster3$entrez_annotation,
                            Cluster4=cluster4$entrez_annotation,
                            Cluster5=cluster5$entrez_annotation)
xx <- compareCluster(perturbed_gene_list, fun="enrichGO",
                     OrgDb         = org.Mm.eg.db)
summary(xx)

p1 <- DimPlot(pbmc, reduction = "umap",group.by = "plate")
p2<-FeaturePlot(pbmc,features="Saa1")
p3<-DimPlot(pbmc)
p4<-ggplot(cluster_density,aes(y=cluster,x=plate$plate,fill=count))+ 
  geom_tile()+ theme_classic()+
  scale_fill_gradientn(colours = c("white", "red"))+geom_text(aes(label = count)) 

p5<- dotplot(xx)
p1+p3+p4
#clusterProfiler::dotplot(ego_result)
#clusterProfiler::emapplot(ego_result.simple)
#clusterProfiler::cnetplot(ego_result, categoryS ize="pvalue")

