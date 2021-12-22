source(file.path(rdir,"hunter_entrez.R"))
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

cluster0 <- cluster_marker_entrez(pbmc.markers,ms_ref,0,0.001,1)
cluster1 <- cluster_marker_entrez(pbmc.markers,ms_ref,1,0.001,1)
cluster2 <- cluster_marker_entrez(pbmc.markers,ms_ref,2,0.01,1)
cluster3 <- cluster_marker_entrez(pbmc.markers,ms_ref,3,0.001,1)
cluster4 <- cluster_marker_entrez(pbmc.markers,ms_ref,4,0.001,1)
cluster5 <- cluster_marker_entrez(pbmc.markers,ms_ref,5,0.001,1)
#https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/compareCluster
perturbed_gene_cancer <- list(Cluster0=cluster0$entrez_annotation,
                            Cluster2=cluster2$entrez_annotation)
perturbed_gene_hepa <- list(Cluster1=cluster1$entrez_annotation,
                              Cluster3=cluster3$entrez_annotation,
                              Cluster4=cluster4$entrez_annotation,
                              Cluster5=cluster5$entrez_annotation)
cancer_go <- compareCluster(perturbed_gene_cancer, fun="enrichGO",
                     OrgDb         = org.Mm.eg.db)
hepa_go <- compareCluster(perturbed_gene_hepa, fun="enrichGO",
                     OrgDb         = org.Mm.eg.db)
#summary(xx)

p1 <- DimPlot(pbmc, reduction = "umap",group.by = "plate")
p2<-FeaturePlot(pbmc,features="Saa1")
p3<-DimPlot(pbmc)
p4<-ggplot(cluster_density,aes(y=cluster,x=plate$plate,fill=count))+ 
  geom_tile()+ theme_classic()+
  scale_fill_gradientn(colours = c("white", "red"))+geom_text(aes(label = count)) 

p5<- dotplot(hepa_go)
p6 <- dotplot(cancer_go)
gridExtra::grid.arrange(p5,p6,nrow=1)
#clusterProfiler::dotplot(ego_result)
#clusterProfiler::emapplot(ego_result.simple)
#clusterProfiler::cnetplot(ego_result, categoryS ize="pvalue")
write.xlsx(summary(xx),"/home/samba/public/shintaku/github/hunter2/minegishi/go_compare_summary.xlsx")
