DefaultAssay(pbmc) <- "ADT"
pbmc<-NormalizeData(pbmc, normalization.method = "CLR", margin = 2)
DefaultAssay(pbmc) <- "RNA"
FeaturePlot(pbmc,features=c("mCherry","Azrite","Cyp27a1","Trp53","Pck1","Erbb2","Apoe","Alb"))

hepa <-subset(pbmc,subset=plate=="p02",invert=TRUE)
hepa <-subset(hepa,subset=plate=="P15",invert=TRUE)
hepa <-subset(hepa,idents =c("1","3","4"))

Idents(hepa) <- factor(x = Idents(hepa), levels = c(1,4,3))
VlnPlot(pbmc,features=c("Actb","Gsn","Lmna","Ptk2b"))
#Actb/Gsn/Lmna/Ptk2b

cluster0 <- cluster_marker_entrez(pbmc.markers,ms_ref,0,0.05,0.25)
cluster1 <- cluster_marker_entrez(pbmc.markers,ms_ref,1,0.05,0.25)
cluster2 <- cluster_marker_entrez(pbmc.markers,ms_ref,2,0.05,0.25)
cluster3 <- cluster_marker_entrez(pbmc.markers,ms_ref,3,0.05,0.25)
cluster4 <- cluster_marker_entrez(pbmc.markers,ms_ref,4,0.05,0.25)

perturbed_gene_hepa <- list(Cluster3=cluster3$entrez_annotation,
                            Cluster2=cluster2$entrez_annotation,
                            Cluster0=cluster0$entrez_annotation,
                            Cluster4=cluster4$entrez_annotation,
                            Cluster1=cluster1$entrez_annotation)
hepa_go <- compareCluster(perturbed_gene_hepa, fun="enrichGO",
                          OrgDb         = org.Mm.eg.db)
dotplot(hepa_go)

# MAH vs distal 
perturbed_gene_hepa <- list(Cluster4=cluster4$entrez_annotation,
                            Cluster1=cluster1$entrez_annotation)
hepa_go <- compareCluster(perturbed_gene_hepa, fun="enrichGO",
                          OrgDb         = org.Mm.eg.db)
dotplot(hepa_go)
#
Idents(pbmc)<-factor(x=Idents(pbmc),levels=c(3,2,0,4,1))
cluster1 <- cluster_marker_entrez(pbmc.markers,ms_ref,1,0.001,0.9)
variable_genes <- pbmc.markers[pbmc.markers$p_val_adj<0.001 &
                                 pbmc.markers$avg_log2FC>0.9 &
                                 pbmc.markers$cluster==1,]
VlnPlot(pbmc,features=variable_genes$gene)

variable_genes<- pbmc.markers[pbmc.markers$p_val_adj<0.05 & pbmc.markers$avg_log2FC>0.25
                              & pbmc.markers$cluster==1,]

ego_result <- enrichGO(gene          = ms_ref[ms_ref$gene_short_name %in% variable_genes$gene,]$entrez_annotation, 
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE, showCategory=30)
ego_result_spread <- ego_result@result

#
# MAH gene diffusion map
# 
hepa.data <- pbmc[["RNA"]]@data
hepa.data.zone <- data.frame(t(hepa.data[variable_genes$gene,]))
cellids<- rownames(hepa.data.zone)
hepa.data.zone$plate <- substr(cellids,1,3)

dm <- DiffusionMap(hepa.data.zone)
pbmc[["pseudospace"]]<-dm$DC1
#FeatureScatter(pbmc,feature1="pseudospace",feature2="Acot1")
pheatmap(pbmc[["RNA"]]@data[variable_genes$gene,order(pbmc[["pseudospace"]],decreasing=FALSE)],
         annotation_col = pbmc[["pseudospace"]],
         cluster_cols = F,
         labels_col = NULL)
VlnPlot(pbmc,features="pseudospace")
VlnPlot(pbmc,features=variable_genes$gene)
VlnPlot(pbmc,features=hepa_genes)
FeatureScatter(pbmc,feature1="Venus",feature2="pseudospace",group.by = "plate")
#
#
# hepatocyte marker/transcription factor heatmap
hepa_genes<- pbmc.markers[pbmc.markers$p_val_adj<0.01 & pbmc.markers$avg_log2FC>1.5,]
hepa_genes <-c("Cyp27a1","Pck1","Alb","Ppara","Mug1","Ces3a","Cyp3a25","Abcc2","Abcb11","Apoc4",
               "mt-Nd1","Sspn","Gm44805","Gm19951","Ftl1-ps1","Syt1","Tjp1","Wtap","Gm23240","Dmpk",
               "Clec4f","Cd5l","Ctsb","Tcf7l2","Scp2-ps2","Gm40841","Hnrnpr","Clpb","AC123870.1","Rein",
               "Prkg1","Ecm1","Rbms3","Raph1")
pheatmap(pbmc[["RNA"]]@data[rownames(pbmc) %in% hepa_genes,order(pbmc[["pseudospace"]],decreasing=FALSE)],
         annotation_col = pbmc[["pseudospace"]],
         cluster_cols = F,
         labels_col = NULL)

#
# gse_result heatmap
upregulated <- gse_result@result[gse_result@result$ID=="GO:0048029",]
upregulated_entrez <- strsplit(upregulated$core_enrichment, split = "/")
upregulated_list <-ms_ref[ms_ref$entrez_annotation %in% upregulated_entrez[[1]],]
pheatmap(pbmc[["RNA"]]@data[rownames(pbmc) %in% upregulated_list$gene_short_name,order(pbmc[["pseudospace"]],decreasing=FALSE)],
         annotation_col = pbmc[["pseudospace"]],
         cluster_cols = F,
         labels_col = NULL)
#
#
Idents(pbmc) <- cluster$cluster
VlnPlot(pbmc,features = g2m_genes, group.by = "cluster")

pbmc[["zone"]] <-hepa.data.zone$norm
FeatureScatter(pbmc,feature1="zone",feature2="pseudospace",group.by = "plate")
FeaturePlot(pbmc,features = c("zone","pseudospace"))

# correlated genes
# all variable gene heatmap
variable_genes <- pbmc.markers[pbmc.markers$p_val_adj<0.01 &
                                 pbmc.markers$avg_log2FC>1.5,]
tree <-pheatmap(pbmc[["RNA"]]@data[variable_genes$gene,order(pbmc[["pseudospace"]],decreasing=FALSE)],
         annotation_col = pbmc[["pseudospace"]],
         cluster_cols = F,
         labels_col = NULL,
         clustering_distance_rows="correlation")

tree$tree_row$order
cluster_attr <- cutree(tree$tree_row, k = 3)
mah_gene_cor <- rownames(data.frame(cluster_attr[cluster_attr==2]))
pheatmap(pbmc[["RNA"]]@data[mah_gene_cor,order(pbmc[["pseudospace"]],decreasing=FALSE)],
         annotation_col = pbmc[["pseudospace"]],
         cluster_cols = F,
         labels_col = NULL,
         clustering_distance_rows="correlation")
ego_result <- enrichGO(gene          = ms_ref[ms_ref$gene_short_name %in% mah_gene_cor,]$entrez_annotation, 
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE, showCategory=20)
