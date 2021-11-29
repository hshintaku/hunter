#allcell <- pbmc
#allcell.markers <- pbmc.markers
DefaultAssay(allcell) <- "ADT"
allcell<-NormalizeData(allcell, normalization.method = "CLR", margin = 2)
DefaultAssay(allcell) <- "RNA"
FeaturePlot(allcell,features=c("mCherry","Azrite","Cyp27a1","Trp53","Pck1","Erbb2","Apoe","Alb"))
hepa_genes <-c("Cyp27a1","Pck1","Alb","Ppara","Mug1","Ces3a","Cyp3a25","Abcc2","Abcb11","Apoc4",
               "mt-Nd1","Sspn","Gm44805","Gm19951","Ftl1-ps1","Syt1","Tjp1","Wtap","Gm23240","Dmpk",
               "Clec4f","Cd5l","Ctsb","Tcf7l2","Scp2-ps2","Gm40841","Hnrnpr","Clpb","AC123870.1","Rein",
               "Prkg1","Ecm1","Rbms3","Raph1","Trp53","Pck1","Erbb2","Apoe","Twist2")
pheatmap(allcell[["RNA"]]@data[rownames(pbmc) %in% hepa_genes,],
         annotation_col = allcell[["gate"]],
         labels_col = NULL,
         cluster_cols = F)
#
#
# hepatocyte specific analysis
#
#hepa <- pbmc
#hepa.markers <- pbmc.markers
#hepa <-subset(pbmc,subset=plate=="p02",invert=TRUE)
#hepa <-subset(hepa,subset=plate=="P15",invert=TRUE)
#hepa <-subset(hepa,idents =c("1","3","4"))
#Idents(hepa) <- factor(x = Idents(hepa), levels = c(1,4,3))
#VlnPlot(pbmc,features=c("Actb","Gsn","Lmna","Ptk2b"))
#Actb/Gsn/Lmna/Ptk2b
#
# entrez ids of marker genes
#
cluster0 <- cluster_marker_entrez(hepa.markers,ms_ref,0,0.05,0.25)
cluster1 <- cluster_marker_entrez(hepa.markers,ms_ref,1,0.05,0.25)
cluster2 <- cluster_marker_entrez(hepa.markers,ms_ref,2,0.05,0.25)
cluster3 <- cluster_marker_entrez(hepa.markers,ms_ref,3,0.05,0.25)
cluster4 <- cluster_marker_entrez(hepa.markers,ms_ref,4,0.05,0.25)
# go analysis
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
Idents(hepa)<-factor(x=Idents(hepa),levels=c(3,2,0,4,1))
#
#
# visualize marker genes in a specific cluster
cluster1 <- cluster_marker_entrez(hepa.markers,ms_ref,1,0.001,0.9)
variable_genes <- hepa.markers[hepa.markers$p_val_adj<0.05 &
                                 hepa.markers$avg_log2FC>0.25 &
                                 hepa.markers$cluster==1,]
VlnPlot(hepa,features=variable_genes$gene)
# show go enrichment
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
hepa.data <- hepa[["RNA"]]@data
hepa.data.zone <- data.frame(t(hepa.data[variable_genes$gene,]))
cellids<- rownames(hepa.data.zone)
hepa.data.zone$plate <- substr(cellids,1,3)
#
# proximal score
dm <- DiffusionMap(hepa.data.zone)
hepa[["pseudospace"]]<-dm$DC1
#FeatureScatter(pbmc,feature1="pseudospace",feature2="Acot1")
variable_genes <- hepa.markers[hepa.markers$p_val_adj<0.0001 &
                                 hepa.markers$avg_log2FC>1.5 | hepa.markers$avg_log2FC< -1.5 &
                                 hepa.markers$cluster==1,]
pheatmap(hepa[["RNA"]]@data[variable_genes$gene,order(hepa[["pseudospace"]],decreasing=FALSE)],
         annotation_col = hepa[["pseudospace"]],
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

pheatmap(hepa[["RNA"]]@data[rownames(hepa) %in% hepa_genes,],
         annotation_col = allcell[["gate"]],
         cluster_rows = T,
         cluster_cols = F,
         labels_col = NULL)
#
# gse_result heatmap
upregulated <- gse_result@result[gse_result@result$ID=="GO:0048029",]
upregulated_entrez <- strsplit(upregulated$core_enrichment, split = "/")
upregulated_list <-ms_ref[ms_ref$entrez_annotation %in% upregulated_entrez[[1]],]
pheatmap(hepa[["RNA"]]@data[rownames(hepa) %in% upregulated_list$gene_short_name,
                            order(hepa[["pseudospace"]],decreasing=FALSE)],
         annotation_col = hepa[["pseudospace"]],
         cluster_cols = F,
         labels_col = NULL)
#
#
Idents(hepa) <- cluster$cluster
VlnPlot(hepa,features = g2m_genes, group.by = "cluster")
source(file.path(rdir,"hunter_Seurat_diffusionmap.R"))
hepa[["zone"]] <-hepa.data.zone$norm
FeatureScatter(hepa,feature1="zone",feature2="pseudospace",group.by = "plate")
FeaturePlot(hepa,features = c("zone","pseudospace"))

# correlated genes
# all variable gene heatmap
variable_genes <- hepa.markers[hepa.markers$p_val_adj<0.01 &
                                 hepa.markers$avg_log2FC>1.5,]
tree <-pheatmap(hepa[["RNA"]]@data[variable_genes$gene,order(hepa[["pseudospace"]],decreasing=FALSE)],
         annotation_col = hepa[["pseudospace"]],
         cluster_cols = F,
         labels_col = NULL,
         clustering_distance_rows="correlation")

tree$tree_row$order
cluster_attr <- cutree(tree$tree_row, k = 3)
mah_gene_cor <- rownames(data.frame(cluster_attr[cluster_attr==2]))
pheatmap(hepa[["RNA"]]@data[mah_gene_cor,order(hepa[["pseudospace"]],decreasing=FALSE)],
         annotation_col = hepa[["pseudospace"]],
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
View(ego_result@result)
