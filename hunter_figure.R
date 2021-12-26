library(destiny)
library(ggplot2)
library(conflicted)
suppressPackageStartupMessages(library(scran))
library(purrr)
library(wordspace)
library(pheatmap)
conflict_prefer("mutate", "dplyr")
hepa1 <- subset(hepa1,subset = gate ==c("g1"),invert=TRUE)
hepa1 <- subset(hepa1,subset = gate ==c("g4"),invert=TRUE)
hepa1 <- subset(hepa1,subset = plate ==c("p04"),invert=TRUE)
hepa1 <- subset(hepa1,subset = plate ==c("p05"),invert=TRUE)
hepa1[["percent.mt"]] <- PercentageFeatureSet(hepa1, pattern = "^mt-")
hepa1 <- subset(hepa1, subset= percent.mt<5)
hepa1 <- NormalizeData(hepa1, normalization.method = "LogNormalize", scale.factor = 1e5)
hepa1 <- FindVariableFeatures(hepa1, selection.method = "vst", nfeatures = 1000)

hepa2 <- subset(hepa2,subset = gate ==c("g1"),invert=TRUE)
hepa2 <- subset(hepa2,subset = gate ==c("g4"),invert=TRUE)
hepa2 <- subset(hepa2,subset = plate ==c("p04"),invert=TRUE)
hepa2 <- subset(hepa2,subset = plate ==c("p05"),invert=TRUE)
hepa2[["percent.mt"]] <- PercentageFeatureSet(hepa2, pattern = "^mt-")
hepa2 <- subset(hepa2, subset= percent.mt<5)
hepa2 <- NormalizeData(hepa2, normalization.method = "LogNormalize", scale.factor = 1e5)
hepa2 <- FindVariableFeatures(hepa2, selection.method = "vst", nfeatures = 1000)

hepa3[["percent.mt"]] <- PercentageFeatureSet(hepa3, pattern = "^mt-")
hepa3 <- subset(hepa3, subset= percent.mt<5)
hepa3 <- NormalizeData(hepa3, normalization.method = "LogNormalize", scale.factor = 1e5)
hepa3 <- FindVariableFeatures(hepa3, selection.method = "vst", nfeatures = 1000)
#
# integrate data with batch correction
# 
batch2<- merge(hepa2, y = hepa3,  project = "hunter")
integrated <- merge(hepa1,batch2,project="hunter")
#hepa.list <-list(hepa1,allcell)
#anchors <- FindIntegrationAnchors(object.list = hepa.list)
#integrated <- IntegrateData(anchorset = anchors)

pbmc <- integrated
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000)

source(file.path(rdir,"hunter_Seurat_biomart.R"))
source(file.path(rdir,'hunter_Seurat_clustering.R'))
allcell <- pbmc
allcell.markers <- pbmc.markers
DimPlot(allcell)+DimPlot(allcell,group.by = "plate")
DefaultAssay(allcell) <- "ADT"
allcell<-NormalizeData(allcell, normalization.method = "CLR", margin = 2)
DefaultAssay(allcell) <- "RNA"
FeaturePlot(allcell,features=c("mCherry","Azrite","Cyp27a1","Trp53","Pck1","Erbb2","Apoe","Alb"))
hepa_genes <-c("Cyp27a1","Pck1","Alb","Ppara","Mug1","Ces3a","Cyp3a25","Abcc2","Abcb11","Apoc4",
               "mt-Nd1","Sspn","Gm44805","Gm19951","Ftl1-ps1","Syt1","Tjp1","Wtap","Gm23240","Dmpk",
               "Clec4f","Cd5l","Ctsb","Tcf7l2","Scp2-ps2","Gm40841","Hnrnpr","Clpb","AC123870.1","Rein",
               "Prkg1","Ecm1","Rbms3","Raph1","Trp53","Pck1","Erbb2","Apoe","Twist2")
pheatmap(allcell[["RNA"]]@data[rownames(allcell) %in% hepa_genes,],
         annotation_col = allcell[["gate"]],
         labels_col = NULL,
         cluster_cols = F)
source(file.path(rdir,"hunter_SingleCellSignal.R"))
source(file.path(rdir,"hunter_clusterProfiler_GO.R"))
#
#
# hepatocyte specific analysis
#
hepa <- subset(pbmc,subset = plate ==c("p02"),invert=TRUE)
hepa <- subset(hepa,subset = plate ==c("P15"),invert=TRUE)
hepa <- subset(hepa,subset = plate ==c("p20"),invert=TRUE)
pbmc <- hepa
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 1e5)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000)

cellids <- colnames(pbmc)
source(file.path(rdir,'hunter_Seurat_clustering.R'))

hepa <- pbmc
hepa.markers <- pbmc.markers
DefaultAssay(hepa) <- "ADT"
hepa<-NormalizeData(hepa, normalization.method = "CLR", margin = 2)
DefaultAssay(hepa) <- "RNA"
DimPlot(hepa)+DimPlot(hepa,group.by = "plate")
#Idents(hepa)<-factor(x=Idents(hepa),levels=c(2,1,0))
#
# entrez ids of marker genes
#
cluster0 <- cluster_marker_entrez(hepa.markers,ms_ref,0,0.01,0.3)
cluster1 <- cluster_marker_entrez(hepa.markers,ms_ref,1,0.01,0.3)
cluster2 <- cluster_marker_entrez(hepa.markers,ms_ref,2,0.01,0.3)
#cluster3 <- cluster_marker_entrez(hepa.markers,ms_ref,3,0.01,0.3)
#cluster4 <- cluster_marker_entrez(hepa.markers,ms_ref,4,0.05,0.25)
# go analysis
perturbed_gene_hepa <- list(Cluster1=cluster1$entrez_annotation,
                            Cluster0=cluster0$entrez_annotation)
                            #Cluster3=cluster3$entrez_annotation,
                            #Cluster0=cluster0$entrez_annotation)
hepa_go <- compareCluster(perturbed_gene_hepa, fun="enrichGO",
                          OrgDb         = org.Mm.eg.db)
dotplot(hepa_go)
# MAH vs distal 
perturbed_gene_hepa <- list(Cluster1=cluster1$entrez_annotation,
                            Cluster0=cluster0$entrez_annotation)
hepa_go <- compareCluster(perturbed_gene_hepa, fun="enrichGO",
                          OrgDb         = org.Mm.eg.db)
dotplot(hepa_go)
hepa_go@compareClusterResult[hepa_go@compareClusterResult$Cluster=="Cluster0",]$Description
hepa_go@compareClusterResult[hepa_go@compareClusterResult$Cluster=="Cluster0",]$geneID
#
#
#
# visualize marker genes in a specific cluster
#cluster1 <- cluster_marker_entrez(hepa.markers,ms_ref,1,0.001,0.9)
variable_genes <- hepa.markers[hepa.markers$p_val_adj<0.001 &
                                 hepa.markers$avg_log2FC>1.5 &
                                 hepa.markers$cluster==0,]
# show go enrichment
ego_result <- enrichGO(gene          = ms_ref[ms_ref$gene_short_name %in% variable_genes$gene,]$entrez_annotation, 
                       OrgDb         = org.Mm.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
barplot(ego_result, drop=TRUE, showCategory=30)
View(ego_result@result)
#VlnPlot(hepa,features=variable_genes$gene)
#
# MAH gene diffusion map
# 
hepa.data <- hepa[["RNA"]]@data
hepa.data.zone <- data.frame(t(hepa.data[c(variable_genes$gene),]))
cellids<- rownames(hepa.data.zone)
hepa.data.zone$plate <- substr(cellids,1,3)
#
# proximal score
dm <- DiffusionMap(hepa.data.zone)
hepa[["pseudospace"]]<-dm$DC1
pheatmap(hepa[["RNA"]]@data[variable_genes$gene,order(hepa[["pseudospace"]],decreasing=FALSE)],
         annotation_col = hepa[["pseudospace"]],
         cluster_cols = F,
         labels_col = NULL)
FeatureScatter(hepa,feature1="Venus",feature2="pseudospace",group.by = "plate")
#FeatureScatter(pbmc,feature1="pseudospace",feature2="Acot1")
#
# visualize result with all other marker genes
hepa_genes <- hepa.markers[hepa.markers$p_val_adj<0.0001 &
                                 hepa.markers$avg_log2FC>0.8,]
pheatmap(hepa[["RNA"]]@data[hepa_genes$gene,order(hepa[["pseudospace"]],decreasing=FALSE)],
         annotation_col = hepa[["pseudospace"]],
         cluster_cols = F,
         labels_col = NULL)
#
# gse_result heatmap
#
# pathway analysis
source(file.path(ridir,"hunter_clusterProfiler_GSEA.R"))
#
#
# zonation vs pseudospace
#Idents(hepa) <- cluster$cluster
#VlnPlot(hepa,features = g2m_genes, group.by = "cluster")
source(file.path(rdir,"hunter_Seurat_diffusionmap.R"))
pbmc[["zone"]] <-hepa.data.zone$norm
FeatureScatter(pbmc,feature1="zone",feature2="pseudospace",group.by = "plate")
FeaturePlot(hepa,features = c("zone","pseudospace"))

# correlated genes
# all variable gene heatmap
#
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
ego_result_spread <- ego_result@result
