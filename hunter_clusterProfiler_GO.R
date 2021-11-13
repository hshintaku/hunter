
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

cluster.marker <- pbmc.markers[pbmc.markers$cluster==2 & pbmc.markers$p_val_adj<0.01 & pbmc.markers$avg_log2FC< -2,]
perturbed_gene <- ms_ref[ms_ref$gene_short_name %in% cluster.marker$gene,]
perturbed_gene <- perturbed_gene[!is.na(perturbed_gene$entrez_annotation),]
#perturbed_gene_RG.entrez_annotation <- ms_ref_subset[ms_ref_subset$gene_short_name %in% rownames(perturbed_gene_RG),]
#perturbed_gene_elp.entrez_annotation <- ms_ref_subset[ms_ref_subset$gene_short_name %in% rownames(perturbed_gene_elp),]

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

cluster.marker.n <- pbmc.markers[pbmc.markers$cluster==2 & pbmc.markers$p_val_adj<0.01 & pbmc.markers$avg_log2FC< -1,]
perturbed_gene.n <- ms_ref[ms_ref$gene_short_name %in% cluster.marker.n$gene,]
perturbed_gene.n <- perturbed_gene.n[!is.na(perturbed_gene.n$entrez_annotation),]

#https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/compareCluster
perturbed_gene_list <- list(positive=perturbed_gene$entrez_annotation,negative=perturbed_gene.n$entrez_annotation)
xx <- compareCluster(perturbed_gene_list, fun="groupGO",
                     OrgDb         = org.Mm.eg.db)
summary(xx)
#clusterProfiler::dotplot(ego_result)
#clusterProfiler::emapplot(ego_result.simple)
#clusterProfiler::cnetplot(ego_result, categoryS ize="pvalue")

