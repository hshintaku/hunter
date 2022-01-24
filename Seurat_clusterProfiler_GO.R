#
# GOenrichmentAnalysis (experimental)
#BiocManager::install("org.Mm.eg.db")
#library("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library("org.Rn.eg.db")
library(clusterProfiler)
#
ensg <- data.frame(unlist(as.list(org.Hs.egENSEMBL2EG)))
ensrg <- data.frame(unlist(as.list(org.Rn.egENSEMBL2EG)))
colnames(ensg)<-"gene"
colnames(ensrg)<-"gene"
# entrez annotation
#ms_ref$entrez_annotation <- ensg[ms_ref$ensembl_gene_id,]
#ms_ref[ms_ref$ensembl_gene_id,]$entrez_annotation <- ensrg[ms_ref$ensembl_gene_id,]
en_ref <- rbind(ensg,ensrg)
# add entrez_annotation
ms_ref$entrez_annotation<-en_ref[ms_ref$ensembl_gene_id,]
# remove genes not found in entrez
ms_ref<-ms_ref[!is.na(ms_ref$entrez_annotation),]
# de_genes_enteez
de_genes_entrez <- ms_ref[ms_ref$gene_short_name %in% rownames(de_genes),]$entrez_annotation
#
# https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
#
# clusterProfiler
# 
#gene_module_go <- allLLIDs[moduleColors_subset=="pink"]

ego_result <- enrichGO(gene          = de_genes_entrez, 
                       OrgDb         = org.Rn.eg.db, # rattus, human etc
                       ont           = "CC",
                         pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
head(as.data.frame(ego_result))
ego_result.simple<-simplify(ego_result)
head(as.data.frame(ego_result.simple))
barplot(ego_result, drop=TRUE, showCategory=30)
#https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/compareCluster
perturbed_gene_list <- list(RG=perturbed_gene_RG.entrez_annotation$entrez_annotation,HEA=perturbed_gene_HEA.entrez_annotation$entrez_annotation)
xx <- compareCluster(perturbed_gene_list, fun="groupGO",
                     OrgDb         = org.Hs.eg.db)
summary(xx)
#clusterProfiler::dotplot(ego_result)
#clusterProfiler::emapplot(ego_result.simple)
#clusterProfiler::cnetplot(ego_result, categoryS ize="pvalue")
goplot(ego_result.simple)

