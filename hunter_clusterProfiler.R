#
# GOenrichmentAnalysis (experimental)
library("org.Mm.eg.db")
library(clusterProfiler)
#
ensmusg <- data.frame(unlist(as.list(org.Mm.egENSEMBL2EG)))
# entrez annotation
ms_ref$entrez_annotation <- ensmusg[ms_ref$ensembl_gene_id,]
#
#
datExpr_names <- data.frame(names(datExpr)) # list gene short names
datExpr_ref <- ms_ref[match(datExpr_names$names.datExpr.,ms_ref$gene_short_name),] # create reference from ms_ref regardless of entrez annotation
datExp_annot_index <- data.frame(which(!is.na(datExpr_ref$entrez_annotation))) #

ms_ref_subset <- datExpr_ref[which(!is.na(datExpr_ref$entrez_annotation)),]#ms_ref[which(!is.na(ms_ref$entrez_annotation)),]

datExpr_subset <- datExpr[,unlist(datExp_annot_index)]
moduleColors_subset <- moduleColors[unlist(datExp_annot_index)]
allLLIDs <- ms_ref_subset$entrez_annotation

GOenr <- GOenrichmentAnalysis(moduleColors_subset,allLLIDs, organism = "mouse", ontologies = c("BP", "CC", "MF"), nBestP = 10)
tab <- GOenr$bestPTerms[[4]]$enrichment

module_genes <- names(datExpr_subset)[moduleColors_subset=="pink"]
write.csv(module_genes,file.path(wdir,"module_genes_pink.csv"))
#
# https://bioc.ism.ac.jp/packages/3.3/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
#
# clusterProfiler
# 
gene_module_go <- allLLIDs[moduleColors_subset=="pink"]
ego_result <- enrichGO(gene          = gene_module_go,
                       OrgDb         = org.Mm.eg.db,
                       ont           = "CC",
                         pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05, 
                       readable      = TRUE)
head(as.data.frame(ego_result))
ego_result.simple<-simplify(ego_result)
head(as.data.frame(ego_result.simple))
barplot(ego_result, drop=TRUE, showCategory=30)
clusterProfiler::dotplot(ego_result)
#clusterProfiler::emapplot(ego_result.simple)
clusterProfiler::cnetplot(ego_result, categorySize="pvalue", foldChange=allLLIDs)
goplot(ego_result.simple)
#
# gse
# 
gse_result<- gseGO(geneList     = gene_module_go,
                   OrgDb        = org.Mm.eg.db,
                   ont          = "BP",
                     nPerm        = 1000,
                   minGSSize    = 120,
                   pvalueCutoff = 0.05,
                   verbose      = FALSE)
head(as.data.frame(gse_result))
ridgeplot(gse_result,showCategory = 8)
