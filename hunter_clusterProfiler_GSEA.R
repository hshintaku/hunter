#
# https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/scRNAseq_workshop_3.html
#
#install.packages("msigdbr")
#library(msigdbr)
#BiocManager::install("enrichplot")
#BiocManager::install("fgsea")
#
# gse
# 
gene_list <- function(perturbed_gene_HEA,ms_ref){
  perturbed_gene_HEA_ref <-match(rownames(perturbed_gene_HEA),ms_ref$gene_short_name)
  perturbed_gene_HEA$entrez <- ms_ref[perturbed_gene_HEA_ref,]$entrez_annotation
  perturbed_gene_HEA <- perturbed_gene_HEA[!is.na(perturbed_gene_HEA$entrez),]
  gene_list_log2fc<- unlist(perturbed_gene_HEA$avg_log2FC)
  names(gene_list_log2fc) <-as.character(perturbed_gene_HEA$entrez)
  gene_list_log2fc <- gene_list_log2fc[order(gene_list_log2fc,decreasing = T)]
  return(gene_list_log2fc)
}
#pbmc.markers$entrez <- ms_ref[ms_ref$gene_short_name %in% pbmc.markers$gene,]$entrez_annotation

#variable_genes_entrez <- ms_ref[ms_ref$gene_short_name %in% pbmc.markers$gene,]
#pbmc.markers[variable_genes_entrez$gene_short_name,]$entrez <-
#  variable_genes_entrez[variable_genes_entrez$gene_short_name,]$entrez_annotation

gene_list_log2fc <- gene_list(hepa.markers[hepa.markers$cluster==0 
#                                             &hepa.markers$avg_log2FC > 1 |
#                                             hepa.markers$avg_log2FC < -1
                                           ,],ms_ref)
#
# try BP: biological process, CC: cellular component, or MF: molecular function
gse_result<- gseGO(geneList     = gene_list_log2fc,
                   OrgDb        = org.Mm.eg.db,
                   ont          = "BP",
                   minGSSize    = 12,
                   pvalueCutoff = 0.4,
                   pAdjustMethod = "BH",
                   verbose      = FALSE)

ridgeplot(gse_result,showCategory = 30)
View(gse_result@result)

#d <- GOSemSim::godata("org.Hs.eg.db", ont = "BP")    
#compare_cluster_GO_emap <- enrichplot::pairwise_termsim(gse_result, semData = d,  method="Wang")
#emapplot(compare_cluster_GO_emap, showCategory = 8)

#
# kegg
#
kegg_organism = "mmu"
kk <- gseKEGG(geneList     = gene_list_log2fc,
              organism     = kegg_organism,
              nPerm=10000,
              minGSSize    = 3,
              maxGSSize    = 800,
              pvalueCutoff = 0.05,
              pAdjustMethod = "none",
              keyType       = "ncbi-geneid")
ridgeplot(kk)
library("enrichplot")
gseaplot2(kk, geneSetID =1, title = kk$Description[1])
gseaplot2(kk, geneSetID =1, title = kk$Description[6])
#
#
#
path_index <- 1
View(kk@result)
kk$Description[path_index]
#
upregulated_entrez <- strsplit(kk$core_enrichment[path_index], split = "/")
ms_ref[ms_ref$entrez_annotation %in% upregulated_entrez[[1]],]$gene_short_name
#browseKEGG(kk, kk$ID[path_index])
library("pathview")
mah <- pathview(gene.data  = unlist(upregulated_entrez),
                     pathway.id = kk$ID[path_index],
                     species    = kegg_organism,
                     limit      = list(gene=max(abs(as.numeric(unlist(upregulated_entrez)))), cpd=1))
#
#
#
#
kk <- gseMKEGG(gene_list_log2fc,
               organism     = kegg_organism, nPerm=10000)
ridgeplot(kk)
#gseaplot(kk, geneSetID = 1, by = "runningScore", title = kk$Description[1])
gseaplot2(kk, geneSetID =1, title = kk$Description[1])
#
#
#
#browseKEGG(kk,)
