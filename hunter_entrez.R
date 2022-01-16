cluster_marker_entrez <- function(pbmc.markers,ms_ref,cluster.num,p.thres,log2FC.thres){
  cluster.marker.n <- pbmc.markers[pbmc.markers$cluster==cluster.num &
                                     pbmc.markers$p_val_adj<p.thres &
                                     pbmc.markers$avg_log2FC>log2FC.thres,]
  perturbed_gene.n <- ms_ref[ms_ref$gene_short_name %in% cluster.marker.n$gene,]
  perturbed_gene.n <- perturbed_gene.n[!is.na(perturbed_gene.n$entrez_annotation),]
  return(perturbed_gene.n)
  #return(cluster.marker.n)
}
#
# GOenrichmentAnalysis (experimental)
#BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")
#BiocManager::install("org.Hs.eg.db")
#library("org.Hs.eg.db")
library(clusterProfiler)
#
#ensmusg <- data.frame(unlist(as.list(org.Hs.egENSEMBL2EG)))
ensmusg <- data.frame(unlist(as.list(org.Mm.egENSEMBL2EG)))
# entrez annotation
ms_ref$entrez_annotation <- ensmusg[ms_ref$ensembl_gene_id,]

