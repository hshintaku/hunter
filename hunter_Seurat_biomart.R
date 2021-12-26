
symbol="mgi_symbol"
#symbol="hgnc_symbol"
#symbol="rgd_symbol"

#filter="ensembl_gene_id"
filter="mgi_symbol"
# download reference data from ensembl with biomaRt
gene_list <- data.frame(unique(rownames(integrated)))
colnames(gene_list) <- "gene"
#
source(file.path(rdir,'util/hunter_biomart_ref.R'))
#hs_ref <- func.biomart.ref(hs_mart,gene_list,"hgnc_symbol")

if (symbol=="mgi_symbol"){
  ms_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(ms_mart,gene_list,filter,symbol))
}else if(symbol=="hgnc_symbol"){
  hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(hs_mart,gene_list,filter,symbol))
}else{
  rgd_mart <- useMart(biomart="ensembl", dataset="rnorvegicus_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(rgd_mart,gene_list,filter,symbol))
}

#
#pig_mart <- useMart(biomart="ensembl", dataset="sscrofa_gene_ensembl")


missing_ref <- subset(gene_list,!(gene %in% ms_ref$ensembl_gene_id))
adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
colnames(adding_ref) <- colnames(ms_ref)
rownames(adding_ref) <- adding_ref$ensembl_gene_id
ms_ref <- rbind(adding_ref,ms_ref)
rm(missing_ref,adding_ref)

source(file.path(rdir,"hunter_entrez.R"))