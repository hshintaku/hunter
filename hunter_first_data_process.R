# load functions for barcode decoding
source(file.path(rdir,"util/whitelist_encode.R"))
# laod whitelist and check the batch effect
source(file.path(rdir,'preprocess/preprocess_whitelist.R'))
active_barcode <- barcode[sort(unique(allencoded$first_index)),]

# preprocess the count data and load reference
source(file.path(rdir,'preprocess/preprocess_RNAseq_data.R'))
#
# download reference data from ensembl with biomaRt
gene_list <- unique(data.frame(str_replace(allData$gene,"_intron","")))
colnames(gene_list) <- "gene"
#
source(file.path(rdir,'util/hunter_biomart_ref.R'))
#hs_ref <- func.biomart.ref(hs_mart,gene_list,"hgnc_symbol")

filter="ensembl_gene_id"
if (symbol=="mgi_symbol"){
  ms_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(ms_mart,gene_list,filter,symbol))
}else{
  hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  ms_ref <- unique(func.biomart.ref(hs_mart,gene_list,filter,symbol))
}
  
#
#pig_mart <- useMart(biomart="ensembl", dataset="sscrofa_gene_ensembl")


missing_ref <- subset(gene_list,!(gene %in% ms_ref$ensembl_gene_id))
adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
colnames(adding_ref) <- colnames(ms_ref)
rownames(adding_ref) <- adding_ref$ensembl_gene_id
ms_ref <- rbind(adding_ref,ms_ref)
rm(missing_ref,adding_ref)

# save count data with 10x format
source(file.path(rdir, 'preprocess/preprocess_save_10x_format.R'))

