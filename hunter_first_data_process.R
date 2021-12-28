# load functions for barcode decoding
source(file.path(rdir,"util/whitelist_encode.R"))
# laod whitelist and check the batch effect
source(file.path(rdir,'preprocess/preprocess_whitelist.R'))

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
symbol="rgd_symbol"
rdg_mart <- useMart(biomart="ensembl", dataset="rnorvegicus_gene_ensembl")
rdg_ref <- unique(func.biomart.ref(rdg_mart,gene_list,filter,symbol))
symbol="hgnc_symbol"
hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
hs_ref <- unique(func.biomart.ref(hs_mart,gene_list,filter,symbol))


ms_ref <- rbind(hs_ref,rdg_ref)

missing_ref <- subset(gene_list,!(gene %in% ms_ref$ensembl_gene_id))
adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
colnames(adding_ref) <- colnames(ms_ref)
rownames(adding_ref) <- adding_ref$ensembl_gene_id
ms_ref <- rbind(adding_ref,ms_ref)
rm(missing_ref,adding_ref)

# save count data with 10x format
source(file.path(rdir, 'preprocess/preprocess_save_10x_format.R'))

