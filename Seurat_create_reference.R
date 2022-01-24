
source(file.path(rdir,'util/hunter_biomart_ref.R'))
gene_list <- data.frame(rownames(pbmc))
colnames(gene_list) <- "gene"
filter="rgd_symbol"
symbol="rgd_symbol"
rdg_mart <- useMart(biomart="ensembl", dataset="rnorvegicus_gene_ensembl")
rdg_ref <- unique(func.biomart.ref(rdg_mart,gene_list,filter,symbol))
filter="hgnc_symbol"
symbol="hgnc_symbol"
hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
hs_ref <- unique(func.biomart.ref(hs_mart,gene_list,filter,symbol))


ms_ref <- rbind(hs_ref,rdg_ref)

missing_ref <- subset(gene_list,!(gene %in% ms_ref$gene_short_name))
adding_ref <- data.frame(cbind(missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene,missing_ref$gene))
colnames(adding_ref) <- colnames(ms_ref)
rownames(adding_ref) <- adding_ref$ensembl_gene_id
ms_ref <- rbind(adding_ref,ms_ref)
rm(missing_ref,adding_ref)
