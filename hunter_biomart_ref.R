library(biomaRt)
#hs_mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ms_mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
#pig_mart <- useMart(biomart="ensembl", dataset="sscrofa_gene_ensembl")

func.biomart.ref <- function(hs_mart, gene_list,symbol){
  reference=getBM(attributes=c("ensembl_gene_id","description",symbol,"gene_biotype","chromosome_name"),
               filters="ensembl_gene_id",values=gene_list,mart=hs_mart)
  
  colnames(reference) <- c("ensembl_gene_id","description","gene_short_name","gene_biotype","chromosome_name")
  reference <- reference[!duplicated(reference$ensembl_gene_id),]
  row.names(reference) <- reference$ensembl_gene_id
  empty_logical <- which(reference$gene_short_name=="")
  reference$gene_short_name <- as.character(reference$gene_short_name)
  reference[empty_logical,]$gene_short_name <- as.character(reference[empty_logical,]$ensembl_gene_id)
  
  return(reference)
}
func.ref.intron <- function(hs_ref){
hs_ref_intron <- hs_ref
hs_ref_intron$ensembl_gene_id <- paste0(hs_ref$ensembl_gene_id,"_intron")
rownames(hs_ref_intron)<- hs_ref_intron$ensembl_gene_id
return(hs_ref_intron)
}

gene_list <- unique(data.frame(str_replace(allData$gene,"_intron","")))
colnames(gene_list) <- "gene"

#hs_ref <- func.biomart.ref(hs_mart,gene_list,"hgnc_symbol")
ms_ref <- unique(func.biomart.ref(ms_mart,gene_list,"mgi_symbol"))


#pig_ref <- func.biomart.ref(pig_mart,gene_list,"hgnc_symbol")

#hs_ref$gene_short_name <- paste0("hs-",hs_ref$gene_short_name)
#ms_ref$gene_short_name <- paste0("ms-",ms_ref$gene_short_name)

#pig_ref$gene_short_name <- paste0("pg-",pig_ref$gene_short_name)

#hs_ref_intron <- func.ref.intron(hs_ref)
#ms_ref_intron <- func.ref.intron(ms_ref)
#pig_ref_intron <- func.ref.intron(pig_ref)
#all_ref <- rbind(hs_ref,hs_ref_intron,ms_ref,ms_ref_intron,pig_ref,pig_ref_intron)
#all_ref <- rbind(hs_ref,ms_ref)
