
# compute pseudotime to summarize the cell cycle progression
fucci_fcs_monocle_pseudotime <- function(merge_fl_data,TF){
library(monocle)
  col_name_list <- tolower(colnames(merge_fl_data))
  time_index <- as.numeric(match("time",col_name_list))
  exp_index <- as.numeric(match("exp",col_name_list))
  gate_index <- as.numeric(match("gate",col_name_list))
  cell_metadata_fl <- new("AnnotatedDataFrame", data = merge_fl_data[,c(time_index,exp_index,gate_index)])

  gene_annotation_fl <- data.frame(c("R","G"))
  colnames(gene_annotation_fl) <- "gene_short_name"
  rownames(gene_annotation_fl) <- c("mCherry","Venus")
  gene_annotation_fl <- new("AnnotatedDataFrame", data = gene_annotation_fl)


  mcherry_index <- as.numeric(match("mcherry.a",col_name_list))
  venus_index <- as.numeric(match("venus.a",col_name_list))
  # if (str_detect(norm_method,"linear")){
    expression_matrix_fl <- t(merge_fl_data[,c(mcherry_index,venus_index)])
  # }else{
  #   expression_matrix_fl <- data.frame(merge_fl_data[,c(mcherry_index,venus_index)])
  #   expression_matrix_fl <- log2(expression_matrix_fl)
  #   expression_matrix_fl <- t(expression_matrix_fl)
  #   
  # }
  
  rownames(expression_matrix_fl) <- c("mCherry","Venus")
  #expression_matrix_flmCherry <- t(scale(expression_matrix_fl$mCherry, center=FALSE)*1000)
  #expression_matrix_fl$Venus <- t(scale(expression_matrix_fl$Venus, center=FALSE)*1000)

  fucci <- newCellDataSet(as.matrix(expression_matrix_fl),
                        phenoData = cell_metadata_fl,
                        featureData = gene_annotation_fl,
                        expressionFamily = negbinomial.size())
  fucci <- estimateSizeFactors(fucci)
  fucci <- estimateDispersions(fucci)
  disp_table <- dispersionTable(fucci)
  ordering_genes <- subset(disp_table, mean_expression >= 0.1)
  fucci <- setOrderingFilter(fucci, ordering_genes)
  fucci <- reduceDimension(fucci)
  fucci <- orderCells(fucci,reverse=TF)
#plot_genes_in_pseudotime(fucci,color_cells_by = 'gate')
return(fucci)
}
