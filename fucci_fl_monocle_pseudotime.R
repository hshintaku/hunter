
# compute pseudotime to summarize the cell cycle progression
library(monocle)
cell_metadata_fl <- new("AnnotatedDataFrame", data = merge_fl_data[,c(1,6,7)])
gene_annotation_fl <- data.frame(c("R","G"))
colnames(gene_annotation_fl) <- "gene_short_name"
rownames(gene_annotation_fl) <- c("mCherry","Venus")
gene_annotation_fl <- new("AnnotatedDataFrame", data = gene_annotation_fl)
expression_matrix_fl <- t(merge_fl_data[,c(4,5)])
rownames(expression_matrix_fl) <- c("mCherry","Venus")

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
fucci <- orderCells(fucci)
plot_cell_trajectory(fucci, color_by = "gate")
plot_genes_in_pseudotime(fucci, color_by = "time")
