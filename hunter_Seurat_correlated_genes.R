
#
#
# correlation at gene level
#
#
library('corrr')
fulldatExpr <- cbind(datExpr,datTraits)
x <- correlate(fulldatExpr)

corr_row <- data.frame(x[,"normGFP"])
genes <- colnames(fulldatExpr)
order_index <- order(corr_row)
corr_order <- data.frame(corr_row[order_index,])
rownames(corr_order) <- genes[order_index]
corr_order <- na.omit(corr_order)
write.csv(corr_order, file.path(wdir,'correlated_genes.csv'))

acorr_gene <- head(corr_order,20)
corr_gene <- tail(corr_order,20)

gene_list <- list(normGFP=normGFP_module,acorr=rownames(acorr_gene),corr=rownames(corr_gene))

venn.diagram(gene_list,filename = file.path(wdir,"gene.jpg"), fill=c(1,2,3), alpha=0.4, lty=3)

