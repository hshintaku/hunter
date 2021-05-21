library(monocle)


pbmc_monocle <- as.CellDataSet(pbmc)


pbmc_monocle <- estimateSizeFactors(pbmc_monocle)
pbmc_monocle <- estimateDispersions(pbmc_monocle)
disp_table <- dispersionTable(pbmc_monocle)
ordering_genes <- subset(disp_table, mean_expression >= 1)
colnames(ordering_genes) <- "gene_id"

rownames(perturbed_gene_HEA)

pbmc_monocle <- setOrderingFilter(pbmc_monocle, ordering_genes)
pbmc_monocle <- reduceDimension(pbmc_monocle)

pbmc_monocle <- orderCells(pbmc_monocle)
p1<-plot_cell_trajectory(pbmc_monocle, color_by = "cell")
p2<-plot_cell_trajectory(pbmc_monocle, color_by = "gate")
p1+p2

plot_genes_in_pseudotime(pbmc_monocle, color_by = "cell")

pbmc[["Pseudotime"]] <- pbmc_monocle$Pseudotime

FeatureScatter(pbmc,feature1 = "Pseudotime",feature2 = "Venus",group.by = "cell")+scale_y_log10()
FeatureScatter(pbmc,feature1 = "Pseudotime",feature2 = "FLD070",group.by = "cell")+scale_y_log10()

FeatureScatter(pbmc,feature1 = "GC",feature2 = "FLD004",group.by = "cell")+scale_y_log10()
FeatureScatter(pbmc,feature1="Venus",feature2 = "FLD010",group.by = "cell")+scale_y_log10()+scale_x_log10()
