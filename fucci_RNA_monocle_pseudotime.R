library(monocle)
#cell_metadata_fl <- new("AnnotatedDataFrame", data = merge_fl_data[,c(1,6)])
#gene_annotation_fl <- data.frame(c("R","G"))
#colnames(gene_annotation_fl) <- "gene_short_name"
#rownames(gene_annotation_fl) <- c("mCherry","Venus")


#gene_annotation_fl <- new("AnnotatedDataFrame", data = gene_annotation_fl)
#expression_matrix_fl <- t(merge_fl_data[,c(2,3)])*1000
#rownames(expression_matrix_fl) <- c("mCherry","Venus")

all.genes <- cell_cycle_markers$gene_short_name #rownames(fucci)
fucci <- ScaleData(fucci, features = all.genes)
fucci <- RunPCA(fucci, npcs=15, features = VariableFeatures(object = fucci))

# fucci <- JackStraw(fucci, num.replicate = 100)
# fucci <- ScoreJackStraw(fucci, dims = 1:7)
# JackStrawPlot(fucci, dims = 1:7)
# ElbowPlot(fucci)
# fucci <- FindNeighbors(fucci, dims = 1:3)
# fucci <- FindClusters(fucci, resolution = 0.01)
# cluster = as.numeric(Idents(fucci))
# fucci <- RunUMAP(fucci, dims = 1:2)

DimPlot(fucci, reduction = "pca",group.by = "exp")
fucci_monocle <- as.CellDataSet(fucci)

# fucci <- newCellDataSet(as.matrix(expression_matrix_fl),
#                         phenoData = cell_metadata_fl,
#                         featureData = gene_annotation_fl,
#                         expressionFamily = negbinomial.size())


fucci_monocle <- estimateSizeFactors(fucci_monocle)
fucci_monocle <- estimateDispersions(fucci_monocle)
disp_table <- dispersionTable(fucci_monocle)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)
ordering_genes <- data.frame(c(cell_cycle_markers$gene_short_name))
colnames(ordering_genes) <- "gene_id"

fucci_monocle <- setOrderingFilter(fucci_monocle, ordering_genes)
fucci_monocle <- reduceDimension(fucci_monocle)

fucci_monocle <- orderCells(fucci_monocle)
plot_cell_trajectory(fucci_monocle, color_by = "exp")
plot_genes_in_pseudotime(fucci_monocle, color_by = "exp")
