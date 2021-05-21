source(file.path(rdir,"fucci_cellcycle_genes.R"))
sub_ref <- ms_ref %>%
  dplyr::filter(gene_short_name %in% rownames(pbmc))
genes <- fucci_cellcycle_genes(sub_ref)
cell_cycle_markers<-genes[[1]]
s_genes <- genes[[2]]
g2m_genes <- genes[[3]]
pbmc <- CellCycleScoring(pbmc,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score"), features = c(s_genes, g2m_genes))
features=c("S.Score","G2M.Score")
#FeaturePlot(pbmc, features = features,reduction = "pca")
p1 <- FeatureScatter(pbmc,feature1 = "G2M.Score",feature2 = "Venus",group.by = "cell")+scale_y_log10()
p2 <- FeatureScatter(pbmc,feature1 = "S.Score",feature2 = "Venus",group.by = "cell")+scale_y_log10()
p1+p2
FeatureScatter(pbmc,feature1 = "S.Score",feature2 = "G2M.Score",group.by = "gate")
