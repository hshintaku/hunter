source(file.path(rdir,"util/fucci_cellcycle_genes.R"))
sub_ref <- ms_ref %>%
  dplyr::filter(gene_short_name %in% rownames(pbmc))
genes <- fucci_cellcycle_genes(sub_ref)
cell_cycle_markers<-genes[[1]]
s_genes <- genes[[2]]
g2m_genes <- genes[[3]]
tig <- CellCycleScoring(tig,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
tig<- ScaleData(tig, vars.to.regress = c("S.Score", "G2M.Score"), features = c(s_genes, g2m_genes))
features=c("S.Score","G2M.Score")
#FeaturePlot(pbmc, features = features,reduction = "pca")
p1 <- FeatureScatter(tig,feature1 = "G2M.Score",feature2 = "FLD004",group.by = "cell")+scale_y_log10()
p2 <- FeatureScatter(tig,feature1 = "S.Score",feature2 = "FLD004",group.by = "cell")+scale_y_log10()
p1+p2
FeatureScatter(tig,feature1 = "S.Score",feature2 = "G2M.Score",group.by = "gate")
rm(p1,p2,sub_ref,genes,cell_cycle_markers)

