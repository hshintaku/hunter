source(file.path(rdir,"util/fucci_cellcycle_genes.R"))
sub_ref <- ms_ref %>%
  dplyr::filter(gene_short_name %in% rownames(pbmc))
mm_url<-"https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv"
#hs_url<-"https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv"
genes <- fucci_cellcycle_genes(sub_ref,mm_url)
cell_cycle_markers<-genes[[1]]
s_genes <- genes[[2]]
g2m_genes <- genes[[3]]
pbmc <- CellCycleScoring(pbmc,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
#pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score"), features = c(s_genes, g2m_genes))
features=c("S.Score","G2M.Score")
#FeaturePlot(pbmc, features = features,reduction = "pca")
#p1 <- FeatureScatter(pbmc,feature1 = "G2M.Score",feature2 = "Venus",group.by = "cell")+scale_y_log10()
#p2 <- FeatureScatter(pbmc,feature1 = "S.Score",feature2 = "Venus",group.by = "cell")+scale_y_log10()
#p1+p2
FeatureScatter(pbmc,feature1 = "S.Score",feature2 = "G2M.Score")
rm(p1,p2,sub_ref,genes,cell_cycle_markers)
VlnPlot(pbmc,features=c("G2M.Score","S.Score"),group.by = "plate")

p1<-DimPlot(pbmc)
p2<-FeaturePlot(pbmc,features = "G2M.Score")
p3<-FeaturePlot(pbmc,features = "S.Score")
p1+p2+p3

