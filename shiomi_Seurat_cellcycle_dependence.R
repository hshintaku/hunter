source(file.path(rdir,"util/fucci_cellcycle_genes.R"))

tig <-pbmc#seurat_object

sub_ref <- ms_ref %>%
  dplyr::filter(gene_short_name %in% rownames(tig))
genes <- fucci_cellcycle_genes(sub_ref)
cell_cycle_markers<-genes[[1]]
s_genes <- genes[[2]]
g2m_genes <- genes[[3]]
tig <- CellCycleScoring(tig,g2m.features = g2m_genes,s.features = s_genes,set.ident = TRUE)
tig$CC.Difference <- tig$S.Score - tig$G2M.Score
#tig<- ScaleData(tig, vars.to.regress = c("S.Score", "G2M.Score"), features = c(s_genes, g2m_genes))

tig<- ScaleData(tig, vars.to.regress = "CC.Difference", features = c(s_genes, g2m_genes))
#features=c("S.Score","G2M.Score")
#FeaturePlot(pbmc, features = features,reduction = "pca")
p1 <- FeatureScatter(tig,feature1 = "G2M.Score",feature2 = "FLD004")+scale_y_log10()
p2 <- FeatureScatter(tig,feature1 = "S.Score",feature2 = "FLD004")+scale_y_log10()
p1+p2

FeatureScatter(tig,feature1 = "S.Score",feature2 = "G2M.Score")

tig <- FindVariableFeatures(tig, selection.method = "vst", nfeatures = 1000)
tig <- RunPCA(tig, npcs=50, features = VariableFeatures(object = tig))
plot1 <- VariableFeaturePlot(tig)
plot1 <- LabelPoints(plot = plot1, points = top10)
plot1
p1<-DimPlot(tig,reduction="umap")
p2<-FeaturePlot(tig,features="fld_FLD500",reduction="umap")
p1+p2
#pbmc<-tig

VlnPlot(tig,features='fld_FLDtotal',group.by = 'Phase')

t20ctl <- subset(tig,subset = condition=="T20CTL")
t50ctl <- subset(tig,subset = condition=="T50CTL")
tazctl <- subset(tig,subset=condition=="TAZCTL")

p1<-DimPlot(t20ctl,group.by = "Phase")
p2<-DimPlot(t50ctl,group.by = "Phase")
p3<-DimPlot(tazctl,group.by = "Phase")
p1+p2+p3
