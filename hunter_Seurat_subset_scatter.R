

p1 <- DimPlot(AML, reduction = "pca",group.by = "dish")
p2 <- DimPlot(AML, reduction = "pca",group.by = "gate")
p3 <- FeaturePlot(AML,features=c("Venus"),min.cutoff = 'q05',max.cutoff = 'q95',reduction="pca")
p4 <- FeaturePlot(AML,features=c("normGFP"),min.cutoff = 'q05',max.cutoff = 'q95',reduction="pca")
p5 <- FeaturePlot(AML,features=c("Ifi27l2a"),min.cutoff = 'q05',max.cutoff = 'q95',reduction="pca")
p6 <- FeaturePlot(AML,features=c("Spp1"),min.cutoff = 'q05',max.cutoff = 'q95',reduction="pca")
p7 <- FeaturePlot(AML,features=c("Lmcd1"),min.cutoff = 'q05',max.cutoff = 'q95',reduction="pca")
p8 <- FeaturePlot(AML,features=c("Ifit3"),min.cutoff = 'q05',max.cutoff = 'q95',reduction="pca")
grid.arrange(p1, p2, p3, p4,p5,p6,p7,p8, nrow = 4)
grid.arrange(p3,p4,nrow=1)

# correlation in PC_2 vs GFP
p1 <- FeatureScatter(object=AML, feature1='Venus',feature2='Lmcd1',group.by = 'dish')
p2 <- FeatureScatter(object=AML, feature1='Venus',feature2='Ccn2',group.by = 'dish')
p3 <- FeatureScatter(object=AML, feature1='Venus',feature2='Ankrd1',group.by = 'dish')
p4 <- FeatureScatter(object=AML, feature1='Venus',feature2='Ifit3',group.by = 'dish')
p5 <- FeatureScatter(object=AML, feature1='Venus',feature2='Ifit1',group.by = 'dish')
p6 <- FeatureScatter(object=AML, feature1='Venus',feature2='Ifit3b',group.by = 'dish')
grid.arrange(p1, p2, p3, p4,p5,p6, nrow = 2)

# correlation in PC_2 vs GFP
p1 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Lmcd1',group.by = 'dish')
p2 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Ccn2',group.by = 'dish')
p3 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Ankrd1',group.by = 'dish')
p4 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Ifit3',group.by = 'dish')
p5 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Ifit1',group.by = 'dish')
p6 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Ifit3b',group.by = 'dish')
grid.arrange(p1, p2, p3, p4,p5,p6, nrow = 2)

# correlation in PC_1 vs GFP
p1 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Ifi27l2a',group.by = 'dish')
p2 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Ube2c',group.by = 'dish')
p3 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Scd1',group.by = 'dish')
p4 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Spp1',group.by = 'dish')
p5 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Peg3',group.by = 'dish')
p6 <- FeatureScatter(object=AML, feature1='normGFP',feature2='Tmem176a',group.by = 'dish')
grid.arrange(p1, p2, p3, p4,p5,p6, nrow = 2)

# correlation in PC_1 vs GFP
p1 <- FeatureScatter(object=AML, feature1='Venus',feature2='Ifi27l2a',group.by = 'dish')
p2 <- FeatureScatter(object=AML, feature1='Venus',feature2='Ube2c',group.by = 'dish')
p3 <- FeatureScatter(object=AML, feature1='Venus',feature2='Scd1',group.by = 'dish')
p4 <- FeatureScatter(object=AML, feature1='Venus',feature2='Spp1',group.by = 'dish')
p5 <- FeatureScatter(object=AML, feature1='Venus',feature2='Peg3',group.by = 'dish')
p6 <- FeatureScatter(object=AML, feature1='Venus',feature2='Tmem176a',group.by = 'dish')
grid.arrange(p1, p2, p3, p4,p5,p6, nrow = 2)

# correlation in PC_1 vs mCherry
p1 <- FeatureScatter(object=AML, feature1='mCherry',feature2='Ifi27l2a',group.by = 'dish')
p2 <- FeatureScatter(object=AML, feature1='mCherry',feature2='Ube2c',group.by = 'dish')
p3 <- FeatureScatter(object=AML, feature1='mCherry',feature2='Scd1',group.by = 'dish')
p4 <- FeatureScatter(object=AML, feature1='mCherry',feature2='Spp1',group.by = 'dish')
p5 <- FeatureScatter(object=AML, feature1='mCherry',feature2='Peg3',group.by = 'dish')
p6 <- FeatureScatter(object=AML, feature1='mCherry',feature2='Tmem176a',group.by = 'dish')
grid.arrange(p1, p2, p3, p4,p5,p6, nrow = 2)