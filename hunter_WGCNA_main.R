# WGCNA package
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/InstallationInstructions.html
# "matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"
# "GO.db", "preprocessCore", "impute"
#BiocManager::install("WGCNA")
library(WGCNA)
library(VennDiagram)

options(stringsAsFactors = FALSE)
datExpr <- data.frame(t(GetAssayData(object=AML[["RNA"]], slot="scale.data")))
dim(datExpr)

sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


traitData <- t(GetAssayData(object=AML[["ADT"]]))


datTraits <- data.frame(traitData[,c(2:7)])
dim(datTraits)

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitColors),
                    main = "Sample dendrogram and trait heatmap")

#
# search gene network modules
source("/home/watson/public/shintaku/HUNTER/hunter_WGCNA_network_auto.R")
#
# trait vs gene network modules
source("/home/watson/public/shintaku/HUNTER/hunter_WGCNA_relateModsToExt.R")
#
# visualize the hierarchy of modules
source("/home/watson/public/shintaku/HUNTER/hunter_WGCNA_network_visualize.R")


corr_module_normGFP <- names(datExpr)[moduleColors=="purple"]
acorr_module_normGFP <- names(datExpr)[moduleColors=="black"]
normGFP_module <- c(corr_module_normGFP,acorr_module_normGFP)

corr_moudle_mCherry <- names(datExpr)[moduleColors=="yellow"]
acorr_module_mCherry <- names(datExpr)[moduleColors=="green"]
mCherry_module <- c(corr_moudle_mCherry,acorr_module_mCherry)

PC_1_gene <- PCASigGenes(object=AML,pcs.use=1,pval.cut=0.1)
PC_2_gene <- PCASigGenes(object=AML,pcs.use=2,pval.cut=0.1)

gene_list<- list(normGFP=corr_module_normGFP, mCherry=mCherry_module,PC_1=PC_1_gene, PC_2=PC_2_gene)
venn.diagram(gene_list,filename = file.path(wdir,"gene.jpg"), fill=c(2,3,4,5), alpha=0.4, lty=3)

