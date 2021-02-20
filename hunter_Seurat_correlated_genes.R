# WGCNA package
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/InstallationInstructions.html
# "matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"
# "GO.db", "preprocessCore", "impute"
#BiocManager::install("WGCNA")
library(WGCNA)
library(Seurat)
library(VennDiagram)

options(stringsAsFactors = FALSE)
datExpr <- data.frame(t(GetAssayData(object=AML[["RNA"]], slot="scale.data")))
dim(datExpr)

sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


traitData <- t(GetAssayData(object=AML[["ADT"]], slot="scale.data"))
dim(traitData)

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(traitData, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(traitColors),
                    main = "Sample dendrogram and trait heatmap")

datTraits <- data.frame(traitData[,c(2:7)])
#rownames(datTraits) <- rownames(traitData)

source("/home/watson/public/shintaku/HUNTER/hunter_Seurat_WGCNA_network_auto.R")

source("/home/watson/public/shintaku/HUNTER/hunter_Seurat_WGCNA_relateModsToExt.R")

corr_module_normGFP <- names(datExpr)[moduleColors=="red"]
acorr_module_normGFP <- names(datExpr)[moduleColors=="pink"]
normGFP_module <- c(corr_module_normGFP,acorr_module_normGFP)

corr_moudle_mCherry <- names(datExpr)[moduleColors=="yellow"]
acorr_module_mCherry <- names(datExpr)[moduleColors=="greenyellow"]
mCherry_module <- c(corr_moudle_mCherry,acorr_module_mCherry)

PC_1_gene <- PCASigGenes(object=AML,pcs.use=1,pval.cut=0.1)
PC_2_gene <- PCASigGenes(object=AML,pcs.use=2,pval.cut=0.1)

gene_list <- list(normGFP=corr_module_normGFP, mCherry=mCherry_module,PC_1=PC_1_gene, PC_2=PC_2_gene)


venn.diagram(gene_list,filename = file.path(wdir,"gene.jpg"), fill=c(2,3,4,5), alpha=0.4, lty=3)

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

acorr_gene <- head(corr_order,100)
corr_gene <- tail(corr_order,100)

gene_list <- list(normGFP=normGFP_module,acorr=rownames(acorr_gene),corr=rownames(corr_gene))

venn.diagram(gene_list,filename = file.path(wdir,"gene.jpg"), fill=c(1,2,3), alpha=0.4, lty=3)

GOenr = GOenrichmentAnalysis(moduleColors, colnames(datExpr), organism = "mouse", nBestP = 10)
