#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
#            "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
# https://github.com/aertslab/SCENIC
library(SCopeLoomR)
library(SCENIC)
library(loomR)
library(conflicted)
conflict_prefer("first", "S4Vectors")
conflict_prefer("finalize", "SCopeLoomR")
exprMat <- as.matrix(hepa_all[["RNA"]]@counts)
cellInfo <- data.frame(seuratCluster=Idents(hepa_all))
#exprMat <- Seurat::Read10X(data.dir="/home/samba/public/shintaku/20211124HiSeqX006_hunter/")
loom <- build_loom("hunter_hepa.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

### Initialize settings
#library(SCENIC)
org <- "mgi" # or hgnc, or dmel
dbDir <- "/home/samba/sanger/genome/cisTarget_mm10/" # RcisTarget databases location
myDatasetTitle <- "SCENIC hunter hepatocytes" # choose a name for your analysis
dbs <- c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
         "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(dbs)<-c("500bp","10kb")
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=16)
saveRDS(scenicOptions, file="scenicOptions.Rds") 

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="scenicOptions.Rds")

# Optional: Binarize activity
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="scenicOptions.Rds") 

### Exploring output 
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Elk3"]
viewMotifs(tableSubset, options=list(pageLength=1)) 

# output/Step2_regulonTargetsInfo.tsv in detail: 
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Elk3" & highConfAnnot==TRUE]
viewMotifs(tableSubset, options=list(pageLength=1)) 

# Cell-type specific regulators (RSS): 
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "seuratCluster"], )
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
