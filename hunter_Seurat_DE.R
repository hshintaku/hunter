library(Seurat)
p16<-subset(hepa,subset=plate=="P16")
p19<-subset(hepa,subset=plate=="P19")
gfpm<-subset(hepa,subset=plate=="P18")
cntl<-subset(hepa,subset=plate=="P17")
p15<- subset(allcell,subset=plate=="P15")
p20<-subset(allcell,subset=plate=="p20")
Idents(object=p15) <- "Cancer_vivo"
Idents(object=p20) <- "Cancer_vitro"
cancer<- merge(p15,y=p20)
gfpp <- merge(p16,y=p19)
Idents(object = gfpp) <- "GFP+"
Idents(object = gfpm) <- "GFP-"
Idents(object = cntl) <- "cntl"
gfp_pm <- merge(gfpp,y=gfpm)
hepa_all <- merge(gfp_pm,y=cntl)
cell_all <- merge(hepa_all,y=cancer)
# Find differentially expressed features 
hepa.de.markers <- FindMarkers(hepa_all, ident.1 = "GFP+", ident.2 = "cntl")
# view results
hepa.de.markers$test <- "FALSE"
hepa.de.markers[abs(hepa.de.markers$avg_log2FC)>1 & hepa.de.markers$p_val_adj<0.01,]$test <-"TRUE"
hepa.de.markers$gene <- rownames(hepa.de.markers)
#View(hepa.de.markers)
#ggplot(hepa.de.markers,aes(y=-log10(p_val_adj),x=avg_log2FC,color=test))+geom_point()
library(pheatmap)
pheatmap(hepa_all[["RNA"]]@data[hepa.de.markers[hepa.de.markers$test==TRUE,]$gene,],
         annotation_col = hepa_all[["plate"]],
         cluster_cols = F)

