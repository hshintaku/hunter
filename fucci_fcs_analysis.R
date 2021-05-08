#BiocManager::install("flowUtils")
library(tidyverse)
library(flowCore)
library(flowUtils)
library(ggplot2)
library(reshape2)
library(ggbeeswarm)
library(dplyr)
library(monocle)
library(ggpubr)
library(cowplot)
source("./fucci_fcs_io.R")
source("./geom_split_violin.R")

# load timelapse fcs data taken by Kaneko
source('./fucci_fcs_load_microscopy_data.R')

# compute pseudotime with monocle 
fl_fucci2 <- fucci_fl_monocle_pseudotime(merge_fl_data,"F")
merge_fl_data$pseudotime <- fl_fucci2$Pseudotime
rm(fl_fucci2)

cellcycle_prog_mean_sd <- merge_fl_data %>%
  dplyr::group_by(time,exp) %>%
  dplyr::summarize(mean_p = mean(pseudotime), sd_p = sd(pseudotime))

#output computed data
write.csv(cellcycle_prog_mean_sd,"/home/samba/storage0/shintaku/20210324MiSeq016_Kaneko/Fucci2_microscopy_ctl_pseudotime.csv")
# plot gate vs 

merge_fl_data_melt <-
  melt(merge_fl_data, id.vars=c("condition","gate"),
       measure.vars = c("pseudotime"))

# P0 <-ggplot(merge_fl_data,aes(x=factor(gate),y=pseudotime, fill=gate))+
#   geom_violin()+
#   scale_x_discrete(limits = c("lateSG2M", "MG1", "G1","S"))+
#   scale_fill_manual(values=c(gate="red","green","grey","orange"))+
#   ylab("cell cycle progression")
yplot<-ggdensity(merge_fl_data, "pseudotime", fill = "gate")+
  scale_fill_manual(values=c(gate="red","green","grey","orange"))+
  rotate()+theme(axis.text.x=element_blank(),
                 axis.text.y = element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank())+
  xlim(c(-1,4))


ggplot(merge_fl_data,aes(x=condition,y=pseudotime, color=gate))+
  geom_quasirandom(size=0.5)
ggplot(merge_fl_data,aes(x=condition,y=pseudotime, color=gate))+
  geom_split_violin()+
  scale_color_manual(values=c(gate="red","green","grey","orange"))
  
# plot time evolution of cell cycle progression defined by pseudotime
P1 <- ggplot(cellcycle_prog_mean_sd,aes(x=time,y=mean_p,color=exp))+geom_point(size=3)+geom_line()+
  ylab("Cell cycle progression")+ylim(c(-1,4))+
  geom_errorbar(aes(ymax = (mean_p + sd_p), ymin = (mean_p - sd_p)), width=2,
                position=position_dodge(0))  + 
  theme(legend.position="top")

P1 <- plot_grid(NULL, NULL, P1, yplot, ncol = 2, align = "hv", 
          rel_widths = c(2,0.5), rel_heights = c(0,1))

# summarize the data for plotting
cellcycle_protein_summary <- merge_fl_data %>%
  dplyr::group_by(condition,time,exp,gate) %>%
  dplyr::summarise(number = length(pseudotime))

cellcycle_protein_spread <- cellcycle_protein_summary %>%
  tidyr::spread(key=gate,value=number,fill=0)
cellcycle_protein_spread$total <- 
  cellcycle_protein_spread$G1+cellcycle_protein_spread$lateSG2M+cellcycle_protein_spread$MG1+cellcycle_protein_spread$S
cellcycle_protein_spread$G1 <- 
  cellcycle_protein_spread$G1/cellcycle_protein_spread$total
cellcycle_protein_spread$lateSG2M <- 
  cellcycle_protein_spread$lateSG2M/cellcycle_protein_spread$total
cellcycle_protein_spread$MG1 <- 
  cellcycle_protein_spread$MG1/cellcycle_protein_spread$total
cellcycle_protein_spread$S <- 
  cellcycle_protein_spread$S/cellcycle_protein_spread$total

cellcycle_protein_summary <- melt( cellcycle_protein_spread,id.vars=c("time","exp"),measure.vars = c("G1","MG1","lateSG2M","S"))

# cellcycle_protein_summary$condition <- "ELP"
# cellcycle_protein_summary[cellcycle_protein_summary$exp %in% c("CTL1","CTL2"),]$condition <- "CTL"
cellcycle_protein_mean_sd <- cellcycle_protein_summary %>%
  dplyr::group_by(time,condition,variable) %>%
  dplyr::summarize(mean = mean(value), sd = sd(value))

P2 <- ggplot(cellcycle_protein_mean_sd,aes(x=time,y=mean,shape=condition,color=variable))+geom_point(size=3)+geom_line()+
  ylab("Fraction (%)")+ylim(c(-0.1,1))+
  geom_errorbar(aes(ymax = (mean + sd), ymin = (mean - sd)), width=2,
                position=position_dodge(0))+
    scale_color_manual(values=c(variable="red","gray","green","orange"))

rm(cellcycle_protein_mean_sd,cellcycle_protein_spread,cellcycle_protein_summary,cellcycle_protein_total)

gridExtra::grid.arrange(P1,P2,nrow = 2)

# old version of cell phase progression annotated by fcs
#library(openxlsx)
#cellcycle_protein <- read.xlsx("/home/samba/pihome/2021/Kaneko/1_Data/4_RNAseq/20210323_fl_cellcycle_score.xlsx")
#cellcycle_protein <- melt(cellcycle_protein, id.vars=c("time","exp"),measure.vars = c("G1","G1-S","lateS-G2-M"))
#ggplot(cellcycle_protein,aes(x=time,y=value,shape=exp,color=variable))+geom_point(size=3)+geom_line()+
#  scale_color_manual(values=c(variable="red","orange","green"))

#P2 <- ggplot(cellcycle_protein_summary,aes(x=time,y=value,color=variable,shape=exp))+geom_point(size=3)+geom_line()+
#  scale_color_manual(values=c(variable="red","gray","green","orange"))


# fcs_dir <- "/home/samba/pihome/2021/Shintaku/FucciSA2/"
# phase <- c("F0","F3","F4","F6","F2","F7","F1","F5")
# fucci2_fcsdata <- fucci_fcs_io(fcs_dir,"Fucci2",phase)
# fucci2_fcsdata <- subset(fucci2_fcsdata, gate!="F0")
# ggplot(fucci2_fcsdata,aes(x=Venus.A,y=mCherry.A,color=gate))+geom_point()+scale_x_log10()+scale_y_log10()
# #sawano_fcsdata <- read.FCS(fcs_dir)
# #write.FCS(sawano_fcsdata,'/home/samba/pihome/2021/Shintaku/150225_Fucci32_8rev.fcs')
# #sawano_expdata <- data.frame(sawano_fcsdata@exprs)


