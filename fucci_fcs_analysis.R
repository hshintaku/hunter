#BiocManager::install("flowUtils")
library(tidyverse)
library(flowCore)
library(flowUtils)
library(ggplot2)
library(reshape2)
library(ggbeeswarm)
library(dplyr)
library(monocle)
fcs_dir <- "/home/samba/pihome/2021/Kaneko/1_Data/2_Microscopy/20210316/20210421_Output/No1_Control"
ctl_fcsdata <- fucci_fcs_io(fcs_dir,"CTL")
ctl_cellnames <- paste0("CTL-HeLaFucci-",ctl_fcsdata$ImageNumber,"-",ctl_fcsdata$ObjectNumber)
rownames(ctl_fcsdata) <- ctl_cellnames

fcs_dir <- "/home/samba/pihome/2021/Kaneko/1_Data/2_Microscopy/20210316/20210421_Output/No2_ELP/"
epl_fcsdata <- fucci_fcs_io(fcs_dir,"EPL")
epl_cellnames <- paste0("EPL-HeLaFucci-",epl_fcsdata$ImageNumber,"-",epl_fcsdata$ObjectNumber)
rownames(epl_fcsdata) <- epl_cellnames

p1 <- ggplot(ctl_fcsdata,aes(x=G,y=R,color=time))+geom_point()+
  scale_y_log10(limits=c(10,1000))+scale_x_log10(limits=c(10,1000))
ggplot(ctl_fcsdata,aes(x=factor(time),y=R,color=gate))+geom_violin()

merge_fl_data <- rbind(ctl_fcsdata,epl_fcsdata)

source('./fucci_fl_monocle_pseudotime.R')

fl_fucci <- fucci_fl_monocle_pseudotime(merge_fl_data)
plot_cell_trajectory(fl_fucci, color_by = "time")


merge_fl_data$pseudotime <- fl_fucci$Pseudotime
#ggplot(merge_fl_data,aes(x=time,y=pseudotime,color=gate))+geom_point()

merge_fl_data_melt <-
  melt(merge_fl_data, id.vars=c("time","gate"),
       measure.vars = c("pseudotime"))


p2 <- ggplot(merge_fl_data_melt,aes(x=factor(time),y=value, color=gate))+ geom_quasirandom(size=0.5)
#p1+p2
gridExtra::grid.arrange(p1,p2,  nrow = 2)
#end of 

#
# scatter plot
#
cellcycle_prog_mean_sd <- merge_fl_data %>% group_by(time,exp) %>% summarize(mean = mean(pseudotime), sd = sd(pseudotime))

ggplot(cellcycle_prog_mean_sd,aes(x=time,y=mean,color=exp))+geom_point(size=3)+geom_line()+
  ylab("Cell cycle progression")+ylim(c(-1,4))+
  geom_errorbar(aes(ymax = (mean + sd), ymin = (mean - sd)), width=2,
                position=position_dodge(0))


#gridExtra::grid.arrange(p1,p2, p3, p4, nrow = 4)