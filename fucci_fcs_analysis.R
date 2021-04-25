#BiocManager::install("flowUtils")
library(flowCore)
library(flowUtils)
library(ggplot2)
library(reshape2)
library(ggbeeswarm)
fcs_dir <- "/home/samba/pihome/2021/Kaneko/1_Data/2_Microscopy/20210316/20210421_Output/No1_Control"
ctl_fcsdata <- fucci_fcs_io(fcs_dir)
ctl_cellnames <- paste0("CTL-HeLaFucci-",ctl_fcsdata$ImageNumber,"-",ctl_fcsdata$ObjectNumber)
rownames(ctl_fcsdata) <- ctl_cellnames

ggplot(ctl_fcsdata,aes(x=G,y=R,color=time))+geom_point()+
  scale_y_log10(limits=c(10,1000))+scale_x_log10(limits=c(10,1000))
ggplot(ctl_fcsdata,aes(x=factor(time),y=R,color=gate))+geom_violin()

merge_fl_data <- ctl_fcsdata

source('./fucci_fl_monocle_pseudotime.R')



merge_fl_data$pseudotime <- fucci$Pseudotime
ggplot(merge_fl_data,aes(x=time,y=pseudotime,color=gate))+geom_point()

merge_fl_data_melt <-
  melt(merge_fl_data, id.vars=c("time","gate"),
       measure.vars = c("pseudotime"))


ggplot(merge_fl_data_melt,aes(x=factor(time),y=value, color=gate))+ geom_quasirandom(size=0.5)
#  geom_point(position = position_jitter(width = 0.4,height=0),size=0.01)
  
  #geom_violin()+ylab("pseudotime")

#beeswarm()

#end of 

#
# scatter plot
#
cellcycle_prog_mean_sd <- merge_fl_data %>% group_by(time) %>% summarize(mean = mean(pseudotime), sd = sd(pseudotime))

ggplot(cellcycle_prog_mean_sd,aes(x=time,y=mean))+geom_point(size=3)+geom_line()+
  ylab("Cell cycle progression")+ylim(c(-1,4))+
  geom_errorbar(aes(ymax = (mean + sd), ymin = (mean - sd)), width=2,
                position=position_dodge(0))


gridExtra::grid.arrange(p1,p2, p3, p4, nrow = 4)