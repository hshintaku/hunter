

fcs_dir <- "/home/samba/pihome/2021/Shintaku/Fucci3.2/Live/"
phase <- c("F0","F6","F2","F5","F3","F4","F7","F1")
ctl_fcsdata <- fucci_fcs_io(fcs_dir,"Fucci3",phase)
ctl_fcsdata <- subset(ctl_fcsdata, gate!="F0")
ggplot(ctl_fcsdata,aes(x=Venus.A,y=mCherry.A,color=gate))+geom_point()+scale_x_log10()+scale_y_log10()

fl_fucci3 <- fucci_fl_monocle_pseudotime(ctl_fcsdata)
ctl_fcsdata$pseudotime <- fl_fucci3$Pseudotime
cellcycle_prog_mean_sd <- ctl_fcsdata %>%
  dplyr::group_by(gate) %>%
  dplyr::summarize(mean_p = mean(pseudotime), sd_p = sd(pseudotime))
write.csv(cellcycle_prog_mean_sd,"/home/samba/storage0/shintaku/20210324MiSeq016_Kaneko/Fucci3.2_facs_pseudotime.csv")


p1 <- ggplot(ctl_fcsdata,aes(x=G,y=R,color=time))+geom_point()+
  scale_y_log10(limits=c(10,1000))+scale_x_log10(limits=c(10,1000))
ggplot(ctl_fcsdata,aes(x=factor(time),y=R))+geom_violin()

merge_fl_data <- rbind(ctl_fcsdata[,c(7,17,19,20,21)],fucci2_fcsdata[,c(7,17,19,20,21)])

ctl_fcsdata <- subset(merge_fl_data,exp=="Fucci3")
fl_fucci32 <- subset(fl_fucci3,exp=="Fucci3")




plot_cell_trajectory(fl_fucci, color_by = "gate")
plot_genes_in_pseudotime(fucci, color_by = "gate")

#p1 <- ggplot(merge_fl_data,aes(y=mCherry.A,x=Venus.A,color=gate))+geom_point()+scale_y_log10()+scale_x_log10()


p1 <- ggplot(merge_fl_data,aes(x=pseudotime,y=mCherry.A,color=gate))+geom_point()+scale_y_log10()
p2 <- ggplot(merge_fl_data,aes(x=pseudotime,y=Venus.A,color=gate))+geom_point()+scale_y_log10()
p3 <- ggplot(merge_fl_data,aes(y=mCherry.A,x=Venus.A,color=gate))+geom_point()+scale_y_log10()+scale_x_log10()
#p4 <- ggplot(merge_fl_data,aes(y=AmCyan.A,x=pseudotime2,color=gate))+geom_point()+scale_y_log10()
gridExtra::grid.arrange(p1,p2, p3, nrow = 3)




fl_fucci2 <- fucci_fl_monocle_pseudotime(fucci2_fcsdata)
fucci2_fcsdata$pseudotime <- fl_fucci2$Pseudotime
#p1 <- ggplot(fucci2_fcsdata,aes(x=pseudotime,y=mCherry.A,color=gate))+geom_point()+scale_y_log10()
#p2 <- ggplot(fucci2_fcsdata,aes(x=pseudotime,y=Venus.A,color=gate))+geom_point()+scale_y_log10()
ggplot(ctl_fcsdata,aes(y=mCherry.A,x=Venus.A,color=gate))+geom_point()+scale_y_log10()+scale_x_log10()

ggplot(ctl_fcsdata,aes(x=factor(gate),y=pseudotime))+geom_violin()

p3 <- ggplot(fucci2_fcsdata,aes(y=mCherry.A,x=Venus.A,color=gate))+geom_point()+scale_y_log10()+scale_x_log10()
p4 <- ggplot(fucci2_fcsdata,aes(x=factor(gate),y=pseudotime))+geom_violin()
#p4 <- ggplot(merge_fl_data,aes(y=AmCyan.A,x=pseudotime2,color=gate))+geom_point()+scale_y_log10()
gridExtra::grid.arrange(p1,p2, p3, p4,nrow = 2)






ggplot(merge_fl_data,aes(x=factor(gate),y=pseudotime))+geom_violin()
gridExtra::grid.arrange(P1,P2, nrow = 2)
cellcycle_prog_mean_sd <- merge_fl_data %>% group_by(gate) %>% summarize(mean = mean(pseudotime), sd = sd(pseudotime))






merge_fl_data_melt <-
  melt(merge_fl_data, id.vars=c("time","gate"),
       measure.vars = c("pseudotime"))


ggplot(merge_fl_data_melt,aes(x=factor(time),y=value, color=gate))+ geom_quasirandom(size=0.5)
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