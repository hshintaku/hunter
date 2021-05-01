#
# microscopy data 
#
phase <- c("MG1","S","G1","lateSG2M")
fcs_dir <- "/home/samba/pihome/2021/Kaneko/1_Data/2_Microscopy/20210316/20210421_Output/No1_Control"
ctl_fcsdata <- fucci_fcs_io(fcs_dir,"CTL1",phase)
ctl_cellnames <- paste0("CTL-HeLaFucci-",ctl_fcsdata$image,"-",ctl_fcsdata$object)
rownames(ctl_fcsdata) <- ctl_cellnames
ggplot(ctl_fcsdata,aes(x=Venus.A,y=mCherry.A,color=gate))+geom_point()

phase <- c("MG1","lateSG2M","S","G1")
fcs_dir <- "/home/samba/pihome/2021/Kaneko/1_Data/2_Microscopy/20210316/20210421_Output/No2_ELP/"
elp_fcsdata <- fucci_fcs_io(fcs_dir,"ELP1",phase)
elp_cellnames <- paste0("ELP-HeLaFucci-",elp_fcsdata$image,"-",elp_fcsdata$object)
rownames(elp_fcsdata) <- elp_cellnames
ggplot(elp_fcsdata,aes(x=Venus.A,y=mCherry.A,color=gate))+geom_point()


phase <- c("MG1","S","G1","lateSG2M")
fcs_dir <- "/home/samba/pihome/2021/Kaneko/1_Data/2_Microscopy/20210316/20210421_Output/No3/"
elp3_fcsdata <- fucci_fcs_io(fcs_dir,"ELP2",phase)
elp3_cellnames <- paste0("ELP-HeLaFucci-",elp3_fcsdata$image,"-",elp3_fcsdata$object)
rownames(elp3_fcsdata) <- elp3_cellnames
ggplot(elp3_fcsdata,aes(x=Venus.A,y=mCherry.A,color=gate))+geom_point()


phase <- c("MG1","lateSG2M","S","G1")
fcs_dir <- "/home/samba/pihome/2021/Kaneko/1_Data/2_Microscopy/20210316/20210421_Output/No4/"
elp4_fcsdata <- fucci_fcs_io(fcs_dir,"CTL2",phase)
elp4_cellnames <- paste0("CTL-HeLaFucci-",elp4_fcsdata$image,"-",elp4_fcsdata$object)
rownames(elp3_fcsdata) <- elp4_cellnames
ggplot(elp4_fcsdata,aes(x=Venus.A,y=mCherry.A,color=gate))+geom_point()


phase <- c("MG1","S","G1","lateSG2M")
fcs_dir <- "/home/samba/pihome/2021/Kaneko/1_Data/2_Microscopy/20210316/20210421_Output/No5/"
elp5_fcsdata <- fucci_fcs_io(fcs_dir,"ELP4",phase)
elp5_cellnames <- paste0("ELP-HeLaFucci-",elp5_fcsdata$image,"-",elp5_fcsdata$object)
rownames(elp5_fcsdata) <- elp5_cellnames
ggplot(elp5_fcsdata,aes(x=Venus.A,y=mCherry.A,color=gate))+geom_point()



phase <- c("MG1","S","lateSG2M","G1")
fcs_dir <- "/home/samba/pihome/2021/Kaneko/1_Data/2_Microscopy/20210316/20210421_Output/No6/"
elp6_fcsdata <- fucci_fcs_io(fcs_dir,"ELP5",phase)
elp6_cellnames <- paste0("ELP-HeLaFucci-",elp6_fcsdata$image,"-",elp6_fcsdata$object)
rownames(elp6_fcsdata) <- elp6_cellnames
ggplot(elp6_fcsdata,aes(x=Venus.A,y=mCherry.A,color=gate))+geom_point()



merge_fl_data <- rbind(ctl_fcsdata,elp_fcsdata,elp3_fcsdata,elp4_fcsdata,elp5_fcsdata,elp6_fcsdata)


merge_fl_data$condition <- "ELP"
merge_fl_data[merge_fl_data$exp %in% c("CTL1","CTL2"),]$condition <- "CTL"
