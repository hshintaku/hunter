fucci_fcs_io <- function (fcs_dir,exp,phase){
  #exp <- "Fucci32"
  fcs_file <- list.files(fcs_dir,pattern=".fcs$")
  gate_file <- list.files(fcs_dir,pattern=".fcs_gates.xml")
  
  flowEnv <- new.env()
  fcFile <-file.path(fcs_dir,fcs_file[1])
  gate_dir <- file.path(fcs_dir,gate_file[1])
  
  
  fcs <- read.FCS(fcFile)
  
  read.gatingML(gate_dir,flowEnv)
  parameters<- ls(flowEnv)
  ctl_fcsdata <- data.frame(fcs@exprs)
  #colnames(ctl_fcsdata) <- c("time","ImageNumber","ObjectNumber","R","G")
  ctl_fcsdata$exp <- exp
  ctl_fcsdata$gate <- phase[1]
  filter <- flowCore::filter
  nphase <- length(phase)-1
  for (icnt in 1:nphase){
    result = filter(fcs, flowEnv[[parameters[icnt+1]]])
    summary(result)
    g1_index <- result[[paste0(parameters[icnt+1],"+")]]@subSet
    #g1_data <- fcs@exprs[g1_index,]
    ctl_fcsdata[g1_index,]$gate <- phase[icnt+1]
  }
  filter <- dplyr::filter
  return(ctl_fcsdata)
} 
