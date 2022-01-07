datadir <-"/home/samba/public/shintaku/20211026HiSeqX005_hunter/downsample/"
read_summary1 <- read.table(file.path(datadir,"read_summary.txt"),header = TRUE)
datadir <-"/home/samba/public/shintaku/20211124HiSeqX006_hunter//downsample/"
read_summary2 <- read.table(file.path(datadir,"read_summary.txt"),header = TRUE)
read_summary <- rbind(read_summary1,read_summary2)
read_summary_ggplot <- reshape2::melt(read_summary,id.vars = c("plate","sample","frac_filter","frac_align","frac_duplicate"))
ggplot(read_summary_ggplot,aes(x=factor(variable),y=value,color=plate))+geom_jitter()+scale_y_log10()
# fraction
read_summary$frac_filter <- read_summary$filtered_reads/read_summary$raw_reads
read_summary$frac_align <- read_summary$aligned_reads/read_summary$filtered_reads
read_summary$frac_duplicate <- read_summary$UMI_count/read_summary$aligned_reads
read_summary$plate <- substr(read_summary$sample,1,3)
read_summary_ggplot <- reshape2::melt(read_summary,id.vars = c("sample","raw_reads","aligned_reads","filtered_reads","UMI_count","plate"))
ggplot(read_summary_ggplot,aes(x=factor(variable),y=value,color=sample))+geom_jitter()#+scale_y_log10()
#
read_summary$plate <- substr(read_summary$sample,1,3)
rownames(read_summary) <- read_summary$sample
read_summary[read_summary$frac_align<0.4,]$raw_reads
read_summary_ggplot <- reshape2::melt(read_summary,id.vars = 
                                        c("sample","aligned_reads","filtered_reads","UMI_count","plate"))
ggplot(read_summary_ggplot,aes(x=factor(variable),y=value,color=plate))+
  geom_jitter()#+scale_y_log10()

read_summary[read_summary$frac_align>0.4,] %>%
  dplyr::group_by(plate) %>%
  dplyr::summarise(raw_reads.mean=mean(raw_reads))
