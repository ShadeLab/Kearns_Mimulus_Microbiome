## OTU cumulative abundance
```
#read data
mimulus_cumul_abun <- read.delim("~/mimulus_cumul_abun.txt")
library(ggplot2)

#determine how many are more abundant than 30 reads/sample and in more than 25 samples
count25<-which(mimulus_cumul_abun$No_Samples<5 & mimulus_cumul_abun$Log_Mean_Abun<1.4)
str(count25)

#extract OTUs of interest (>30 samples, >30 reads/sample)
rows_mim<-which(mimulus_cumul_abun$No_Samples>30 & mimulus_cumul_abun$Log_Mean_Abun>1.5)
rows_mim2<-mimulus_cumul_abun[rows_mim,]

#plot it, colouring OTUs of interest in black and all other points grey
ggplot(mimulus_cumul_abun, aes(Log_Mean_Abun, No_Samples))+
  geom_point(colour='grey')+
  geom_point(data=rows_mim2, aes(Log_Mean_Abun, No_Samples), colour='black')+
  theme_bw()+
  ylab("Number of Samples")+
  xlab("Log10 Mean Abundance")+
  theme(axis.text=element_text(size=15, colour="black"), axis.title=element_text(size=20, colour="black"), strip.text = element_text(size = 12))


#extract and write 'core' OTUs to table
mim_cuml_core<-subset(mimulus_cumul_abun, Log_Mean_Abun>1.5 & No_Samples>30)
write.csv(mim_cuml_core, "mim_core.csv")
```
