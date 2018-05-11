## Dot plot of change in class relative abundance between coastal and inland sites
```
#read in data
mimulus_class_delta <- read.delim("mim_classes_delta.txt")

library(ggplot2)
library(plyr)
library(reshape2)

#shape data and summarize, getting mean and se
mim_class_delta_m<-melt(mimulus_class_delta)
mim_delta_sum<-ddply(mim_class_delta_m, c("variable"), summarize, N=length(value), sd=sd(value), mean=mean(value), se=sd/sqrt(N))

#remove things other than names from class names
mim_delta_sum$variable2<-mim_delta_sum$variable
mim_delta_sum$variable<-gsub('k__Bacteria.', '', mim_delta_sum$variable)
mim_delta_sum$variable<-gsub('k__Archaea.', '', mim_delta_sum$variable)
mim_delta_sum$variable<-gsub('p__', '', mim_delta_sum$variable)
mim_delta_sum$variable<-gsub('c__', '', mim_delta_sum$variable)
mim_delta_sum$variable<-gsub('[.]', '-', mim_delta_sum$variable)
mim_delta_sum$variable<-gsub('--', '-', mim_delta_sum$variable)
mim_delta_sum$variable<-gsub('Cyanobacteria-', 'Cyanobacteria', mim_delta_sum$variable)

#make sure all names are fixed
unique(mim_delta_sum$variable)

#dot plot with mean/se
ggplot(mim_delta_sum, aes(variable, mean, colour=variable), show.legend=F)+
  geom_point(aes(cex=1.2), show.legend = F)+
  ylab("Delta Relative Abundance (Coastal-Inland)")+
  theme_bw()+
  coord_flip()+
  xlab("")+
  theme(text = element_text(size=20))+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se, width=0.4))+
  scale_colour_manual(values=c("Red", "Grey", "Red", "Red", "Red", "Red", "Red", "Red", "Grey", "grey", "Red", "Red", "Red", "Red", "Red", "black", "Red", "Red", "Red", "Red", "Red", "Red", 'Grey', 'grey', "Red", 'grey'))+
  theme(legend.position="none")
  ```
