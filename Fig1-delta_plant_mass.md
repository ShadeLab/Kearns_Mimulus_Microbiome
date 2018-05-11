## Box and whisker plot of change in plant mass between coastal and inland site
```
#read in data
library(readr)
delta_plant_mass <- read_delim("delta_plant_mass.txt",  "\t", escape_double = FALSE, trim_ws = TRUE)
#melt and plot data
library(reshape2)
delta_m<-melt(delta_plant_mass)
library(ggplot2)
#replance 'Delta'
delta_m$variable<-gsub('Delta_', '', delta_m$variable)
delta_m$variable<-gsub('shoots', 'Aboveground', delta_m$variable)
delta_m$variable<-gsub('Roots', 'Belowground', delta_m$variable)


ggplot(delta_m, aes(Genotype, value, fill=Origin))+
  theme_bw()+
  geom_boxplot()+
  facet_wrap(~variable, scales='free')+
  xlab("")+
  ylab("Change in Plant Mass (Inland-Coastal)")+
  coord_flip()
```
