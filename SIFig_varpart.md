## Partition variance due to environmental variables and make stacked barplot of variance explained
```
library(ggplot2)
library(vegan)

#read in rarefied OTU table
mim_rare<-read.table("single_rare.txt", header=T, row.names=1)

#transpose data
mim_rare_t<-t(mim_rare)

#write to .csv and reload the data
write.csv(mim_rare_t, 'mim_rare_t.csv')
mim_rare_t<-read.csv('mim_rare_t.csv', header=T, row.names = 1)

#read environmental variables
mim_envir<-read.table('mim_envir.txt', header=T, row.names = 1)

#partition variance
mim.varpart1<-varpart(mim_rare_t, ~Genotype, ~Origin, ~Site_planted,  data=mim_envir)
mim.varpart2<-varpart(mim_rare_t, Nitrate + TotalN +Ammonium, data=mim_envir)
mim.varpart3<-varpart(mim_rare_t, pH + Phosphorous + Potassium + Calcium,  data=mim_envir)
mim.varpart4<-varpart(mim_rare_t, Magnesium + Copper + Percent_Organic_Matter + Sodium,  data=mim_envir)
mim.varpart5<-varpart(mim_rare_t, Percent_Moisture + Sulfur + Soil_tempC + Air_TempC,  data=mim_envir)
```

## Plot data
```
mim_var_tot<-read.table("mim_tot_var.txt", header=T)
library(reshape2)
mim_var_m<-melt(mim_var_tot)
tot_var<-ggplot(mim_var_m, aes(variable, value, fill=Variable, colour=variable))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_colour_manual(values=c('black', 'black', 'black'))+
  scale_fill_manual(values=c("black", "white", "grey"))+
  xlab("")+
  scale_y_continuous(expand=c(0,0))+
  ylab("Percent Variation Explained")+
  ggtitle("A) Total Variance")+
  theme(axis.text=element_text(size=15, colour="black"), axis.title=element_text(size=20, colour="black"), strip.text = element_text(size = 12))+
  coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

#read in table of site-based variance, parsed as % of variation explained by Site
mim_site_var<-read.table("mim_site_var.txt", header=T)
mim_site_m<-melt(mim_site_var)
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")
sitey_var<-ggplot(mim_site_m, aes(variable, value, fill=Variable, colour=variable))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_colour_manual(values=c('black', 'black', 'black'))+
  xlab("")+
  ylab("Percent Variation Explained")+
  scale_fill_manual(values=pal)+
  scale_y_continuous(expand=c(0,0))+
  ggtitle("B) Site Variance")+
  theme(axis.text=element_text(size=15, colour="black"), axis.title=element_text(size=20, colour="black"), strip.text = element_text(size = 12))+
  coord_cartesian(ylim=c(0,1))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

#read in table of plant-based variance, parsed as % of variation explained by Genotype
mim_plant_var<-read.table("mim_plant_var.txt", header=T)
plant_var_m<-melt(mim_plant_var)
planty_var<-ggplot(plant_var_m, aes(variable, value, fill=Var, colour=variable))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_colour_manual(values=c('black', 'black', 'black'))+
  xlab("")+
  scale_y_continuous(expand=c(0,0))+
  ylab("Percent Variation Explained")+
  coord_cartesian(ylim=c(0,1))+
  ggtitle("C) Plant Variance")+
  theme(axis.text=element_text(size=15, colour="black"), axis.title=element_text(size=20, colour="black"), strip.text = element_text(size = 12))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

library(gridExtra)
grid.arrange(tot_var, sitey_var, planty_var, ncol=3, widths=c(1,1.2,1))
```
