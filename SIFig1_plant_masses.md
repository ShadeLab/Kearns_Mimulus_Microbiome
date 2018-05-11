## Plant mass box-and-whisker and t-tests
```
#read in data
full_mimulus_data <- read.delim("~/full_mimulus_data.txt")

#remove non-plant things
mimulus_plants<-full_mimulus_data[,c(-2,-3,-9,-10,-11,-12,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,-27)]

#melt and plot data
library(reshape2)
library(ggplot2)
mimulus_plants_m<-melt(mimulus_plants)
ggplot(mimulus_plants_m, aes(Genotype, value, fill=Origin))+
  geom_boxplot()+
  facet_grid(variable ~ Site_planted, scales='free')+
  theme_bw()+
  ylab("Plant Mass")+
  theme(axis.text=element_text(size=15, colour="black"), axis.title=element_text(size=20, colour="black"), strip.text = element_text(size = 12))+
  xlab("")

##Variance test
bartlett.test(Root_Mass_mg ~ GenotypeSite, data=mimulus_plants)  
#p=0.004
bartlett.test(Shoot_Mass_g ~ GenotypeSite, data=mimulus_plants)  
#p=1.116e-12
bartlett.test(Root_Mass_mg ~ Site_planted, data=mimulus_plants)
#p=0.003
bartlett.test(Shoot_Mass_g ~ Site_planted, data=mimulus_plants)
#p=1.36e-12
#will test differences in plant performance with a Welch's t-test :(

#roots between sites
t.test(mimulus_plants$Root_Mass_mg ~ mimulus_plants$Site_planted)
#t=4.71, p=6.15e-5, df=27.752

#shoots between sites
t.test(mimulus_plants$Shoot_Mass_g ~ mimulus_plants$Site_planted)
#t=2.75, df=19.652, p=0.01
#biomass is signifcantly lower at inland site, compared to coastal site

#all combos of pairwise t-tests, w/BH correction
root_t.tests<-pairwise.t.test(mimulus_plants$Root_Mass_mg, mimulus_plants$GenotypeSite, p.adjust.method = 'hochberg')
View(root_t.tests$p.value)

shoot_t.tests<-pairwise.t.test(mimulus_plants$Shoot_Mass_g, mimulus_plants$GenotypeSite, p.adjust.method = 'hochberg')
View(shoot_t.tests$p.value)
```
