## Shannon Diversity and phylogentic diversity box-and-whisker plot
```
mimulus_coords<-read.table("mimulus_coords.txt", header=T)
mimulus_coords$PD_whole_tree_alpha<-as.numeric(mimulus_coords$PD_whole_tree_alpha)
mimulus_coords$shannon_alpha<-as.numeric(mimulus_coords$shannon_alpha)

#purge the PCoA axis values
mimulus_diversity<-mimulus_coords[,c(-2,-3)]
#reshape the data
library(reshape2)
mimulus_diversity_m<-melt(mimulus_diversity)
#plot ze data
library(ggplot2)
ggplot(mimulus_diversity_m, aes(Genotype, value, fill=Origin))+
  facet_grid(variable ~ Site_planted, scales='free')+
  theme_bw()+
  ylab("Bacterial Diversity")+
  geom_boxplot()+
  theme(axis.text=element_text(size=15, colour="black"), axis.title=element_text(size=20, colour="black"), strip.text = element_text(size = 12))+
  xlab("")

###maths
#test variance
bartlett.test(shannon_alpha ~ GenotypeSite, data= mimulus_diversity)
#p=0.1139
library(car)
leveneTest(shannon_alpha ~ GenotypeSite, data= mimulus_diversity)
#p=0.68
bartlett.test(PD_whole_tree_alpha ~ GenotypeSite, data= mimulus_diversity)
#p=0.27
leveneTest(PD_whole_tree_alpha ~ GenotypeSite, data= mimulus_diversity)
#0.31

#will use a ANOVA w/Tukey to assess differences in diversity

#######Effects of site
#GenotypeSite is a genotype+site it is planted
shan.aov.site<-aov(shannon_alpha ~Site_planted+GenotypeSite, data= mimulus_diversity)
summary(shan.aov.site)
#F=7.61, p(site)=0.0021, p(GT)=0.00211
TukeyHSD(shan.aov.site)


pd.aov.site<-aov(PD_whole_tree_alpha ~Site_planted+GenotypeSite, data= mimulus_diversity)
summary(pd.aov.site)
#F=21, p(site)=0.0088, p(GT)=0.0000344
pd.tukey<-TukeyHSD(pd.aov.site)
View(pd.tukey$GenotypeSite)
```
