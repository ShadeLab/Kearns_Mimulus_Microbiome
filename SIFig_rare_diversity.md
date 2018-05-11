## Box-and-whisker plot of breakdown of phylogenetic diversity
```
#read in diversity for numerically abundant taxa
top_diversity<-read.table("top_taxa_map.txt", header=T)
#plot data
library(ggplot2)
top.pd<-ggplot(top_diversity, aes(Genotype, PD_whole_tree_alpha, fill=Origin))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Site_planted)+
  ylab("Phylogenetic Distance")+
  xlab("")+
  theme(text = element_text(size=20))+
  ggtitle("(A) Most Abundant Taxa (>1% rel. abun.)")

#read in diversity for numerically rare taxa
rare_diversity<-read.table("rare_taxa_map.txt", header=T)

#plot data
rare.pd<-  ggplot(rare_diversity, aes(Genotype, as.numeric(PD_whole_tree_alpha), fill=Origin))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Site_planted)+
  ylab("Phylogenetic Distance")+
  xlab("")+
    theme(text = element_text(size=20))+
  ggtitle("(B) Rare Taxa (<1% rel. abun.)")

#plot both together on same plot
library(gridExtra)
grid.arrange(top.pd, rare.pd, ncol=1)

#stats for top taxa
bartlett.test(top_diversity$PD_whole_tree_alpha, top_diversity$Site_planted)
#p=0.56
bartlett.test(top_diversity$PD_whole_tree_alpha, top_diversity$GenotypeSite)
#0.551
top.aov<-aov(PD_whole_tree_alpha ~ Site_planted+GenotypeSite, data=top_diversity)
summary(top.aov)
#p=0.97 (site), 0.71 (gt)
#no signifcant effect of GT or site

#stats for rare taxa
bartlett.test(rare_diversity$PD_whole_tree_alpha, rare_diversity$Site_planted)
#p=0.15
bartlett.test(rare_diversity$PD_whole_tree_alpha, rare_diversity$GenotypeSite)
#p=0.56
rare.aov<-aov(PD_whole_tree_alpha ~ Site_planted+GenotypeSite, data=rare_diversity)
summary(rare.aov)
#p=0.02, F=2.933
TukeyHSD(rare.aov)
```
