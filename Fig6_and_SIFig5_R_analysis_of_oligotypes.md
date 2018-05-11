# R Analysis of Mimulus Oligotypes (graph, PERMANOVA)

```
############Betaproteobacteria
setwd("./mimulus")
beta_oligo <- read.delim("~/beta_oligo.txt")
library(ggplot2)
library(reshape2)
beta_oligo_m<-melt(beta_oligo)

beta_colour<-rainbow(198, s=.6, v=.9)[sample(1:198,198)]

  ggplot(beta_oligo_m, aes(x=Sample, y=value, fill=variable))+
    geom_bar(stat='identity', colour='black')+
    theme_bw()+
    facet_wrap(~Site_planted, scales='free')+
    scale_y_continuous(expand=c(0,0))+
    ylab("Relative Abundance")+
    xlab("")+
      #coord_cartesian((ylim=c(0,100)))+
    scale_fill_manual(values=beta_colour)+
    theme(legend.position="none") +
    theme(text = element_text(size=20))+
    theme(axis.text.x = element_text(angle = 90))

#beta PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
beta_dist<-vegdist(beta_oligo[7:198], method='bray')
  
#by site
adonis(beta_dist ~ Site_planted, data=beta_oligo, permutations=100)
#F=5.58(1,37), p=0.0009, R2=0.41

#by genotype/site
adonis(beta_dist ~ Genotype, data=beta_oligo, permutations=100)
#F=0.62(3,35), p=0.8212, R2=0.05

#by genotype/site
adonis(beta_dist ~ Genotype+Site_planted, data=beta_oligo, permutations=100)
#F=1.09(4,34), p=0.8212, R2=0.05


############Actinobacteria 
actino_oligo <- read.delim("~/actino_oligo.txt")
library(ggplot2)
library(reshape2)
actino_oligo_m<-melt(actino_oligo)
  
actino_colour<-rainbow(198, s=.6, v=.9)[sample(1:198,198)]
  
ggplot(actino_oligo_m, aes(x=Sample, y=value, fill=variable))+
    geom_bar(stat='identity', colour='black')+
    theme_bw()+
    facet_wrap(~Site_planted, scales='free')+
    scale_y_continuous(expand=c(0,0))+
    ylab("Relative Abundance")+
    xlab("")+
    #coord_cartesian((ylim=c(0,100)))+
    scale_fill_manual(values=actino_colour)+
    theme(legend.position="none") +
    theme(text = element_text(size=20))+
    theme(axis.text.x = element_text(angle = 90))


#Actino PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
actino_dist<-vegdist(actino_oligo[7:160], method='bray')

#by site
adonis(actino_dist ~ Site_planted, data=actino_oligo, permutations=100)
#F=12.556(1,37), p=0.0009, R2=0.25

#by genotype/site
adonis(actino_dist ~ Genotype, data=actino_oligo, permutations=100)
#F=0.75(3,35), p=0.8218, R2=0.06

#by genotype/site
adonis(actino_dist ~ Genotype+Site_planted, data=actino_oligo, permutations=100)
#F=0.99(4,34), p=0.445, R2=0.06

  
##############Acidobacteria
acido_oligo <- read.delim("~/acido_oligo.txt")
library(ggplot2)
library(reshape2)
acido_oligo_m<-melt(acido_oligo)

acido_colour<-rainbow(398, s=.6, v=.9)[sample(1:398,398)]

ggplot(acido_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=acido_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))

#Acido PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
acido_dist<-vegdist(acido_oligo[7:160], method='bray')
#acido_mds<-metaMDS(acido_dist)

#adonis by site
adonis(acido_dist ~ Site_planted, data=acido_oligo, permutations=10000)
#F=17.986(1,38), p=9.9e-5, R2=0.32

#adonis by genotype/site
adonis(acido_dist ~ Genotype, data=acido_oligo, permutations=10000)
#F=0.89(3,35), p=0.56, R2=0.06

#adonis by genotype/site
adonis(acido_dist ~ Genotype+Site_planted, data=acido_oligo, permutations=10000)
#F=0.99(4,34), p=0.156, R2=0.06

################Alphaproteobacteria
alpha_oligo <- read.delim("~/alpha_oligo.txt")
library(ggplot2)
library(reshape2)
alpha_oligo_m<-melt(alpha_oligo)

alpha_colour<-rainbow(160, s=.6, v=.9)[sample(1:160,160)]

ggplot(alpha_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=alpha_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))

#Alpha PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
alpha_dist<-vegdist(alpha_oligo[7:166], method='bray')

#adonis by site
adonis(alpha_dist ~ Site_planted, data=alpha_oligo, permutations=10000)
#F=47.27(1,38), p=9.9e-5, R2=0.55

#adonis by genotype/site
adonis(alpha_dist ~ Genotype, data=alpha_oligo, permutations=10000)
#F=0.60(3,35), p=0.74, R2=0.04

#adonis by genotype/site
adonis(alpha_dist ~ Genotype+Site_planted, data=alpha_oligo, permutations=10000)
#F=1.39(4,34), p=0.199, R2=0.04

##################Cyanobacteria
cyano_oligo <- read.delim("~/cyano_oligo.txt")
library(ggplot2)
library(reshape2)
cyano_oligo_m<-melt(cyano_oligo)

cyano_colour<-rainbow(160, s=.6, v=.9)[sample(1:160,160)]

ggplot(cyano_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=cyano_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))

#Cyano PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
cyano_dist<-vegdist(cyano_oligo[7:139], method='bray')

#adonis by site
adonis(cyano_dist ~ Site_planted, data=cyano_oligo, permutations=10000)
#F=14.62(1,38), p=9.9e-5, R2=0.30

#adonis by genotype/site
adonis(cyano_dist ~ Genotype, data=cyano_oligo, permutations=10000)
#F=0.81(3,35), p=0.61, R2=0.07

#adonis by genotype/site
adonis(cyano_dist ~ Genotype+Site_planted, data=cyano_oligo, permutations=10000)
#F=1.18(4,34), p=0.28, R2=0.07

#################Verrucumicrobia
verruco_oligo <- read.delim("~/verruco_oligo.txt")
library(ggplot2)
library(reshape2)
verruco_oligo_m<-melt(verruco_oligo)

verruco_colour<-rainbow(319, s=.6, v=.9)[sample(1:319,319)]

ggplot(verruco_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=verruco_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))

#Verruco PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
verruco_dist<-vegdist(verruco_oligo[7:325], method='bray')

#adonis by site
adonis(verruco_dist ~ Site_planted, data=verruco_oligo, permutations=10000)
#F=28.61(1,38), p=9.9e-5, R2=0.42

#adonis by genotype/site
adonis(verruco_dist ~ Genotype, data=verruco_oligo, permutations=10000)
#F=0.65(3,35), p=0.61, R2=0.05

#adonis by genotype/site
adonis(verruco_dist ~ Genotype+Site_planted, data=verruco_oligo, permutations=10000)
#F=1.15(4,34), p=0.28, R2=0.05


####################Chloroflexi
chloro_oligo <- read.delim("~/chloroflexi.txt")
library(ggplot2)
library(reshape2)
chloro_oligo_m<-melt(chloro_oligo)

chloro_colour<-rainbow(319, s=.6, v=.9)[sample(1:319,319)]

ggplot(chloro_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=chloro_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))


#Chloroflexi PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
chloro_dist<-vegdist(chloro_oligo[7:190], method='bray')

#adonis by site
adonis(chloro_dist ~ Site_planted, data=chloro_oligo, permutations=10000)
#F=19.27(1,38), p=9.9e-5, R2=0.333

#adonis by genotype/site
adonis(chloro_dist ~ Genotype, data=chloro_oligo, permutations=10000)
#F=1.07(3,35), p=0.36, R2=0.08

#adonis by genotype/site
adonis(verruco_dist ~ Genotype+Site_planted, data=verruco_oligo, permutations=10000)
#F=1.16(4,34), p=0.284, R2=0.051


#############Gammaproteobacteria
gamma_oligo <- read.delim("~/gamma_oligo.txt")
library(ggplot2)
library(reshape2)
gamma_oligo_m<-melt(gamma_oligo)

gamma_colour<-rainbow(319, s=.6, v=.9)[sample(1:319,319)]

ggplot(gamma_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=gamma_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))


#gamma PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
gamma_dist<-vegdist(gamma_oligo[7:146], method='bray')

#adonis by site
adonis(gamma_dist ~ Site_planted, data=gamma_oligo, permutations=10000)
#F=18.818(1,38), p=0.9.99e-5, R2=0.33

#adonis by genotype/site
adonis(gamma_dist ~ Genotype, data=gamma_oligo, permutations=10000)
#F=90(3,35), p=0.54, R2=0.06

#adonis by genotype/site
adonis(gamma_dist ~ Genotype+Site_planted, data=gamma_oligo, permutations=10000)
#F=1.35(4,34), p=0.155, R2=0.06


#####################################################################################################################
##########################################Non-signifcant below#######################################################
#####################################################################################################################

#############Bacillus 
bacillus_oligo <- read.delim("~/bacillus_oligo.txt")
library(ggplot2)
library(reshape2)
bacillus_oligo_m<-melt(bacillus_oligo)

bacillus_colour<-rainbow(319, s=.6, v=.9)[sample(1:319,319)]

ggplot(bacillus_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=bacillus_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))


#Bacillus PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
bacillus_dist<-vegdist(bacillus_oligo[7:20], method='bray')

#adonis by site
adonis(bacillus_dist ~ Site_planted, data=bacillus_oligo, permutations=10000)
#F=0.93(1,38), p=0.698, R2=0.02

#adonis by genotype/site
adonis(bacillus_dist ~ Genotype, data=bacillus_oligo, permutations=10000)
#F=1.01(3,35), p=0.446, R2=0.07

#adonis by genotype/site
adonis(bacillus_dist ~ Genotype+Site_planted, data=bacillus_oligo, permutations=10000)
#F=1.16(4,34), p=0.284, R2=0.07

#############Bacteroidetes
bact_oligo <- read.delim("~/bact_oligo.txt")
library(ggplot2)
library(reshape2)
bact_oligo_m<-melt(bact_oligo)

bact_colour<-rainbow(319, s=.6, v=.9)[sample(1:319,319)]

ggplot(bact_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=bact_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))


#Bacteroidetes PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
bact_dist<-vegdist(bact_oligo[7:129], method='bray')

#adonis by site
adonis(bact_dist ~ Site_planted, data=bact_oligo, permutations=10000)
#F=1.66(1,38), p=0.09, R2=0.04

#adonis by genotype/site
adonis(bact_dist ~ Genotype, data=bact_oligo, permutations=10000)
#F=1.38(3,35), p=0.10, R2=0.10

#adonis by genotype/site
adonis(bact_dist ~ Genotype+Site_planted, data=bact_oligo, permutations=10000)
#F=1.41(4,34), p=0.09, R2=0.10


#############Archaea
archa_oligo <- read.delim("~/arch_oligo.txt")
library(ggplot2)
library(reshape2)
archa_oligo_m<-melt(archa_oligo)

archa_colour<-rainbow(319, s=.6, v=.9)[sample(1:319,319)]

ggplot(archa_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=archa_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))


#Archaea PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
archa_dist<-vegdist(archa_oligo[7:10], method='bray')

#adonis by site
adonis(archa_dist ~ Site_planted, data=archa_oligo, permutations=10000)
#F=0.40(1,38), p=0.71, R2=0.01

#adonis by genotype/site
adonis(archa_dist ~ Genotype, data=archa_oligo, permutations=10000)
#F=0.77(3,35), p=0.61, R2=0.06

#adonis by genotype/site
adonis(archa_dist ~ Genotype+Site_planted, data=archa_oligo, permutations=10000)
#F=.76(4,34), p=0.61, R2=0.06


#############Deltaproteobacteria
delta_oligo <- read.delim("~/delta_oligo.txt")
library(ggplot2)
library(reshape2)
delta_oligo_m<-melt(delta_oligo)

delta_colour<-rainbow(319, s=.6, v=.9)[sample(1:319,319)]

ggplot(delta_oligo_m, aes(x=Sample, y=value, fill=variable))+
  geom_bar(stat='identity', colour='black')+
  theme_bw()+
  facet_wrap(~Site_planted, scales='free')+
  scale_y_continuous(expand=c(0,0))+
  ylab("Relative Abundance")+
  xlab("")+
  #coord_cartesian((ylim=c(0,100)))+
  scale_fill_manual(values=delta_colour)+
  theme(legend.position="none") +
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(angle = 90))


#delta PERMANOVA
library(vegan)
#calculate Bray-Curtis Distance
delta_dist<-vegdist(delta_oligo[7:195], method='bray')

#adonis by site
adonis(delta_dist ~ Site_planted, data=delta_oligo, permutations=10000)
#F=0.93(1,38), p=0.5, R2=0.10

#adonis by genotype/site
adonis(delta_dist ~ Genotype, data=delta_oligo, permutations=10000)
#F=1.30(3,35), p=0.14, R2=0.10

#adonis by genotype/site
adonis(delta_dist ~ Genotype+Site_planted, data=delta_oligo, permutations=10000)
#F=1.41(4,34), p=0.09, R2=0.10
```
