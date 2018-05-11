## PCoA
```
library(readr)
library(vegan)
#otu table
mim_otu<-read.table("mim_otu_table.txt", header=T, row.names=1)
#environmental data
mim_envir<-read.table("mim_envir.txt", header=T, row.names=1)
#class data
mim_class<-read.table("mim_class.txt", header=T, row.names=1)
#pcoa
mim.pcoa<-capscale(mim_otu, distance='bray')


#fit environment
mim.envfit<-envfit(mim.pcoa, mim_envir)
#fit class taxonomy
mim.classfit<-envfit(mim.dis, mim_class)

#environment scores
mim.scores<-as.data.frame(scores(mim.envfit, display='vectors'))

#wite scores to file
write.csv(mim.scores, 'mim.score.csv')

#re-read in scores with environ data
mim.scores<-read.table("mim.scores.txt", header=T, row.names=1, sep="")

#remove p>0.1
mim.scores<-mim.scores[-c(2,3,9,11,12,14,15),]

#read in PCoA coords
library(readr)
mimulus_coords <- read_delim("mimulus_coords.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

#plot PCoA with envfit
library(ggplot2)
ggplot(mimulus_coords, aes(WU_Axis1, WU_axis2))+
    theme_bw()+
  geom_point(aes(cex=1.2, shape=Site_planted, colour=Genotype))+
  xlab("PC1- 49.21%")+
  ylab("PC2- 15.21%")+
  geom_segment(data=mim.scores, aes(x=0, xend=MDS1, y=0, yend=MDS2), arrow = arrow(length = unit(0.1, "cm")), colour = "black")+
  geom_text(data = mim.scores, aes(x = MDS1, y = MDS2, label = variable), size = 6)+
  theme(text = element_text(size=20))
```
