## Stacked bar plot of bacterial classes
```
#read data
mimulus_classes <- read.delim("~/mimulus_classes.txt")

#melt data
library(reshape2)
mimulus_classes_m<-melt(mimulus_classes)

#26 colour palett
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

#plot it
library(ggplot2)
ggplot(mimulus_classes_m, aes(x=Order, y=value, fill=variable))+
  geom_bar(stat='identity')+
  coord_cartesian(ylim=c(0,1))+
  theme_bw()+
  scale_fill_manual(values=pal)+
  scale_y_continuous(expand = c(0, 0))+
  theme(text = element_text(size=20))+
  ylab("Relative Abundance")+
  xlab("")+
  guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x = element_text(angle = 90))
  ```
