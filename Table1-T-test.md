## Change in geochemistry between sites
```
library(ggplot2)
library(reshape2)
library(readr)
mimulus_envir <- read_delim("mimulus_envir.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

mimulus_m<-melt(mimulus_envir)
str(mimulus_m)

ggplot(mimulus_m, aes(Site, value, fill=Site))+
  geom_boxplot()+
  facet_wrap(~variable, scales="free")+
  theme_bw()+
  ylab("Concentration")+
  xlab("Site")

#summary stats for table
library(plyr)
mim_geo_sum<-ddply(mimulus_m, c("Site", "variable"), summarise, N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
write.csv(mim_geo_sum, "mim_geo_sum.csv")  

#maths
t.test(pH ~ Site, data=mimulus_envir)
t.test(Phosphorous ~ Site, data=mimulus_envir)
t.test(Potassium ~ Site, data=mimulus_envir)
t.test(Calcium ~ Site, data=mimulus_envir)
t.test(Magnesium ~ Site, data=mimulus_envir)
t.test(Copper ~ Site, data=mimulus_envir)
t.test(Percent_Organic_Matter ~ Site, data=mimulus_envir)
t.test(Sodium ~ Site, data=mimulus_envir)
t.test(Nitrate ~ Site, data=mimulus_envir)
t.test(Ammonium ~ Site, data=mimulus_envir)
t.test(Percent_Moisture ~ Site, data=mimulus_envir)
t.test(`Total-N` ~ Site, data=mimulus_envir)
t.test(Sulfur ~ Site, data=mimulus_envir)
t.test(Soil_tempC ~ Site, data=mimulus_envir)
t.test(Air_TempC ~ Site, data=mimulus_envir)
#most things significantly differ between sites
```
