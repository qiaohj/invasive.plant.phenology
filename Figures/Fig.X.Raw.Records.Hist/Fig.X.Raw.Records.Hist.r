library(data.table)
library(ggplot2)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")

dfff<-readRDS("../Figures/Fig.X.Raw.Records.Hist/Fig.X.Raw.Records.Hist.rda")
p<-ggplot(dfff)+geom_histogram(aes(x=month, fill=hemi), bins=50)+
  geom_vline(xintercept = 7, linetype=2, color="grey")+
  facet_grid(growthform~phenology.l)+
  labs(fill="Hemisphere", x="Month", y="")+
  theme_bw()
p
ggsave(p, filename="../Figures/Fig.X.Raw.Records.Hist/Fig.X.Raw.Records.Hist.png", width=7, height=6)
