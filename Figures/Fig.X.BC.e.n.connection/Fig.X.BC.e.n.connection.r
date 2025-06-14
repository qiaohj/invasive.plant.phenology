library(ggalluvial)
library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
paired.differ<-readRDS("../Data/paired.differ.rda")

alluv.df<-paired.differ[,.(N=.N), by=.(e.invasive, e.bot_country, n.scientific_name)]

allev.p<-ggplot(alluv.df, aes(axis1 = e.bot_country,
                     axis2 = e.invasive,
                     axis3 = n.scientific_name,
                     y = N)) +
  geom_alluvium(aes(fill = e.invasive), width = 1/12, alpha = 0.7) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Invasive Species", "Country", "Native Species"), expand = c(.05, .05)) +
  theme(legend.position = "none")
ggsave(allev.p, filename="../Figures/Fig.X.BC.e.n.connection/Fig.X.BC.e.n.connection.png",
       width=8, height=6)
