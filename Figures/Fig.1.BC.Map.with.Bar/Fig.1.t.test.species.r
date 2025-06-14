library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")


all.species.df<-readRDS("../Figures/Fig.1.BC.Map.with.Bar/Fig.1.t.test.species.rda")


t.test.species<-ggplot(all.species.df)+
  geom_vline(aes(xintercept=0), linetype=2, color="black", alpha=0.5)+
  geom_errorbarh(aes(xmin=ci.low, xmax=ci.high, y=species, color=growthform, linetype=lt),
                 height=0.5)+
  geom_point(aes(x=mean.difference, y=species, color=growthform))+
  scale_y_discrete (position = "right", expand = expansion(mult = c(0.01, 0.01)))+
  theme_bw()+
  theme(legend.position="none",
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face = "italic", size=12),
        title = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size=15),
        strip.background = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5, 0, 5, 0))+
  facet_wrap(~phenology)
t.test.species
ggsave(t.test.species, filename="../Figures/Fig.1.BC.Map.with.Bar/Fig.1.t.test.species.png",
       width=7, height=6)
