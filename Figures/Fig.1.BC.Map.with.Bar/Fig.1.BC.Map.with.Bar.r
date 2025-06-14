library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
dist<-10
bot_country <- st_read("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/bot_country_robin.shp") ## Transformed boundary of botanical country (Robinson projection)
bbox <- st_read("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/worldmap_boundary_box.shp") ## World map boundary
bbox_robin <- st_transform(bbox, crs="+proj=robin")
all.df<-readRDS("../Figures/Fig.1.BC.Map.with.Bar/Fig.1.BC.Map.with.Bar.A.rda")

bot_country_v<-merge(bot_country, all.df, by.x="LEVEL3_COD", by.y="bot_country")

p<-ggplot()+
  geom_sf(data=bbox_robin, fill = "#E0FFFF", alpha=0.3)+
  geom_sf(data=bot_country, fill="#F7F7F7", color="#A0A0A0", alpha=0.3)+
  geom_sf(data=bot_country_v, aes(fill=mean.difference))+
  coord_sf(expand = FALSE)+
  scale_fill_gradient2(midpoint = 0,
                       low = "#FC8D59",  # Red for negative difference
                       mid = "#FFFFBF",  # Yellow for neutral (0 difference)
                       high = "#91CF60",  # Green for positive difference
                       space = "Lab",
                       name = "Phenological difference")+
  facet_wrap(~phenology, nrow=2, strip.position = "left")+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.position = "bottom",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 15),
        plot.margin = margin(5, 0, 5, 0))#top, right, bottom, left

p


bar <- ggplot(all.df, aes(x = phenology)) +
  geom_bar(aes(fill = direction.factor.sc), 
           width = 0.5, color = "black", position = "stack") +  
  geom_text(stat = "count", 
            aes(group=direction.factor.sc, label = after_stat(count)), 
            position = position_stack(vjust = 0.5), color = "black", size = 5)+
  scale_fill_manual(values = c("neutral" = "#FFFFBF",
                               "positive" = "#91CF60",
                               "negative" = "#FC8D59")) +
  facet_wrap(~phenology, nrow=2, scale="free", strip.position = "right")+
  theme(legend.position="none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        title = element_blank(),
        #axis.text.y = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5, 0, 5, 0))

bar

combined_plot <- p + bar  + plot_layout(ncol = 2, widths = c(6, 1))
combined_plot
ggsave(combined_plot, filename="../Figures/Fig.1.BC.Map.with.Bar/Fig.1.BC.Map.with.Bar.png",
       width=8, height=6)
