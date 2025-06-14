library(data.table)
library(ggplot2)
library(patchwork)

library(sf)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
dist<-10
bot_country <- st_read("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/bot_country_robin.shp") ## Transformed boundary of botanical country (Robinson projection)
bbox <- st_read("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/worldmap_boundary_box.shp") ## World map boundary
bbox_robin <- st_transform(bbox, crs="+proj=robin")

if (F){
  
  paired.records<-readRDS(sprintf("../Data/paired.records.%dkm.rda", dist))
  paired.records.filtered<-paired.records[e.phenology==n.phenology |
                                            e.phenology==12 |
                                            n.phenology==12]
  paired.records.filtered<-paired.records.filtered[e.phenology!=0 & n.phenology!=0]
  paired.records.filtered<-paired.records.filtered[!((e.phenology==1 & n.phenology==2) | 
                                                       (e.phenology==2 & n.phenology==1))]
  herbaceous <- c("Reynoutria japonica", 
                  "Hedychium gardnerianum",
                  "Pueraria montana",
                  "Mikania micrantha", 
                  "Lythrum salicaria",
                  "Chromolaena odorata",
                  "Sphagneticola trilobata")
  paired.records.filtered$growthform <- ifelse(paired.records.filtered$e.invasive %in% herbaceous, "Herbaceous", "Woody")
  paired.records.filtered$e.pheno_arc<-as.numeric(format(paired.records.filtered$e.observed_on, format = "%j")) * 360 / 365
  paired.records.filtered$n.pheno_arc<-as.numeric(format(paired.records.filtered$n.observed_on, format = "%j")) * 360 / 365
  paired.records.filtered$pheno_diff_abs<-abs(paired.records.filtered$e.pheno_arc - paired.records.filtered$n.pheno_arc)
  paired.records.filtered$e.pheno_arc_fixed<-paired.records.filtered$e.pheno_arc
  paired.records.filtered$n.pheno_arc_fixed<-paired.records.filtered$n.pheno_arc
  paired.records.filtered[e.pheno_arc>n.pheno_arc & pheno_diff_abs>180, e.pheno_arc_fixed:=
                            paired.records.filtered[e.pheno_arc>n.pheno_arc & pheno_diff_abs>180, "e.pheno_arc"]-360]
  paired.records.filtered[e.pheno_arc<n.pheno_arc & pheno_diff_abs>180, n.pheno_arc_fixed:=
                            paired.records.filtered[e.pheno_arc<n.pheno_arc & pheno_diff_abs>180, "n.pheno_arc"]-360]
  paired.records.filtered$pheno_arc<-paired.records.filtered$e.pheno_arc_fixed - 
    paired.records.filtered$n.pheno_arc_fixed
  saveRDS(paired.records.filtered, sprintf("../Data/paired.records.clean.%dkm.rda", dist))
}

paired.records.filtered<-readRDS(sprintf("../Data/paired.records.clean.%dkm.rda", dist))

#The figure for N of records per hemisphere by days of year
if (F){
  dfff.e<-paired.records.filtered[, c("e.days_in_year", "e.decimalLatitude",
                                    "e.phenology", "growthform")]
  colnames(dfff.e)<-c("days_in_year", "decimalLatitude",
                      "phenology", "growthform")
  dfff.e$type<-"exotic"
  dfff.n<-paired.records.filtered[, c("n.days_in_year", "n.decimalLatitude",
                                      "n.phenology", "growthform")]
  colnames(dfff.n)<-c("days_in_year", "decimalLatitude",
                      "phenology", "growthform")
  dfff.n$type<-"native"
  
  dfff<-rbindlist(list(unique(dfff.n), unique(dfff.e)))
  
  dfff1<-dfff[phenology==12]
  dfff1.1<-dfff1
  dfff1.1$phenology<-1
  dfff1.2<-dfff1
  dfff1.2$phenology<-2
  
  dfff2<-dfff[phenology!=12]
  dfff<-rbindlist(list(dfff1.1, dfff1.2, dfff2))
  dfff$hemi<-ifelse(dfff$decimalLatitude>0, "N", "S")
  dfff$month<-dfff$days_in_year/30+1
  dfff$phenology.l<-ifelse(dfff$phenology==1, "Flowering", "Fruiting")
  p<-ggplot(dfff)+geom_histogram(aes(x=month, fill=hemi), bins=50)+
    geom_vline(xintercept = 7, linetype=2, color="grey")+
    facet_grid(growthform~phenology.l)+
    labs(fill="Hemisphere", x="Month", y="")+
    theme_bw()
  p
  ggsave(p, filename="../Figures/records.hist.png", width=7, height=6)
}

paired.differ<-unique(paired.records.filtered[, 
                                              c("e.invasive", "e.status", "e.scientific_name",
                                                "e.phenology", "e.bot_country", 
                                                "n.invasive", "n.status", "n.scientific_name",
                                                "n.phenology", "n.bot_country", 
                                                "growthform", "e.pheno_arc_fixed",
                                                "n.pheno_arc_fixed")])

coms<-data.table(expand.grid(phenology=c("Flowering", "Fruiting"),
                             growthform=c("All")))
i=1
all.df<-list()
for (i in c(1:nrow(coms))){
  print(paste(i, nrow(coms)))
  com<-coms[i]
  item.paired.differ<-paired.differ
  if (com$phenology=="Flowering"){
    item.paired.differ<-item.paired.differ[n.phenology %in% c(1, 12) &
                                             e.phenology %in% c(1, 12)]
  }
  if (com$phenology=="Fruiting"){
    item.paired.differ<-item.paired.differ[n.phenology %in% c(2, 12) &
                                             e.phenology %in% c(2, 12)]
  }
  if (com$growthform=="Herbaceous"){
    item.paired.differ<-item.paired.differ[growthform=="Herbaceous"]
  }
  if (com$growthform=="Woody"){
    item.paired.differ<-item.paired.differ[growthform=="Woody"]
  }
  
  bcs<-item.paired.differ[, .(N=.N), by=c("e.bot_country")]
  bc<-bcs[1]
  j=1
  #table(item.paired.differ$growthform)
  all<-list()
  for (j in c(1:nrow(bcs))){
    bc<-bcs[j]
    item<-item.paired.differ[e.bot_country==bc$e.bot_country]
    if (nrow(item)<=1){
      next()
    }
    for (direction in c("two.sided")){
      model<-t.test(item$e.pheno_arc_fixed, item$n.pheno_arc_fixed, alternative=direction, paired = T)
      result<-data.table(p.value=model$p.value,
                         ci.low=model$conf.int[1],
                         ci.high=model$conf.int[2],
                         t=model$statistic,
                         df=model$parameter,
                         mean.difference=model$estimate,
                         stderr=model$stderr,
                         direction=direction,
                         bot_country=bc$e.bot_country,
                         N=bc$N)
      all[[length(all)+1]]<-result
    }
    
  }
  
  all<-rbindlist(all)
  all$direction.factor<-"No difference"
  all[p.value<=0.05 & mean.difference>0, direction.factor:="Greater"]
  all[p.value<=0.05 & mean.difference<0, direction.factor:="Less"]
  
  all$direction.factor<-factor(all$direction.factor, 
                               levels=c("Greater", "No difference", "Less"), 
                               labels=c("Greater", "No difference", "Less"))
  all$phenology<-com$phenology
  all.df[[length(all.df)+1]]<-all
  
  if (F){
    bot_country_v<-merge(bot_country, all, by.x="LEVEL3_COD", by.y="bot_country")
    p<-ggplot()+
      geom_sf(data=bot_country, fill=NA, color="lightgrey", alpha=0.3)+
      geom_sf(data=bot_country_v, aes(fill=direction.factor))+
      labs(title=paste(com$phenology, com$growthform, sep="+"))+
      theme_bw()
    p
    ggsave(p, filename=sprintf("../Figures/bc.map/bc.%s.%s.png", com$phenology, com$growthform), width=9, height=7)
    setorderv(all, c("direction.factor", "mean.difference"), c(1,1))
    all[bot_country=="SRL"]
    all$bot_country <- factor(all$bot_country, levels = all$bot_country)
    p<-ggplot(all)+
      geom_vline(aes(xintercept=0), linetype=2, color="black", alpha=0.5)+
      geom_errorbarh(aes(xmin=ci.low, xmax=ci.high, y=bot_country, color=direction.factor))+
      geom_point(aes(x=mean.difference, y=bot_country))+
      labs(title=paste(com$phenology, com$growthform, sep="+"))+
      theme_bw()
    ggsave(p, filename=sprintf("../Figures/bc.map/t.test.%s.%s.png", com$phenology, com$growthform), width=6, height=10)
  }
}
all.df<-rbindlist(all.df)
bot_country_v<-merge(bot_country, all.df, by.x="LEVEL3_COD", by.y="bot_country")


coms.item<-paired.differ[, .(N=.N), by=c("e.invasive",
                        "growthform")]
coms<-rbind(
  copy(coms.item)[, phenology := "Flowering"],
  copy(coms.item)[, phenology := "Fruiting"]
)
i=1
all.species<-list()
for (i in c(1:nrow(coms))){
  print(paste(i, nrow(coms)))
  com<-coms[i]
  item.paired.differ<-paired.differ
  if (com$phenology=="Flowering"){
    item.paired.differ<-item.paired.differ[n.phenology %in% c(1, 12) &
                                             e.phenology %in% c(1, 12)]
  }
  if (com$phenology=="Fruiting"){
    item.paired.differ<-item.paired.differ[n.phenology %in% c(2, 12) &
                                             e.phenology %in% c(2, 12)]
  }
  item.paired.differ<-item.paired.differ[e.invasive==com$e.invasive]
  if (nrow(item.paired.differ)<=0){
    next()
  }
  model<-t.test(item.paired.differ$e.pheno_arc_fixed, 
                item.paired.differ$n.pheno_arc_fixed, alternative="two.sided", paired = T)
  result<-data.table(p.value=model$p.value,
                     ci.low=model$conf.int[1],
                     ci.high=model$conf.int[2],
                     t=model$statistic,
                     df=model$parameter,
                     mean.difference=model$estimate,
                     stderr=model$stderr)
  
  result$direction.factor<-"No difference"
  result[p.value<=0.05 & mean.difference>0, direction.factor:="Greater"]
  result[p.value<=0.05 & mean.difference<0, direction.factor:="Less"]
  
  result$phenology<-com$phenology
  result$species<-com$e.invasive
  result$growthform<-com$growthform
  
  all.species[[length(all.species)+1]]<-result
 
}
all.species.df<-rbindlist(all.species)
all.species.df$direction.factor<-factor(all.species.df$direction.factor, 
                                levels=c("Less", "No difference", "Greater"), 
                                labels=c("Less", "No difference", "Greater"))


all.species.df[species=="Mimosa pigra" & phenology=="Flowering", 
               ci.low:=min(all.species.df[species!="Mimosa pigra" & phenology=="Flowering"]$ci.low)]
all.species.df[species=="Mimosa pigra" & phenology=="Flowering", 
               ci.high:=max(all.species.df[species!="Mimosa pigra" & phenology=="Flowering"]$ci.high)]
all.species.df$lt<-"1"
all.species.df[species=="Mimosa pigra" & phenology=="Flowering", lt:="2"]
setorderv(all.species.df, c("growthform", "mean.difference"), c(1,1))
all.species.df$species <- factor(all.species.df$species, levels = unique(all.species.df$species))

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

all.df$direction.factor.sc<-factor(all.df$direction.factor,
                                   levels=c("Less", "No difference", "Greater"),
                                   labels=c("negative", "neutral", "positive"))
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

combined_plot <- p + bar + t.test.species + plot_layout(ncol = 3, 
                                                        widths = c(6, 1, 8))
combined_plot

library(ggalluvial)

alluv.df<-paired.differ[,.(N=.N), by=.(e.invasive, e.bot_country, n.scientific_name)]

ggplot(alluv.df, aes(axis1 = e.bot_country,
               axis2 = e.invasive,
               axis3 = n.scientific_name,
               y = N)) +
  geom_alluvium(aes(fill = e.invasive), width = 1/12, alpha = 0.7) +
  geom_stratum(width = 1/12, fill = "grey90", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Invasive Species", "Country", "Native Species"), expand = c(.05, .05)) +
  theme(legend.position = "none")

