library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
dist<-10

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



paired.records.filtered<-readRDS(sprintf("../Data/paired.records.clean.%dkm.rda", dist))

paired.differ<-unique(paired.records.filtered[, 
                                              c("e.invasive", "e.status", "e.scientific_name",
                                                "e.phenology", "e.bot_country", 
                                                "n.invasive", "n.status", "n.scientific_name",
                                                "n.phenology", "n.bot_country", 
                                                "growthform", "e.pheno_arc_fixed",
                                                "n.pheno_arc_fixed")])
saveRDS(paired.differ, "../Data/paired.differ.rda")
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

all.df$direction.factor.sc<-factor(all.df$direction.factor,
                                   levels=c("Less", "No difference", "Greater"),
                                   labels=c("negative", "neutral", "positive"))
saveRDS(all.df, "../Figures/Fig.1.BC.Map.with.Bar/Fig.1.BC.Map.with.Bar.A.rda")

