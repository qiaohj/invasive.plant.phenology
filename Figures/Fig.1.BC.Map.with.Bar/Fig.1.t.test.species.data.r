library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
paired.differ<-readRDS("../Data/paired.differ.rda")

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
saveRDS(all.species.df, "../Figures/Fig.1.BC.Map.with.Bar/Fig.1.t.test.species.rda")
