library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
library(lme4)
library(performance)

setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
dist=10
paired.records.filtered<-readRDS(sprintf("../Data/paired.records.clean.%dkm.rda", dist))

DT<-data.table(e.id=paired.records.filtered$e.id,
               n.id=paired.records.filtered$n.id,
               type=paired.records.filtered$type, 
               e.days=paired.records.filtered$e.days, 
               bot_country=paired.records.filtered$e.bot_country, 
               e.phenology=paired.records.filtered$e.phenology, 
               n.phenology=paired.records.filtered$n.phenology, 
               decimalLatitude=abs(paired.records.filtered$e.decimalLatitude), 
               decimalLongitude=paired.records.filtered$e.decimalLongitude, 
               mean=paired.records.filtered$e.mean - paired.records.filtered$n.mean, 
               sd=paired.records.filtered$e.sd - paired.records.filtered$n.sd, 
               min=paired.records.filtered$e.min - paired.records.filtered$n.min,
               max=paired.records.filtered$e.max - paired.records.filtered$n.max, 
               range=paired.records.filtered$e.range - paired.records.filtered$n.range, 
               accu_temp=paired.records.filtered$e.accu_temp - paired.records.filtered$n.accu_temp, 
               pheno_arc=paired.records.filtered$e.pheno_arc_fixed-paired.records.filtered$n.pheno_arc_fixed,
               growthform=paired.records.filtered$growthform)

conditions<-data.table(expand.grid(type=unique(DT$type),
                                   e.days=unique(DT$e.days)))
coms<-data.table(expand.grid(phenology=c("Flowering", "Fruiting"),
                             growthform=c("All", "Herbaceous", "Woody")))
j=25
cor.all<-list()
for (j in c(1:nrow(conditions))){
  cond<-conditions[j]
  if (cond$type!="Both"){
    dt.item<-DT[type==cond$type & e.days==cond$e.days]
  }else{
    dt.item<-DT[e.days==cond$e.days]
  }
  i=1
  for (i in c(1:nrow(coms))){
    print(paste(i, nrow(coms), j, nrow(conditions)))
    com<-coms[i]
    item.paired.differ<-dt.item
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
    
    item.paired.differ$sd<-scale(item.paired.differ$sd)
    item.paired.differ$accu_temp<-scale(item.paired.differ$accu_temp)
    item.paired.differ$min<-scale(item.paired.differ$min)
    item.paired.differ$max<-scale(item.paired.differ$max)
    item.paired.differ$range<-scale(item.paired.differ$range)
    item.paired.differ$mean<-scale(item.paired.differ$mean)
    cor(item.paired.differ[, c("sd", "accu_temp", "min", "max", "range", "mean")])
    
    item.paired.differ$lat<-scale(abs(item.paired.differ$decimalLatitude))
   
    model <- lmer(pheno_arc ~ range+accu_temp + sd + min + max + mean + (1 | bot_country), 
                  data = item.paired.differ)
    check_collinearity(model)
    
    r2<-r2(model)
    vif(model)
    result<-summary(model)
    a<-data.table(result$coefficients)
    a$var<-c("intercept", "accu_temp", "sd")
    a$R2_conditional<-r2$R2_conditional
    a$R2_marginal<-r2$R2_marginal
    a$type<-cond$type
    a$e.days<-cond$e.days
    a$phenology<-com$phenology
    a$growthform<-com$growthform
    
    cor.all[[length(cor.all)+1]]<-a
    
  }
}
cor.df<-rbindlist(cor.all)

saveRDS(cor.df, "../Figures/Fig.3.t.test.BC/lmer.rda")
