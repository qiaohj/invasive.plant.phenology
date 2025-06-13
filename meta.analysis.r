library(tidyverse)
library(sf)
library(castor)
library(lubridate)
library(U.PhyloMaker)
library(do)
library(metafor)
library(ggbeeswarm)
library(stringr)
library(knitr)
library(broom)
library(kableExtra)
library(data.table)
library(car)
library(ggplot2)
library(ggcorrplot)
library(lme4)

setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
dist<-10
bot_country <- st_read("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/bot_country_robin.shp") ## Transformed boundary of botanical country (Robinson projection)

#Defining the function to calculate the sd and mean of circular vector

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
paired.records.filtered$growthform <- if_else(paired.records.filtered$e.invasive %in% herbaceous, "Herbaceous", "Woody")
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
paired.records.filtered[pheno_arc==-360]

paired.differ<-unique(paired.records.filtered[, 
                                              c("e.invasive", "e.status", "e.scientific_name",
                                                "e.phenology", "e.bot_country", 
                                                "n.invasive", "n.status", "n.scientific_name",
                                                "n.phenology", "n.bot_country", 
                                                "growthform", "e.pheno_arc_fixed",
                                                "n.pheno_arc_fixed")])

if (F){
  dfff<-paired.records.filtered[, c("e.days_in_year", "e.decimalLatitude",
                                    "e.phenology")]
  dfff<-unique(dfff)
  dfff$hemi<-ifelse(dfff$e.decimalLatitude>0, "N", "S")
  ggplot(dfff)+geom_histogram(aes(x=e.days_in_year, fill=hemi))+
    geom_vline(xintercept = 183, linetype=2)+
    facet_wrap(~e.phenology)
}
coms<-data.table(expand.grid(phenology=c("Phenology", "Flowering", "Fruiting"),
                             growthform=c("All", "Herbaceous", "Woody")))
i=1
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
  table(item.paired.differ$growthform)
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
  
  bot_country_v<-merge(bot_country, all, by.x="LEVEL3_COD", by.y="bot_country")
  p<-ggplot()+
    geom_sf(data=bot_country, fill=NA, color="lightgrey", alpha=0.3)+
    geom_sf(data=bot_country_v, aes(fill=direction.factor))+
    labs(title=paste(com$phenology, com$growthform, sep="+"))+
    theme_bw()
  
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


DT<-data.table(e.id=paired.records.filtered$e.id,
               n.id=paired.records.filtered$n.id,
               type=paired.records.filtered$type, 
               e.days=paired.records.filtered$e.days, 
               bot_country=paired.records.filtered$e.bot_country, 
               e.phenology=paired.records.filtered$e.phenology, 
               n.phenology=paired.records.filtered$n.phenology, 
               decimalLatitude=paired.records.filtered$e.decimalLatitude, 
               decimalLongitude=paired.records.filtered$e.decimalLongitude, 
               mean=paired.records.filtered$e.mean - paired.records.filtered$n.mean, 
               sd=paired.records.filtered$e.sd - paired.records.filtered$n.sd, 
               min=paired.records.filtered$e.min - paired.records.filtered$n.min,
               max=paired.records.filtered$e.max - paired.records.filtered$n.max, 
               range=paired.records.filtered$e.range - paired.records.filtered$n.range, 
               accu_temp=paired.records.filtered$e.accu_temp - paired.records.filtered$n.accu_temp, 
               pheno_arc=paired.records.filtered$e.pheno_arc_fixed-paired.records.filtered$n.pheno_arc_fixed,
               growthform=paired.records.filtered$growthform)

conditions<-data.table(expand.grid(type=c("Both", unique(DT$type)),
                                   lat=c("with", "without"),
                                   bc=c("with", "without"),
                                   e.days=unique(DT$e.days)))
j=41
cor.all<-list()
for (j in c(1:nrow(conditions))){
  cond<-conditions[j]
  if (cond$type!="Both"){
    dt.item<-DT[type==cond$type & e.days==cond$e.days]
  }else{
    dt.item<-DT[e.days==cond$e.days]
  }
  i=2
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
    if (F){
      model <- lm(accu_temp~sd+range , data = item.paired.differ)
      vif_values <- vif(model)
      vif_dt <- data.table(var = names(vif_values), VIF = as.numeric(vif_values))
      
      p1 <- ggplot(vif_dt, aes(x = reorder(var, VIF), y = VIF)) +
        geom_col(fill = "steelblue") +
        geom_hline(yintercept = c(5, 10), linetype = "dashed", color = "red") +
        theme_minimal() +
        labs(title = "VIF barplot", x = "Variable", y = "VIF")
      p1
      
      
      mat_cor <- cor(item.paired.differ[, c("sd", "accu_temp")])
      if (F){
        xxx<-item.paired.differ[sample(nrow(item.paired.differ), 1000)]
        plot(xxx$range, xxx$sd)
      }
      p.mat <- cor_pmat(item.paired.differ[, c("mean", "sd", "range", "accu_temp")])
      p2 <- ggcorrplot(mat_cor, hc.order = TRUE, type = "lower",
                       lab = TRUE, p.mat = p.mat, insig = "blank") +
        ggtitle("Correlation matrix heatmap")
      
      p2
    }
    if (cond$type=="Both"){
      model.item1<-item.paired.differ[type=="2m_temperature"]
      model.item2<-item.paired.differ[type=="skin_reservoir_content"]
      model.item1<-model.item1[, c("decimalLatitude", "decimalLongitude",
                                   "mean", "sd", "range", "accu_temp", "pheno_arc",
                                   "e.id", "n.id", "bot_country")]
      model.item2<-model.item2[, c("decimalLatitude", "decimalLongitude",
                                   "mean", "sd", "range", "accu_temp", "pheno_arc",
                                   "e.id", "n.id", "bot_country")]
      colnames(model.item1)<-c("decimalLatitude", "decimalLongitude",
                               "t.mean", "t.sd", "t.range", "t.accu_temp", "pheno_arc",
                               "e.id", "n.id", "bot_country")
      colnames(model.item2)<-c("decimalLatitude", "decimalLongitude",
                               "p.mean", "p.sd", "p.range", "p.accu_temp", "pheno_arc",
                               "e.id", "n.id", "bot_country")
      model.dt<-merge(model.item1, model.item2, 
                      by=c("e.id", "n.id", "decimalLatitude", "decimalLongitude",
                           "pheno_arc", "bot_country"))
      model.dt$t.sd<-scale(model.dt$t.sd)
      model.dt$t.accu_temp<-scale(model.dt$t.accu_temp)
      model.dt$p.sd<-scale(model.dt$p.sd)
      model.dt$p.accu_temp<-scale(model.dt$p.accu_temp)
      model.dt$decimalLatitude<-scale(abs(model.dt$decimalLatitude))
      model.dt$pheno_arc<-scale(model.dt$pheno_arc)
      if (cond$lat=="with"){
        if (cond$bc=="with"){
          model <- lmer(pheno_arc ~ t.accu_temp + t.sd + p.accu_temp + p.sd + decimalLatitude +(1 | bot_country), 
                        data = model.dt)
        }else{
          model <- lm(pheno_arc ~ t.accu_temp + t.sd + p.accu_temp + p.sd + decimalLatitude, 
                        data = model.dt)
        }
        
      }else{
        if (cond$bc=="with"){
          model <- lmer(pheno_arc ~ t.accu_temp + t.sd + p.accu_temp + p.sd + (1 | bot_country), 
                        data = model.dt)
        }else{
          model <- lm(pheno_arc ~ t.accu_temp + t.sd + p.accu_temp + p.sd, 
                        data = model.dt)
        }
        
      }
      model.dt$pred <- predict(model)
      cor<-cor(model.dt$pheno_arc, model.dt$pred)

    }else{
      item.paired.differ$sd<-scale(item.paired.differ$sd)
      item.paired.differ$accu_temp<-scale(item.paired.differ$accu_temp)
      item.paired.differ$min<-scale(item.paired.differ$min)
      item.paired.differ$max<-scale(item.paired.differ$max)
      
      item.paired.differ$decimalLatitude<-scale(abs(item.paired.differ$decimalLatitude))
      #item.paired.differ$pheno_arc<-scale(item.paired.differ$pheno_arc)
      if (cond$lat=="with"){
        if (cond$bc=="with"){
          model <- lmer(pheno_arc ~ accu_temp + sd + decimalLatitude + (1 | bot_country), 
                    data = item.paired.differ)
        }else{
          model <- lm(pheno_arc ~ accu_temp + sd + decimalLatitude + 1, 
                        data = item.paired.differ)
        }
      }else{
        if (cond$bc=="with"){
          model <- lmer(pheno_arc ~ accu_temp + sd + min + (1 | bot_country), 
                        data = item.paired.differ)
        }else{
          model <- lm(pheno_arc ~ accu_temp + sd, 
                        data = item.paired.differ)
        }
      }
      item.paired.differ$pred <- predict(model)
      #plot(item.paired.differ$pheno_arc, item.paired.differ$pred)
      cor<-cor(item.paired.differ$pheno_arc, item.paired.differ$pred)
      
      result<-summary(model)
      a<-data.table(result$coefficients)
      a<-a[2]$Estimate
      line<-data.table(x=seq(min(item.paired.differ$accu_temp), 
                             max(item.paired.differ$accu_temp), by=1))
      line$y<-line$x * a
      ggplot(item.paired.differ[sample(nrow(item.paired.differ), 1e3)])+
        geom_point(aes(y=pheno_arc, x=accu_temp), color="lightgrey", alpha=0.5)+
        geom_line(data=line, aes(x=x, y=y))
    }
    cor.item<-data.table(type=cond$type, e.days=cond$e.days, lat=cond$lat, bc=cond$bc,
                        phenology=com$phenology, growthform=com$growthform, 
                        cor=cor)
    #summary(model)
    
    
    cor.all[[length(cor.all)+1]]<-cor.item
    if (F){
      ggplot(item.paired.differ[sample(nrow(item.paired.differ), 1000)], 
             aes(y = pred, x = pheno_arc, color = bot_country)) +
        geom_point() +
        #geom_line(aes(y = pred)) +
        theme_minimal() +
        labs(title = "Mixed Effects Model: Random Intercepts by Group")
      
      ggplot(item.paired.differ[sample(nrow(item.paired.differ), 1000)], aes(x = pred, y = pheno_arc)) +
        geom_point(alpha = 0.5) +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        theme_minimal() +
        labs(title = "Observed vs Predicted")
      
      summary(model)
      
      fixef(model)
      summarise(model)
      
    }
  }
}
cor.df<-rbindlist(cor.all)
saveRDS(cor.df, "../Data/lmer.rda")
ggplot(cor.df)+geom_point(aes(x=e.days, y=cor, color=type, shape=phenology))+
  facet_grid(lat~growthform+bc)

ggplot(cor.df[bc=="with" & phenology!="Phenology" &
                type!="Both" & lat=="without"])+
  geom_point(aes(x=e.days, y=cor, color=type, shape=phenology))+
  facet_wrap(~growthform, nrow=3)


df_sampled<-paired.records.filtered[sample(nrow(paired.records.filtered), 1e4)]
df_sampled$phenology<-ifelse(df_sampled$e.phenology==1, "Flowering", "Fruiting")
ggplot(df_sampled)+geom_quasirandom(aes(x=pheno_are, y=phenology, fill = phenology))
