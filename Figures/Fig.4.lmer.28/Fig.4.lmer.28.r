library(lme4)
library(ggplot2)
library(data.table)
library(ggeffects)
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
conditions<-conditions[e.days==28 & type=="2m_temperature"]
coms<-data.table(expand.grid(phenology=c("Flowering", "Fruiting"),
                             growthform=c("Herbaceous", "Woody")))
j=1
point.data.all<-list()
line.data.all<-list()
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
    
    model <- lmer(pheno_arc ~ accu_temp + sd + (1 | bot_country), 
                  data = item.paired.differ)
    check_collinearity(model)
    
    r2<-r2(model)
    
    pred_effects.sd <- ggpredict(model, terms = c("sd"),
                              condition = c(accu_temp = mean(item.paired.differ$accu_temp)))
    line.data.sd<-data.table(var="sd", x=pred_effects.sd$x, y=pred_effects.sd$predicted,
                             ci.low=pred_effects.sd$conf.low, ci.high=pred_effects.sd$conf.high)
    
    pred_effects.accu_temp <- ggpredict(model, terms = c("accu_temp"),
                                 condition = c(sd = mean(item.paired.differ$sd)))
    line.data.accu_temp<-data.table(var="accu_temp", x=pred_effects.accu_temp$x, y=pred_effects.accu_temp$predicted,
                             ci.low=pred_effects.accu_temp$conf.low, ci.high=pred_effects.accu_temp$conf.high)
    
    
    line.data <- rbindlist(list(line.data.sd, line.data.accu_temp))
    sample.size<-ifelse(nrow(item.paired.differ)>1e3, 1e3, nrow(item.paired.differ))
    point.data<-item.paired.differ[sample(nrow(item.paired.differ), sample.size)]
    
    point.data.sd<-point.data[, c("sd", "pheno_arc")]
    colnames(point.data.sd)<-c("x", "y")
    point.data.sd$var<-"sd"
    
    point.data.accu_temp<-point.data[, c("accu_temp", "pheno_arc")]
    colnames(point.data.accu_temp)<-c("x", "y")
    point.data.accu_temp$var<-"accu_temp"
    point.data<-rbindlist(list(point.data.sd, point.data.accu_temp))
    
    point.data$type<-cond$type
    point.data$phenology<-com$phenology
    point.data$growthform<-com$growthform
    point.data.all[[length(point.data.all)+1]]<-point.data
    
    line.data$type<-cond$type
    line.data$phenology<-com$phenology
    line.data$growthform<-com$growthform
    line.data.all[[length(line.data.all)+1]]<-line.data
    
    if (F){
      ggplot(point.data, aes(x = x, y = y)) +
        geom_point(data=point.data, alpha = 0.2, color = "gray40") +
        geom_line(
          data = line.data,
          color = "deepskyblue3",
          linewidth = 1.2
        ) +
        facet_wrap(~ var, scales = "free_x") +
        theme(
          strip.background = element_rect(fill = "gray85"),
          strip.text = element_text(face = "bold")
        )
    }
  }
}
line.data.df[var=="sd" & phenology=="Flowering" & growthform=="Herbaceous"]

line.data.df<-rbindlist(line.data.all)
point.data.df<-rbindlist(point.data.all)
p<-ggplot(point.data.df, aes(x = x, y = y)) +
  geom_point(data=point.data.df, alpha = 0.1, color = "gray40") +
  geom_ribbon(data = line.data.df, aes(x=x, ymin=ci.low, ymax=ci.high), 
              fill="deepskyblue3", alpha=0.3)+
  geom_line(
    data = line.data.df,
    color = "deepskyblue3",
    linewidth = 1
  ) +
  facet_grid(phenology+growthform~ var+type, scales = "free") +
  theme(
    strip.background = element_rect(fill = "gray85"),
    strip.text = element_text(face = "bold")
  )
ggsave(p, filename="../Figures/Fig.4.lmer.28/Fig.4.lmer.28.png", width=6, height=8)
