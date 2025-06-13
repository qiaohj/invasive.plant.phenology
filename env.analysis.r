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


setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
if (F){
  iNadata <- readRDS("../Data/iNadata_ll_with_env_bot_country.rda")
  
  
  iNadata$days_in_year<-as.numeric(format(iNadata$observed_on, "%j"))
  records<-iNadata[,.(N=.N), by=.(id, invasive, status, scientific_name, 
                                  decimalLatitude, decimalLongitude, bot_country)]
  
  native<-records[status=="native"]
  native_p<-st_as_sf(
    native,
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  )
  native_p <- st_transform(native_p, 3857)
  
  
  exotic<-records[status=="exotic"]
  exotic_p<-st_as_sf(
    exotic,
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326
  )
  exotic_p <- st_transform(exotic_p, 3857)
  
  dists<-c(40e3)
  i=1
  all<-list()
  for (dist in dists){
    for (i in c(1:nrow(exotic))){
      print(paste(i, nrow(exotic), dist))
      item<-exotic[i]
      buffer<-st_buffer(exotic_p[i,], dist=dist)
      relevant<-unique(native[invasive==item$scientific_name]$scientific_name)
      native_in_buffer <- native_p[st_contains(buffer, native_p)[[1]],]
      native_in_buffer <- native_in_buffer[which(native_in_buffer$scientific_name %in% relevant),]
      if (nrow(native_in_buffer)==0){
        next()
      }
      exotic_record<-iNadata[id==item$id]
      native_record<-iNadata[id %in% native_in_buffer$id]
      colnames(exotic_record)[c(2, 4:21)]<-
        paste("e", colnames(exotic_record)[c(2, 4:21)], sep=".")
      colnames(native_record)[c(2, 4:21)]<-
        paste("n", colnames(native_record)[c(2, 4:21)], sep=".")
      record_item<-merge(exotic_record, native_record, by=c("type", "group"))
      record_item$dist<-dist
      all[[length(all)+1]]<-record_item
      if (F){
        ggplot(buffer)+geom_sf()+geom_sf(data=native_in_buffer)
      }
    }
    all_df<-rbindlist(all)
    
    saveRDS(all_df, sprintf("../Data/paired.records.%dkm.rda", dist/1000))
  }
  
}

paired.records<-readRDS("../Data/paired.records.rda")

herbaceous <- c("Reynoutria japonica", 
                "Hedychium gardnerianum",
                "Pueraria montana",
                "Mikania micrantha", 
                "Lythrum salicaria",
                "Chromolaena odorata",
                "Sphagneticola trilobata")
paired.records$growthform <- if_else(paired.records$e.invasive %in% herbaceous, "Herbaceous", "Woody")

result.flower <- paired.records[
  e.phenology %in% c(1, 12) & n.phenology %in% c(1,12), 
  {
    valid_idx <- complete.cases(e.mean, n.mean)
    n_valid <- sum(valid_idx)
    
    if (n_valid >= 2) {
    test <- wilcox.test(e.mean, n.mean, paired = F)
    
    list(
      p_value = test$p.value,
      e.mean = mean(e.mean, na.rm = TRUE),
      n.mean = mean(n.mean, na.rm = TRUE),
      direction = if (test$p.value < 0.05) {
        if (mean(e.mean, na.rm = TRUE) > mean(n.mean, na.rm = TRUE)) "e > n"
        else if (mean(e.mean, na.rm = TRUE) < mean(n.mean, na.rm = TRUE)) "e < n"
        else "No Difference"
      } else {
        "Not Significant"
      }
    )
    }else{
      list(
        p_value = as.numeric(NA),
        e.mean = as.numeric(NA),
        n.mean = as.numeric(NA),
        direction = "Too Few Observations"
      )
    }
  }, by = .(type, group, growthform)]


xxx<-result.flower[, .(N=.N), by=c("type", "group", "direction", "growthform")]
xxx[group==1]
g<-0
var<-"2m_temperature"

bot_country <- st_read("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/bot_country_robin.shp") ## Transformed boundary of botanical country (Robinson projection)
bot_country_v<-merge(bot_country, result.flower, by.x="LEVEL3_COD", by.y="bot_country")
ggplot(bot_country_v[which(bot_country_v$group==0),])+geom_sf(aes(fill=direction))+
  facet_grid(growthform~type)+
  theme_bw()

