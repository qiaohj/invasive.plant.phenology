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
library(rinat)
library(terra)
library(httr)
library(jsonlite)



setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")

if (F){
  iNadata <- fread("../Data/Similar_phenology_invasive_code&data_20250218/iNadata.csv")
  iNatRaw<-readRDS("../Data/iNaturalist.rda")
  
  iNadata$occurrenceID<-sprintf("https://www.inaturalist.org/observations/%d", iNadata$id)
  
  nrow(iNadata[occurrenceID %in% iNatRaw$occurrenceID])
  iNadata[!occurrenceID %in% iNatRaw$occurrenceID]
  
  iNatRaw_LL<-iNatRaw[, c("occurrenceID", "decimalLongitude", "decimalLatitude")]
  iNatRaw_LL<-iNatRaw_LL[!is.na(decimalLongitude) & !is.na(decimalLatitude)]
  
  iNadata_ll<-merge(iNadata, iNatRaw_LL, by="occurrenceID", all.x=T)
  get_inat_obs_id(33170675)
  
  iNatRaw
  i=1
  
  iNadata_ll<-readRDS("../Data/iNadata_ll.rda")
  iNadata_ll$is.na<-is.na(iNadata_ll$decimalLongitude)
  ids<-iNadata_ll[is.na(decimalLongitude)]$id
  results <- list()
  
  chunks <- split(ids, ceiling(seq_along(ids) / 200))
  
  for (ii in c(1:length(chunks))) {
    chunk<-chunks[[ii]]
    print(paste(ii, length(chunks)))
    id_str <- paste(chunk, collapse = ",")
    url <- paste0("https://api.inaturalist.org/v1/observations?id=", id_str)
    #print(url)
    resp <- GET(url)
    if (status_code(resp) != 200) {
      warning("Error request ", url)
      next
    }
    
    json <- fromJSON(content(resp, "text", encoding = "UTF-8"))
    
    if (length(json$results) == 0) next
    
    ll <- rbindlist(lapply(json$results$geojson$coordinates, function(x) data.table(lon = x[1], lat = x[2])))
    ll$id<-json$results$id
    
    
    results[[length(results) + 1]] <- ll
  }
  
  lll<-rbindlist(results, fill = TRUE)
  
  saveRDS(iNadata_ll, "../Data/iNadata_ll.rda")
}
iNadata_ll<-readRDS("../Data/iNadata_ll.rda")
i=1
r_cache <- new.env(parent = emptyenv())
get_monthly_rast <- function(tif_path) {
  if (! exists(tif_path, envir = r_cache) ) {
    r_cache[[tif_path]] <- rast(tif_path)
  }
  r_cache[[tif_path]]
}
iNadata_ll<-iNadata_ll[!is.na(decimalLongitude)]
iNadata_ll<-iNadata_ll[sample(nrow(iNadata_ll), nrow(iNadata_ll))]
for (i in c(1:nrow(iNadata_ll))){
  print(paste(i, nrow(iNadata_ll)))
  item<-iNadata_ll[i]
  target<-sprintf("../Data/Env.Items/%d.rda", item$id)
  if (file.exists(target)){
    next()
  }
  saveRDS(NULL, target)
  end_date <- item$observed_on
  date_seq <- seq(end_date - 89, end_date, by = "day")
  
  dt <- CJ(date = date_seq, 
           file = c("/media/huijieqiao/WD22T_11/CDS/GeoTIFF/CDS.derived-era5-land-daily-statistics.y%d.m%d.daily_mean.2m_temperature.tif", 
                    "/media/huijieqiao/WD22T_11/CDS/GeoTIFF/CDS.derived-era5-land-daily-statistics.y%d.m%d.daily_mean.skin_reservoir_content.tif"),
           lon=item$decimalLongitude, lat=item$decimalLatitude)
  
  
  min.date<-min(dt$date)
  dt$layer.index<-as.numeric(as.Date(dt$date)- as.Date(sprintf("%d-%d-1", 
                                                    as.numeric(format(min.date, "%Y")),
                                                    as.numeric(format(min.date, "%m"))))+1)
  
  dt[, `:=`(
    year = as.numeric(format(date, "%Y")),
    month = as.numeric(format(date, "%m")),
    day = as.numeric(format(date, "%d"))
  )]
  
  dt$file<-sprintf(dt$file, dt$year, dt$month)
  files<-unique(dt$file)
  dt$value<--9999
  for (f in files){
    r<-get_monthly_rast(f)
    values<-extract(r, item[, c("decimalLongitude", "decimalLatitude")], raw=T, ID=F)
    dt[file==f, value:=values[dt[file==f]$day]]
  }
  
  dt$type<-ifelse(grepl("temperature", dt$file), "2m_temperature", "skin_reservoir_content")
  dt$file<-NULL
  
  saveRDS(dt, target)
}
