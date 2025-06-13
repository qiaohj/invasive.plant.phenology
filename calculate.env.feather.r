library(data.table)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
iNadata_ll<-readRDS("../Data/iNadata_ll.rda")

all_v<-list()
for (i in c(1:nrow(iNadata_ll))){
  print(paste(i, nrow(iNadata_ll)))
  item<-iNadata_ll[i]
  target<-sprintf("../Data/Env.Items/%d.rda", item$id)
  if (!file.exists(target)){
    next()
  }
  v<-readRDS(target)
  v[type=="2m_temperature", value:=value-273.16]
  v[type=="skin_reservoir_content" & value<0, value:=0]
  v$group<-floor((v$layer.index-min(v$layer.index))/7)
  for (g in c(unique(v$group))){
    item.v<-v[group<=g]
    
    v_cal<-item.v[, .(
      id=item$id,
      group=g,
      days=.N,
      invasive=item$invasive,
      status=item$status,
      observed_on=item$observed_on,
      scientific_name=item$scientific_name,
      phenology=item$phenology,
      decimalLatitude=item$decimalLatitude,
      decimalLongitude=item$decimalLongitude,
      mean=mean(value),
      sd=sd(value),
      min=min(value),
      max=max(value),
      sum=sum(value),
      range=max(value)-min(value),
      accu_temp=sum(value * (value>0))),
      by=list(type)]
    v_cal$effe_accu<-0
    for (t in unique(item.v$type)){
      index<-min(item.v[type==t & value<0]$layer.index)
      item.v2<-item.v[type==t & layer.index<index]
      v_cal[type==t, effe_accu:=sum(item.v2$value)]
    }
    all_v[[length(all_v)+1]]<-v_cal
  }
}
all_v<-rbindlist(all_v)

saveRDS(all_v, "../Data/iNadata_ll_with_env.rda")


iNadata <- readRDS("../Data/iNadata_ll_with_env.rda")
iNadata <- iNadata[!is.na(mean) & phenology!=0 & observed_on >= as.Date("2008-01-01")]
iNadata[status=="exotic"]
iNadata[scientific_name=="Reynoutria japonica" & status=="native"]
colnames(iNadata)
table(iNadata$status)
records<-iNadata[,.(N=.N), by=.(id, invasive, status, scientific_name, 
                                decimalLatitude, decimalLongitude)]
points<-st_as_sf(
  records,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)
world<-read_sf("../Data/Similar_phenology_invasive_code&data_20250218/geographic boundary data/bot_country_robin.shp")
world_buffered <- st_buffer(world, dist = 20000)

iNadata.p<-st_as_sf(
  iNadata,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)
iNadata.p<-st_transform(iNadata.p, crs=st_crs(world))
iNadata$bot_country<-""
index<-st_intersects(world_buffered, iNadata.p)
for (i in c(1:nrow(index))){
  item<-index[[i]]
  if (length(item)>0){
    iNadata[item, bot_country:=world[i,]$LEVEL3_COD]
  }
}
table(iNadata$bot_country)

iNadata.p<-st_as_sf(
  unique(iNadata[, c("decimalLongitude", "decimalLatitude", "bot_country")]),
  coords = c("decimalLongitude", "decimalLatitude"),
  crs = 4326
)
iNadata.p<-st_transform(iNadata.p, crs=st_crs(world))

ggplot(world)+geom_sf()+
  geom_sf(data=iNadata.p[sample(nrow(iNadata.p), 1e4),], aes(color=bot_country))+
  theme(legend.position="none")
saveRDS(iNadata, "../Data/iNadata_ll_with_env_bot_country.rda")
