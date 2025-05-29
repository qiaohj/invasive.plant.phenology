library(terra)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")

netcdf<-list.files("/media/huijieqiao/WD22T_11/CDS/netcdf/", pattern="\\.nc")
i=1
netcdf<-netcdf[sample(length(netcdf), length(netcdf))]
for (i in c(1:length(netcdf))){
  print(paste(i, length(netcdf)))
  item<-netcdf[i]
  target<-sprintf("/media/huijieqiao/WD22T_11/CDS/GeoTIFF/%s", gsub("\\.nc", "\\.tif", item))
  if (file.exists(target)){
    next()
  }
  saveRDS(NULL, target)
  
  r<-rast(sprintf("/media/huijieqiao/WD22T_11/CDS/netcdf/%s", item))
  r<-rotate(r)
  writeRaster(r, target, overwrite=T)
}