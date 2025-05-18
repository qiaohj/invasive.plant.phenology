library(terra)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")

r<-rast("/media/huijieqiao/WD22T_11/CDS/CDS.derived-era5-land-daily-statistics.y2000.m8.daily_mean/2m_dewpoint_temperature_0_daily-mean.nc")

rr<-rotate(r)
plot(rr[[1]])


