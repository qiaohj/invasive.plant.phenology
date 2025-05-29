setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")

source<-"/media/huijieqiao/WD22T_11/CDS"
zips<-list.files(source, pattern = "\\.zip")
f<-zips[1]
for (f in zips){
  print(f)
  label<-gsub("2m_temperature\\.zip", "", f)
  zip<-sprintf("%s/%s", source, f)
  unzip(zipfile = zip, exdir = source)
  file.rename(sprintf("%s/skin_reservoir_content_0_daily-mean.nc", source),
              sprintf("%s/%s%s.nc", source, label, "skin_reservoir_content"))
  file.rename(sprintf("%s/2m_temperature_0_daily-mean.nc", source),
              sprintf("%s/%s%s.nc", source, label, "2m_temperature"))
  file.rename(zip, sprintf("/media/huijieqiao/WD22T_11/CDS.zips/%s", f))
}

