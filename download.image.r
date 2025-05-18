library(data.table)
library(rinat)
library(httr)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
if (F){
  df<-fread("../Data/Similar_phenology_invasive_code&data_20250218/iNadata.csv")
  df[is.na(url)]$url<-sprintf("https://www.inaturalist.org/observations/%d", 
                              c(6221797, 6222017, 6222099, 6222223, 6491383, 6222187, 6491403, 6789124))
  df$image<-regmatches(df$url, regexpr("\\d+$", df$url))
  saveRDS(df, "../Data/iNadata.rda")
}
df<-readRDS("../Data/iNadata.rda")
i=1
df<-df[sample(nrow(df), nrow(df))]

df<-df[grepl("static", image_url)]
for (i in c(1:nrow(df))){
  
  number <- regmatches(df[i]$url, regexpr("\\d+$", df[i]$url))
  target<-sprintf("../Data/Images/%s.jpg", number)
  if (file.exists(target)){
    next()
  }
  url<-gsub("medium", "original", df[i]$image_url)
  if (F){
    if (grepl("amazonaws", url)){
      url<-gsub("https://inaturalist-open-data.s3.amazonaws.com/",
                "https://static.inaturalist.org/", url)
    }
  }
  
  if (grepl("static.inaturalist.org", url)){
    url<-gsub("https://static.inaturalist.org/",
              "https://inaturalist-open-data.s3.amazonaws.com/", url)
  }
  #https://static.inaturalist.org/photos/93481325/original.jpg
  
  print(paste(i, nrow(df), url))
  try({
    
    resp <- GET(url, user_agent("Mozilla/5.0 (compatible; MyRBot/1.0)"))
    
    if (status_code(resp) == 200) {
      writeBin(content(resp, "raw"), target)
      Sys.sleep(runif(1, 0.5, 1.5))  # 随机暂停 0.5 ~ 1.5 秒
    } else {
      message("下载失败：", url)
    }
  }, silent = TRUE)
  
 #download.file(url, target)
}
df[208652]
