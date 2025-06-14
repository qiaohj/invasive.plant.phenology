library(data.table)
library(ggplot2)
library(patchwork)
library(sf)
setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
paired.records<-readRDS(sprintf("../Data/paired.records.%dkm.rda", dist))

#The figure for N of records per hemisphere by days of year
dfff.e<-paired.records.filtered[, c("e.days_in_year", "e.decimalLatitude",
                                    "e.phenology", "growthform")]
colnames(dfff.e)<-c("days_in_year", "decimalLatitude",
                    "phenology", "growthform")
dfff.e$type<-"exotic"
dfff.n<-paired.records.filtered[, c("n.days_in_year", "n.decimalLatitude",
                                    "n.phenology", "growthform")]
colnames(dfff.n)<-c("days_in_year", "decimalLatitude",
                    "phenology", "growthform")
dfff.n$type<-"native"

dfff<-rbindlist(list(unique(dfff.n), unique(dfff.e)))

dfff1<-dfff[phenology==12]
dfff1.1<-dfff1
dfff1.1$phenology<-1
dfff1.2<-dfff1
dfff1.2$phenology<-2

dfff2<-dfff[phenology!=12]
dfff<-rbindlist(list(dfff1.1, dfff1.2, dfff2))
dfff$hemi<-ifelse(dfff$decimalLatitude>0, "N", "S")
dfff$month<-dfff$days_in_year/30+1
dfff$phenology.l<-ifelse(dfff$phenology==1, "Flowering", "Fruiting")
saveRDS(dfff, "../Figures/Fig.X.Raw.Records.Hist/Fig.X.Raw.Records.Hist.rda")

