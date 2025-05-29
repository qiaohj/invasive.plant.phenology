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

setwd("/media/huijieqiao/WD22T_11/invasive.plant.phenology/invasive.plant.phenology")
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