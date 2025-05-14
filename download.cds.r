library(ecmwfr)

year<-"2000"
month<-"8"

for (year in c(2010:2024)){
  for (month in c(1:12)){
    request <- list(
      dataset_short_name = "derived-era5-land-daily-statistics",
      #product_type   = "reanalysis",
      format = "grib",
      variable = c("2m_dewpoint_temperature",
                   "2m_temperature",
                   "skin_temperature",
                   "soil_temperature_level_1",
                   "soil_temperature_level_2",
                   "soil_temperature_level_3",
                   "soil_temperature_level_4"),
      year = as.character(year),
      month = as.character(month),
      day = c("01", "02", "03",
              "04", "05", "06",
              "07", "08", "09",
              "10", "11", "12",
              "13", "14", "15",
              "16", "17", "18",
              "19", "20", "21",
              "22", "23", "24",
              "25", "26", "27",
              "28", "29", "30",
              "31"),
      #time = c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"),
      daily_statistic= c("daily_mean"), #, "daily_minimum", "daily_maximum"
      time_zone= "utc+00:00",
      frequency= "6_hourly",
      target = sprintf("cds.y%s.m%s.grib", year, month),
      download_format= "zip"
    )
    
    wf_request(request, user =user)
  }
}