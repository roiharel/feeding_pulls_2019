library(lubridate)
library(data.table)
library(dplyr)
library(move)
library(readr)

library(geosphere)
library(sf)
library(rgdal)
library(gdata)

#### prepare dataset to extract pulls and anchors using eobs data

#load movement dataset
#load("~/GPSdata.Rdata")
setwd("~/GitHub/feeding_pulls_2019")
load("~/GitHub/feeding_pulls_2019/data/baboon_data_2019.RData")

gps_v0 <- read_csv("data/Papio Anubis Mpala 2019_gps_raw_stage0.csv", 
                   col_types = cols(`event-id` = col_skip(), 
                                    visible = col_skip(), timestamp = col_datetime(format = "%Y-%m-%d %H:%M:%S"), 
                                    `bar:barometric-pressure` = col_skip(), 
                                    `data-decoding-software` = col_skip(), 
                                    `eobs:activity` = col_skip(), `eobs:activity-samples` = col_skip(), 
                                    `eobs:battery-voltage` = col_skip(), 
                                    `eobs:fix-battery-voltage` = col_skip(), 
                                    `eobs:key-bin-checksum` = col_skip(), 
                                    `mag:magnetic-field-raw-x` = col_skip(), 
                                    `mag:magnetic-field-raw-y` = col_skip(), 
                                    `mag:magnetic-field-raw-z` = col_skip(), 
                                    `quaternion-raw-w` = col_skip(), 
                                    `quaternion-raw-x` = col_skip(), 
                                    `quaternion-raw-y` = col_skip(), 
                                    `quaternion-raw-z` = col_skip(), 
                                    `sensor-type` = col_skip(), `individual-taxon-canonical-name` = col_skip(), 
                                    `study-name` = col_skip()))
View(gps_v0)
gps_v0 <- as.data.frame(gps_v0)
gps_v0 <- plyr::compact(gps_v0) #removes null objects

gps_v0 <- gps_v0 %>%
  select(`tag-local-identifier`, `location-lat`, `location-long`, `timestamp`)

gps_v0 <- gps_v0 %>%
  filter(!if_any(everything(), is.na))


#exclude outliers outside the extended study area
#gps_v0 <- gps_v0[which(gps_v0$'location-lat' > 0.260464 & gps_v0$'location-lat' < 0.324190 & gps_v0$'location-long' < 36.926129 & gps_v0$'location-long' >36.873491),]

#   transform positional data to utm values
# distance will therefore be in meters
# geological position of Mpala Kenya in UTM
projectionWant <- "+proj=utm +zone=37 +ellps=WGS84" 
coordinates(gps_v0) <- c("location-long",  "location-lat")
proj4string(gps_v0) <- CRS("+proj=longlat +datum=WGS84")

gps_v0 <- spTransform(gps_v0, CRSobj=projectionWant)

#   variables needed
ids <- unique(gps_v0$'tag-local-identifier')

days <- unique(as.Date(gps_v0$timestamp))
hours <- seq(3,16,1)
minutes <- seq(0,59,1)
seconds <- seq(0,59,1)

times <- expand.grid(day=days,hour=hours,minute=minutes,second=seconds)
times <- paste(times$day," ",times$hour,":",times$minute,":",times$second,sep="")
times <- sort(strptime(times, format="%Y-%m-%d %H:%M:%S", tz="UTC"))
times <- as.POSIXct(times)


#   positional matrix of all individuals
# converting individual position from movebank object
# into two matrices with each individual position 
#
xs <- ys <- matrix(NA, nrow = length(ids), ncol = length(times))

# fill matrices with data from data2 for all individuals
for (i in 1:length(ids)) {
  a <- gps_v0[which(gps_v0$'tag-local-identifier' == ids[i]),]
  xs[i,] <- a@coords[,1][match(times, a$timestamp)]
  ys[i,] <- a@coords[,2][match(times, a$timestamp)]
}
rownames(xs) <- ids
rownames(ys) <- ids







# timeline organisation and matrices

save(times, file="~/times.Rdata")

save(xs, file="~/x_axis.Rdata")
save(ys, file="~/y_axis.Rdata")

rm(list=setdiff(ls(), c("times", "ids", "data2")))


#first timestamp of each day
a <- unique(as.Date(times))
start.day <- c() 

for (i in 1:length(a)){
  b <- times[as.Date(times) %in% a[i]]
  start.day[i] <- as.character(b[1])
}

start.day <- as.POSIXct(start.day, origin = "1970-01-01", tz= "UTC")

day.start.idxs <- which(times %in% start.day)
save(day.start.idxs, file="~/day.start.idxs.Rdata")


