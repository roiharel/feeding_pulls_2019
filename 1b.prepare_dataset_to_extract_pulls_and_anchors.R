library(lubridate)
library(data.table)
library(geosphere)
library(dplyr)
library(move)

#### prepare dataset to extract pulls and anchors using eobs data

#load movement dataset
load("~/GPSdata.Rdata")

a <-mp
rm("mp")
a <- plyr::compact(a) #removes null objects

for (i in 1:length(a)){
  options(digits=10)
  w <- as.data.frame(a[[i]])
  w <- w[,c("local_identifier", "location_lat", "location_long", "timestamps")]
  a[[i]] <- w
  rm(w)
} #this for loop reduces the size of the list -by only keeping the useful columns- which is massive otherwise

library (plyr)
d <- ldply (a, data.frame)
rm(a)


#exclude outliers outside the extended study area
data <- d[which(d$location_lat > 0.260464 & d$location_lat < 0.324190 & d$location_long < 36.926129 & d$location_long >36.873491),]
rm("d")

library(sp)
library(rgdal)
library(gdata)

#   transform positional data to utm values
# distance will therefore be in meters
# geological position of Mpala Kenya in UTM
projectionWant <- "+proj=utm +zone=37 +ellps=WGS84" 
coordinates(data) <- c("location_long",  "location_lat")
proj4string(data) <- CRS("+proj=longlat +datum=WGS84")

data2 <- spTransform(data, CRSobj=projectionWant)

#   variables needed
ids <- unique(data2$local_identifier)

#   positional matrix of all individuals
# converting individual position from movebank object
# into two matrices with each individual position 
#
####

# timeline organisation and matrices
days <- unique(as.Date(data2$timestamps))
hours <- seq(3,16,1)
minutes <- seq(0,59,1)
seconds <- seq(0,59,1)

times <- expand.grid(day=days,hour=hours,minute=minutes,second=seconds)
times <- paste(times$day," ",times$hour,":",times$minute,":",times$second,sep="")
times <- sort(strptime(times, format="%Y-%m-%d %H:%M:%S", tz="UTC"))
times <- as.POSIXct(times)

#save(times, file="~/times.Rdata")

rm(list=setdiff(ls(), c("times", "ids", "data2")))

xs <- ys <- matrix(NA, nrow = length(ids), ncol = length(times))

# fill matrices with data from data2 for all individuals
for (i in 1:length(ids)) {
  a <- data2[which(data2$local_identifier == ids[i]),]
  xs[i,] <- a@coords[,1][match(times, a$timestamps)]
  ys[i,] <- a@coords[,2][match(times, a$timestamps)]
}
rownames(xs) <- ids
rownames(ys) <- ids


save(xs, file="~/x_axis.Rdata")
save(ys, file="~/y_axis.Rdata")


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


