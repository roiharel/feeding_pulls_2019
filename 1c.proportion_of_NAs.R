library(lubridate)
library(tidyverse)
library(viridis)

#load("~/times.Rdata")
#load("~/x_axis.Rdata")

load("~/times.Rdata")
load("~/xs.Rdata")

metadata <- read.csv("~/metadata_Avulturinum.csv")#loading metadata to know when tags were on birds
metadata <- metadata[which(metadata$Animal.ID %in% rownames(xs)),]
metadata$Deploy.off.timestamp <- as.Date(as.character(metadata$Deploy.off.timestamp), format = "%d.%m.%Y")
metadata$Deploy.on.timestamp <- as.Date(as.character(metadata$Deploy.on.timestamp), format = "%d.%m.%Y")
metadata$Deploy.off.timestamp[which(is.na(metadata$Deploy.off.timestamp))] <- as.Date(Sys.time())

ids <- rownames(xs)

tracked <- data.frame(days = unique(as.Date(times)))

for (i in 1:length(ids)){ # this for loop will give time sequences when each bird was tracked
  metadata_ <- metadata[which(metadata$Animal.ID == ids[i]),]
  
  if(length(metadata_$Deploy.on.timestamp) == 1){
    seqd <- seq(metadata_$Deploy.on.timestamp, metadata_$Deploy.off.timestamp,1)}
  
  if(length(metadata_$Deploy.on.timestamp) == 2){
    seqd <- c(seq(metadata_$Deploy.on.timestamp[1], metadata_$Deploy.off.timestamp[1],1),  seq(metadata_$Deploy.on.timestamp[2], metadata_$Deploy.off.timestamp[2],1))}
  
  if(length(metadata_$Deploy.on.timestamp) == 3){
    seqd <- c(seq(metadata_$Deploy.on.timestamp[1], metadata_$Deploy.off.timestamp[1],1),  seq(metadata_$Deploy.on.timestamp[2], metadata_$Deploy.off.timestamp[2],1), seq(metadata_$Deploy.on.timestamp[3], metadata_$Deploy.off.timestamp[3],1))}
  
  tracked[,ncol(tracked) + 1] <- as.numeric(tracked$days %in% seqd)
}
sums <- rowSums(tracked[,c(2:ncol(tracked))])# this gives the number of birds tracked each of the study days

n.tracked <- sums[match(as.Date(times), tracked$days)]# this gives the number of birds tracked in each of the study timestamps
nonas <- colSums(is.na(xs))# this gives the number of NAs in each timestamp

a <- (nonas-nrow(xs)+n.tracked)/n.tracked

nas_prop <- data.frame(t = 1:length(times), na_prop = a)
save(nas_prop, file="~/nas_prop.Rdata")
