library(maptools)

load("~/co-pull_co-anchor_thresh35.Rdata")
mydata <- mydata[order(mydata$event.num),]

load("~/xs.Rdata")
load("~/ys.Rdata")
load("~/times.Rdata")

not.cohesive.event.cluster <- c()
median_time_t2 <- c()

for (i in sort(unique(mydata$event.num))){
  data_ <- mydata[which(mydata$event.num == i), ] #data_ includes data of just one event
  a <- sort(unique(c(data_$leader, data_$follower))) #which are the individuals being leaders and followers in this event ?
  
  centx <- mean(xs[a, median(data_$t2)], na.rm=T) # centroid of each event
  centy <- mean(ys[a, median(data_$t2)], na.rm=T)
  
  #the code below adds a column on whether the cluster of leaders and followers of each event is cohesive or not- each leader or follower has to be within 30 meters from the centroid of these leaders and followers
  
  tempx <- xs[,median(data_$t2)]  
  tempy <- ys[,median(data_$t2)]  
  
  #dist <- sqrt((centx-tempx)^2 + (centy-tempy)^2) it's the same as below
  dist <- as.matrix(dist(cbind(c(centx,tempx),c(centy,tempy))))[1,][-1]
  
  
  median_time_t2 <- c(median_time_t2, rep(median(data_$t2), length(data_$leader)))
  
  w <- as.numeric(any(dist[a] > 30))   #Binary variable: cohesive cluster of leaders and followers? 
  
  not.cohesive.event.cluster <- c(not.cohesive.event.cluster, rep(w, length(data_$leader)))
  if(i %in% seq(3, length(mydata$leader), 1000)){
    print(paste("events", i))}
  
}

mydata$not.cohesive.event.cluster <- not.cohesive.event.cluster 
mydata$median_time_t2 <- median_time_t2 

  # the code herein will add one column on the dataset on the proportion of tracked individuals that have NA on t2, later we will only keep events where this proportion is lower than 0.5

load("~/nas_prop.Rdata")
mydata$na_prop <- nas_prop$na_prop[match(mydata$t2, nas_prop$t)]

save(mydata, file="~/co-pull_co-anchor_thresh35.RData")
