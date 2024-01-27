#library(moments)
library(fields)
library(PerformanceAnalytics)
library(diptest)
library(grDevices)
library(openxlsx)

#FUNCTIONS
angle.btwn.vecs <- function(vs){
  x1 <- vs$x[1]
  x2 <- vs$x[2]
  y1 <- vs$y[1]
  y2 <- vs$y[2]
  
  ang <- acos((x1*x2 + y1*y2) / (sqrt(x1^2+y1^2) * sqrt(x2^2+y2^2)))
  return(ang)
  
}

rotate.vector <- function(v,theta){
  x <- v$x
  y <- v$y
  v1 <- rbind(x,y)
  rot.mat <- rbind(c(cos(theta),-sin(theta)),c(sin(theta),cos(theta)))
  v2 <- rot.mat %*% v1
  return(v2)
  
} 

#rotate the direction vectors such that the x axis aligns with the first individual's initiation direction
#flip the vectors such that the second initiator direction is at an angle between 0 and pi
compute.rotated.vectors <- function(df){
  x.ref <- df$init.direc.x[1]
  y.ref <- df$init.direc.y[1]
  
  theta <- -atan2(y.ref,x.ref)
  
  v1 <- rotate.vector(data.frame(x=df$init.direc.x[1],y=df$init.direc.y[1]),theta)
  v2 <- rotate.vector(data.frame(x=df$init.direc.x[2],y=df$init.direc.y[2]),theta)
  v3 <- rotate.vector(data.frame(x=df$result.direc.x[1],y=df$result.direc.y[1]),theta)
  
  df$init.direc.x.rot <- c(v1[1],v2[1])
  df$init.direc.y.rot <- c(v1[2],v2[2])
  df$result.direc.x.rot <- c(v3[1],v3[1])
  df$result.direc.y.rot <- c(v3[2],v3[2])
  
  if(df$init.direc.y.rot[2] < 0){
    df$init.direc.y.rot[2] <- -df$init.direc.y.rot[2]
    df$result.direc.y.rot <- -df$result.direc.y.rot
    
    
  }
  
  
  return(df)
  
}

#ins is a data frame with 5 columns - x component, y component, cluster id, and result vector (x and y)
which.cluster.wins.by.direction <- function(ins){
  
  if(length(unique(ins$init.cluster))==1){
    return(ins$init.cluster[1])
  } else{
    
    #mean vector of each cluster
    means.x <- sapply(by(ins$init.direc.x,ins$init.cluster,mean),identity)
    means.y <- sapply(by(ins$init.direc.y,ins$init.cluster,mean),identity)
    
    #resulting movement vector
    res.x <- ins$result.direc.x[1]
    res.y <- ins$result.direc.y[1]
    
    #difference between result vector and mean of each cluster
    diffs <- acos(means.x*res.x + means.y*res.y)
    
    #the winning cluster is the one with the minimum distance between its mean and the resulting direction
    winner <- names(means.x)[which.min(diffs)]
    return(as.numeric(unlist(winner)))
    
  }
}
  
  #LOAD DATA
  
  load("~/co_events_gmm_cluster_thres35.RData")
  
  load("~/ys.Rdata")
  load("~/xs.Rdata")
  load("~/leadership_efficiency_thres35.RData")
  pull.ranks <-NA
  
  #SOURCE FUNCTIONS
  source("~/pull_anchor_funcs.R")
  
  #1. AVERAGE VS CHOOSE - 2 individual pullers
  
  sort_by_pull_rank <- FALSE
  sort_by_pull_efficiency <- TRUE
  
  #get cases with 2 pullers where at least 1 is a pull
  events <- data[which(data$n.inits==2 & data$success==1),]
  events <- events[which(!is.na(events$result.direc.x)),]
  events <- events[sample(1:nrow(events)),]
  
  #get angle difference between the two pullers
  d.theta <- sapply(by(data.frame(x=events$init.direc.x,y=events$init.direc.y),events$event.num,FUN=angle.btwn.vecs),identity)
  events$d.theta <- d.theta[match(events$event.num,as.numeric(names(d.theta)))]
  
  #rotate both vectors by the angle between vector 1 and the x axis
  #also flip the vectors such that the second vector is between 0 and pi (and flip the resultant vector too, accordingly)
  new.events <- data.frame(NULL)
  if(sort_by_pull_efficiency){
    load("~/leadership_efficiency_thres35.RData")
  }
  for(i in unique(events$event.num)){
    curr <- events[which(events$event.num==i),]
    if(sort_by_pull_rank){
      if(pull.ranks[curr$leader[1]] < pull.ranks[curr$leader[2]]){
        curr <- curr[c(2,1),]
      }
    }
    if(sort_by_pull_efficiency){
      if(p.succ.alone[curr$leader[1]] < p.succ.alone[curr$leader[2]]){
        curr <- curr[c(2,1),]
      }
    }
    curr <- compute.rotated.vectors(curr)
    new.events <- rbind(new.events, curr)
    
  }
  events <- new.events
  
  events$degrees <- events$d.theta*(180/pi)
  
  events <- events[which(events$degrees > 130 ),] # change thresholds according to group 
  
  df <- data.frame(ids = rownames(xs), n = 1:length(rownames(xs)))
  df$ids <- as.character(df$ids)
  df$n <- as.character(df$n)
  events$leader <- as.character(events$leader)
  
  events$leader_id <- df$ids[match(events$leader, df$n)]
  
  ################this code checks how many pull-pull cases there are 
n.events <- unique(events$event.num)
n.pull.pull <- c()
exclude <- c()
for (i in 1:length(n.events)){
  events_ <- events[which(events$event.num == n.events[i]),]
  if(sum(events_$type == "pull") > 1) {n.pull.pull <- c(n.pull.pull,1)
  exclude <- c(exclude, n.events[i])}
}
 sum(n.pull.pull)
 ################ 
 '%!in%' <- function(x,y)!('%in%'(x,y))
 events <- events[which(events$event.num %!in% exclude),]
 
 sex <- read.xlsx("~/VGF_habituated_groups.xlsx", sheet=1) # change according to group
 cap <- read.xlsx("~/VGF_capture.xlsx", sheet=1)
 sex$id <- NA
 sex$id <- cap$Original.Ring.number[match(sex$Ind, cap$Colour.Bands)]
 sex$id[is.na(sex$id)] <- cap$Original.Ring.number[match(sex$Ind[is.na(sex$id)], cap$Wing.Tag)]
 
events$sex <- sex$Sex[match(events$leader_id, sex$id)]

include <- c()
for (i in 1:length(n.events)){
  events_ <- events[which(events$event.num == n.events[i]),]
  if(sum(events_$sex == "M") == 1) {include <- c(include,n.events[i])}
}

events_fm <- events[which(events$event.num %in% include),]
n.males <- length(events_fm$event.num[events_fm$sex == "M" & events_fm$type == "pull"])   
n.females <- length(events_fm$event.num[events_fm$sex == "F" & events_fm$type == "pull"])   
n.events_fm <- length(unique(events_fm$event.num))

n.males/n.events_fm

################We will now randomize sexes to check if males are more successful pullers than expected by chance
#library(gdata)
#keep(events_fm, sure = TRUE)

randomized.props <- c()

for (i in 1:1000){
  randomized <- events_fm
  randomized$sex <- sample(events_fm$sex)
  n.males.r <- length(randomized$event.num[randomized$sex == "M" & randomized$type == "pull"])   
  n.events_fm.r <- length(unique(randomized$event.num))
  randomized.props <- c(randomized.props, n.males.r/n.events_fm.r)
}
mean(randomized.props)

plot(randomized.props)
points(x= 500, y=n.males/n.events_fm, col="blue", pch = 16, cex=3)
points(x= 500, y=n.females/n.events_fm, col="red", pch = 16, cex = 3)

#calculate p values
length(which(randomized.props > n.males/n.events_fm))/1000









  