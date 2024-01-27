#2-INDIVIDUAL COMPETITIONS

library(fields)
library(ClassDiscovery)

#FUNCTIONS
n.votes.with <- function(df){
	votes.with <- rep(NA,nrow(df))
	for(i in 1:nrow(df)){
		my.clust <- df$init.cluster[i]
		if(!is.na(my.clust)){
			votes.with[i] <- df[i,paste('votes.clust',my.clust,sep='')]
		}
		
	}
	return(votes.with)
	
}

transitivity<-function(direc.mat){
	D<-direc.mat
	triads.tot <- 0
	triads.rightway <- 0
	N<-dim(D)[1]
	for(i in 1:N){
		for(j in 1:N){
			for(k in 1:N){
				if(!is.na(D[i,j]) & !is.na(D[i,k]) & !is.na(D[j,k])){
					if(D[i,j]>0 & D[j,k]>0){
						if(D[i,k]>0){
							triads.rightway <- triads.rightway + 1
							triads.tot <- triads.tot+1
						}
						else{
							if(D[i,k]<0){
								triads.tot <- triads.tot + 1
							}
						}
					}
				}
			}
		}
	}
	transitivity<-triads.rightway / triads.tot
	return(transitivity)
}

#shuffle a directionality matrix, preserving the antisymmetric nature of it, and preserving any NAs
shuffle.direc.mat <- function(direc.mat){
	
	N<-dim(direc.mat)[1]
	if(dim(direc.mat)[2]!=N){stop('ERROR: input matrix must be square')}
	
	#get vector of upper triangle
	vec <- direc.mat[which(upper.tri(direc.mat) & !is.na(direc.mat))]
	
	#shuffle vector
	vec.shuff <- sample(vec)
	
	new.mat <- matrix(0,nrow=N,ncol=N)
	new.mat[which(is.na(direc.mat))]<-NA
	new.mat[which(upper.tri(direc.mat) & !is.na(direc.mat))] <- vec.shuff
	new.mat <- new.mat - t(new.mat)
	diag(new.mat)<-NA
	
	return(new.mat)
	
	
	
}


#LOAD DATA

load("~/xs.Rdata")
load("~/ys.Rdata")
load("~/co-pull_co-anchor_thresh35.Rdata")
events <- mydata[which(mydata$not.cohesive.event.cluster == 0),] # only keeps events in which the leaders and followers where within 30 meters from each other

events <- mydata[which(mydata$na_prop <= 0.5),] # only keeps events in which more than 50% of the group's tags were working


source("~/pull_anchor_funcs.R")



#PREP DATA
n.inits <- aggregate(events$event.num,by=list(events$event.num),length)
events$n.inits <- n.inits$x[match(events$event.num,n.inits$Group.1)]
n.pulls <- sapply(by(events,events$event.num,FUN=function(x){sum(x$type=='pull')}),identity)
events$n.pulls <- n.pulls[match(events$event.num,names(n.pulls))]

#get events where there are 2 initiators and only 1 is followed
data <- events[which(events$n.inits==2 & events$n.pulls==1),]

#get puller characteristics (directedness and speed)
#get speed of puller from t1 to t2
data$puller.speed <- sqrt((xs[cbind(data$leader,data$t2)]-xs[cbind(data$leader,data$t1)])^2 + (ys[cbind(data$leader,data$t2)]-ys[cbind(data$leader,data$t1)])^2) / (data$t2 - data$t1)

#get directedness of puller from t1 to t2
data$puller.path.direc <- NA
for(i in 1:nrow(data)){
	curr <- data[i,]
	dx <- diff(xs[curr$leader,curr$t1:curr$t2])
	dy <- diff(ys[curr$leader,curr$t1:curr$t2])
	ds <- sqrt(dx^2 + dy^2)
	path.len <- sum(ds)
	disp <- sqrt((xs[curr$leader,curr$t2]-xs[curr$leader,curr$t1])^2 + (ys[curr$leader,curr$t2]-ys[curr$leader,curr$t1])^2)
	data$puller.path.direc[i] <- disp / path.len
	
}

#PROBABILITY OF FIRST INDIVIDUAL WINNING VS THE PATH DIREC OF EACH PULLER
path.direc.bins <- seq(0,1,.1)
event.list <- unique(data$event.num)
row.wins <- row.loses <- histo <- array(0,dim=c(length(path.direc.bins),length(path.direc.bins)))
corrects <- wrongs <- 0
for(i in 1:length(event.list)){
	curr <- data[which(data$event.num==event.list[i]),]
	curr <- curr[sample(nrow(curr)),]
	pd1 <- curr$puller.path.direc[1]
	pd2 <- curr$puller.path.direc[2]
	bin1 <- which(path.direc.bins > pd1)[1]
	bin2 <- which(path.direc.bins > pd2)[1]
	histo[bin1,bin2] <- histo[bin1,bin2] + 1
	if(curr$type[1]=='pull'){
		row.wins[bin1,bin2] <- row.wins[bin1,bin2] + 1
		if(pd1 > pd2){
			corrects <- corrects + 1
		}
		else{
			wrongs <- wrongs + 1
		}
	}
	if(curr$type[2]=='pull'){
		row.loses[bin1,bin2] <- row.loses[bin1,bin2] + 1
		if(pd2 > pd1){
			corrects <- corrects + 1
		}
		else{
			wrongs <- wrongs + 1
		}
	}
}

#PROBABILITY OF FIRST INDIVIDUAL WINNING VS THE SPEED OF EACH PULLER
log.speed.bins <- seq(-4,3,.5)
event.list <- unique(data$event.num)
row.wins <- row.loses <- histo <- array(0,dim=c(length(log.speed.bins),length(log.speed.bins)))
corrects <- wrongs <- 0
for(i in 1:length(event.list)){
	curr <- data[which(data$event.num==event.list[i]),]
	curr <- curr[sample(nrow(curr)),]
	speed1 <- log(curr$puller.speed[1])
	speed2 <- log(curr$puller.speed[2])
	bin1 <- which(log.speed.bins > speed1)[1]
	bin2 <- which(log.speed.bins > speed2)[1]
	histo[bin1,bin2] <- histo[bin1,bin2] + 1
	if(curr$type[1]=='pull'){
		row.wins[bin1,bin2] <- row.wins[bin1,bin2] + 1
		if(speed1 < speed2){
			corrects <- corrects + 1
		}
		else{
			wrongs <- wrongs + 1
		}
	}
	if(curr$type[2]=='pull'){
		row.loses[bin1,bin2] <- row.loses[bin1,bin2] + 1
		if(speed2 < speed1){
			corrects <- corrects + 1
		}
		else{
			wrongs <- wrongs + 1
		}
	}
}

#COMBINE PATH DIRECTEDNESS AND SPEED
data$speed.direc.combo <- -log(data$puller.speed)*data$puller.path.direc
combo.bins <- seq(-2,2,.2)
event.list <- unique(data$event.num)
row.wins <- row.loses <- histo <- array(0,dim=c(length(combo.bins),length(combo.bins)))
corrects <- wrongs <- 0
for(i in 1:length(event.list)){
	curr <- data[which(data$event.num==event.list[i]),]
	curr <- curr[sample(nrow(curr)),]
	combo1 <- curr$speed.direc.combo[1]
	combo2 <- curr$speed.direc.combo[2]
	bin1 <- which(combo.bins > combo1)[1]
	bin2 <- which(combo.bins > combo2)[1]
	histo[bin1,bin2] <- histo[bin1,bin2] + 1
	if(curr$type[1]=='pull'){
		row.wins[bin1,bin2] <- row.wins[bin1,bin2] + 1
		if(combo1 > combo2){
			corrects <- corrects + 1
		}
		else{
			wrongs <- wrongs + 1
		}
	}
	if(curr$type[2]=='pull'){
		row.loses[bin1,bin2] <- row.loses[bin1,bin2] + 1
		if(combo2 > combo1){
			corrects <- corrects + 1
		}
		else{
			wrongs <- wrongs + 1
		}
	}
}




#WINNING AND LOSING WHEN THERE ARE 2 INITIATORS

#compute directionality matrix based on these

N<-max(max(events$leader), max(events$follower))

event.mat <- matrix(0,nrow=N,ncol=N)
all.events <- unique(data$event.num)
winners <- losers <- vector()

for(i in 1:length(all.events)){
  curr <- data[which(data$event.num==all.events[i]),]
  winner <- curr$leader[which(curr$type=='pull')]
  loser <- curr$leader[which(curr$type=='anchor')]
  event.mat[winner,loser]<-event.mat[winner,loser]+1
  winners <- c(winners,winner)
  losers <- c(losers,loser)
}

prob.win.2pullers <- event.mat / (event.mat + t(event.mat))

prob.win.2pullers.tot <- rowSums(event.mat)/(rowSums(event.mat)+rowSums(t(event.mat)))
data <- events[which(events$n.inits==1),]
p.succ.alone <- rep(NA,N)
p.follow.alone <- rep(NA,N)
for(i in 1:N){
  p.succ.alone[i] <- length(which(data$leader==i & data$type=='pull')) / length(which(data$leader==i))
  p.follow.alone[i] <- length(which(data$follower==i & data$type=='pull')) / length(which(data$follower==i))
}
p.succ.dyad <- matrix(NA,nrow=N,ncol=N)
for(i in 1:N){
  for(j in 1:N){
    p.succ.dyad[i,j]<-length(which(data$leader==i & data$follower == j & data$type=='pull')) / length(which(data$leader==i & data$follower==j))
  }
}

#save(list=c('p.succ.alone','p.succ.dyad','prob.win.2pullers'),file="~/leadership_efficiency_thres35.RData")

