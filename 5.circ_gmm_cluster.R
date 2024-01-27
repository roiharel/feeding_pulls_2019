#DOES GMM CLUSTERING TO PRODUCE GMM CLUSTERED DATA ON THE INITIATION DIRECTIONS OF CO-PULL-ANCHORERS

library(fields)
library(lme4)
library(mclust)

source("~/pull_anchor_funcs.R")

#FUNCTIONS
#distance between points
dist.func <- function(x1,x2){
	return((x1-x2)^2)
}

#get a matrix of distances for every event
event.dist.matrix <- function(x,y){
	xmat<-outer(x,x,FUN=dist.func)
	ymat<-outer(y,y,FUN=dist.func)
	out<- sqrt(xmat + ymat)
	return(out)
	
}

rotate.vector <- function(v,theta){
	x <- v$x
	y <- v$y
	v1 <- rbind(x,y)
	rot.mat <- rbind(c(cos(theta),-sin(theta)),c(sin(theta),cos(theta)))
	v2 <- rot.mat %*% v1
	return(v2)
	
} 


#read in x and y, return cluster numbers
cluster.nums.gmm <- function(ins){
	x <- ins[,1]
	y <- ins[,2]
	
	#get angles
	angs <- atan2(y,x)
	
	out <- rep(NA,length(x))
	if(nrow(ins)>2){
		#sort by angle
		angs.sorted <- order(angs)
		angs.sorted.wrap <- angs[c(angs.sorted,angs.sorted[1])]
		
		#get distances between sorted angles
		adj.ang.dists <- abs(diff(angs.sorted.wrap))
		
		#get maximum distance
		max.dist.idx <- which.max(adj.ang.dists)
		x.sorted <- cos(angs.sorted.wrap)
		y.sorted <- sin(angs.sorted.wrap)
		cut.x <- (x.sorted[max.dist.idx] + x.sorted[max.dist.idx+1])/2
		cut.y <- (y.sorted[max.dist.idx] + y.sorted[max.dist.idx+1])/2
		
		#get cutting angle
		cut.ang <- atan2(cut.y,cut.x)
		if(cut.ang < 0){
			cut.ang <- cut.ang + 2*pi
		}
		
		#get new angles (rotated so that 0 is between the most distant adjacent vectors)
		new.angs <- rotate.vector(data.frame(x=x,y=y),-cut.ang)
		new.angs <- atan2(new.angs[2,],new.angs[1,])
		
		#cluster the vectors
		out<-Mclust(new.angs)$classification
		
	}
	return(out)
	
	
	
	
}

#ins is a data frame with 5 columns - x component, y component, cluster id, and result vector (x and y)
which.cluster.wins <- function(ins){

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

n.inds.each.cluster<-function(init.clusts,max.clusts=4){
	a<-table(init.clusts)
	out<-rep(0,max.clusts)
	out[as.numeric(names(a))] <- a
	return(out)
	
}

n.inds.in.cluster<-function(ins){
	init.clusts <- ins[[1]]
	which.clust <- ins[[2]]
	out<-sum(init.clusts==which.clust)
	return(out)
}

#mean direction (takes in a matrix of two vectors (columns), x and y, which given the x and y components of direction vectors to average)
mean.direc.x <- function(vecs){
	x <- vecs[,1]
	y <- vecs[,2]
	meanx <- mean(x)
	meany <- mean(y)
	x <- meanx / sqrt(meanx^2 + meany^2)
	out<- x
	return(out)
	
}
#mean direction (takes in a matrix of two vectors (columns), x and y, which given the x and y components of direction vectors to average)
mean.direc.y <- function(vecs){
	x <- vecs[,1]
	y <- vecs[,2]
	meanx <- mean(x)
	meany <- mean(y)
	y <- meany / sqrt(meanx^2 + meany^2)
	out <- y
	return(out)
	
}


#MAIN

#takes in data grouped into events with directions included
load("~/events_with_directions_thresh35.RData")

#SET UP
#measure the level of agreement
agreement.by.event<- sapply(by(cbind(events$init.direc.x,events$init.direc.y),events$event.num,circ.var),identity)
all.events <- as.numeric(names(agreement.by.event))
events$agreement <- 1-agreement.by.event[match(events$event.num,all.events)]

#remove events with only one individual initiating
data<-events[which(events$n.inits>1),]
data<-data[order(data$event.num),]

#for each event, return which cluster each initiator belongs to
clusters <- sapply(by(data=cbind(data$init.direc.x,data$init.direc.y),INDICES=data$event.num,FUN=cluster.nums.gmm),identity)
data$init.cluster <- unlist(clusters)

#make a column with the number of clusters in each event
n.clusters <- sapply(by(data=data$init.cluster,data$event.num,FUN=max),identity)
data$n.clusters <- n.clusters[match(data$event.num,names(n.clusters))]

#numeric event type (pull = 1, anchor = 0)
data$type.num <- 0
data$type.num[which(data$type=='pull')]<-1

#compute which cluster wins for each event (the one whose mean initiation angle is closest to the resultant vector)
winning.cluster <- sapply(by(data.frame(init.cluster=data$init.cluster,init.direc.x=data$init.direc.x,init.direc.y=data$init.direc.y,result.direc.x=data$result.direc.x,result.direc.y=data$result.direc.y),data$event.num,which.cluster.wins),identity)
data$winning.cluster<-winning.cluster[match(data$event.num,names(winning.cluster))]

#number of individuals in each cluster
clust1 <- sapply(by(data=list(data$init.cluster,1),INDICES=data$event.num,FUN=n.inds.in.cluster),identity)
clust2 <- sapply(by(data=list(data$init.cluster,2),INDICES=data$event.num,FUN=n.inds.in.cluster),identity)
clust3 <- sapply(by(data=list(data$init.cluster,3),INDICES=data$event.num,FUN=n.inds.in.cluster),identity)
clust4 <- sapply(by(data=list(data$init.cluster,4),INDICES=data$event.num,FUN=n.inds.in.cluster),identity)
clust5 <- sapply(by(data=list(data$init.cluster,5),INDICES=data$event.num,FUN=n.inds.in.cluster),identity)
clust6 <- sapply(by(data=list(data$init.cluster,6),INDICES=data$event.num,FUN=n.inds.in.cluster),identity)
clust7 <- sapply(by(data=list(data$init.cluster,7),INDICES=data$event.num,FUN=n.inds.in.cluster),identity)
clust8 <- sapply(by(data=list(data$init.cluster,8),INDICES=data$event.num,FUN=n.inds.in.cluster),identity)
clust9 <- sapply(by(data=list(data$init.cluster,9),INDICES=data$event.num,FUN=n.inds.in.cluster),identity)
data$votes.clust1 <- clust1[match(data$event.num,names(clust1))]
data$votes.clust2 <- clust2[match(data$event.num,names(clust2))]
data$votes.clust3 <- clust3[match(data$event.num,names(clust3))]
data$votes.clust4 <- clust4[match(data$event.num,names(clust4))]
data$votes.clust5 <- clust4[match(data$event.num,names(clust5))]
data$votes.clust6 <- clust4[match(data$event.num,names(clust6))]
data$votes.clust7 <- clust4[match(data$event.num,names(clust7))]
data$votes.clust8 <- clust4[match(data$event.num,names(clust8))]
data$votes.clust9 <- clust4[match(data$event.num,names(clust9))]

#save output
save(data, file="~/co_events_gmm_cluster_thres35.RData")
