#This code investigates the situation when there is more than one cluster and makes plots of:

library(fields)
library(lme4)
library(Hmisc)

source("~/pull_anchor_funcs.R")


#FUNCTIONS
#distance between points
dist.func <- function(x1,x2){
	return((x1-x2)^2)
}

angle.btwn.vecs <- function(vs){
	x1 <- vs$x[1]
	x2 <- vs$x[2]
	y1 <- vs$y[1]
	y2 <- vs$y[2]
	
	ang <- acos((x1*x2 + y1*y2) / (sqrt(x1^2+y1^2) * sqrt(x2^2+y2^2)))
	return(ang)
	
}

#get a matrix of distances for every event
event.dist.matrix <- function(x,y){
	xmat<-outer(x,x,FUN=dist.func)
	ymat<-outer(y,y,FUN=dist.func)
	out<- sqrt(xmat + ymat)
	return(out)
	
}

#cluster vectors
cluster.vecs <- function(dist.mat,cut.h=sqrt(2)){
	dist.mat <- as.dist(dist.mat,diag=FALSE)
	a<-hclust(dist.mat)
	clusts<-cutree(a,h=cut.h)
	return(clusts)
	
}

#read in x and y, return cluster numbers
cluster.nums <- function(ins,cut.h=sqrt(2)){
	x <- ins[,1]
	y <- ins[,2]
	if(nrow(ins)>2){
		mat<-event.dist.matrix(x,y)
		clusts<-cluster.vecs(mat,cut.h)
		return(clusts)
	}
	else{
		if(sqrt(dist.func(x[1],x[2])+dist.func(y[1],y[2])) > cut.h){
			return(c(1,2))
		}
		else{
			return(c(1,1))
		}
		
		
		
	}
	
}

which.cluster.wins.old <- function(ins){
	clusters <- ins[,1]
	types <- ins[,2]
	wins <- clusters[which(types==1)]
	if(length(wins)==0){
		return(NA)
	}
	else{
		a<-table(types,clusters)
		return(as.numeric(colnames(a)[which.max(a[which(rownames(a)==1),])]))
	}
}

which.cluster.wins <- function(ins){
	clusters <- ins[,1]
	types <- ins[,2]
	wins <- clusters[which(types==1)]
	if(length(wins)==0){
		return(NA)
	}
	if(length(unique(wins))==1){
		return(wins[1])
	}
	else{
		return(NaN)
	}
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

#load data

load("~/co_events_gmm_cluster_thres35.RData")


#randomise the raws of the dataset
set.seed(42)
rows <- sample(nrow(data))
data <- data[rows,]

first <- min(min(data$leader), min(data$follower))

NB <- max(max(data$leader), max(data$follower))

#number of clusters with pulls in each event (clustering accuracy)  ~93%
pulls<-data[which(data$type=='pull' & data$n.clusters>1),]
n.clust.in.pulls<-sapply(by(data=pulls$init.cluster,INDICES=pulls$event.num,FUN=function(x){return(length(unique(x)))}),identity)
n.clusters.with.pulls<- table(n.clust.in.pulls)/sum(table(n.clust.in.pulls))

#probability of success vs. number of clusters and number of initiators
ind.bins <- 1:NB
n.clusts.bins <- 1:4
p.success <- array(NA,dim=c(length(ind.bins),length(n.clusts.bins)))
n.inds <- sapply(by(data=data$n.inits,data$event.num,FUN=max),identity)
n.clusts <- sapply(by(data=data$init.cluster,data$event.num,FUN=max),identity)
success <- sapply(by(data=data$success,data$event.num,FUN=max),identity)
for(i in ind.bins){
	for(k in 1:length(n.clusts.bins)){
		j<-rev(n.clusts.bins)[k]
		r <- which(n.inds == i & n.clusts == j)
		if (length(r) > 0){
		p.success[i,k] <- mean(success[which(n.inds == i & n.clusts == j)], na.rm=T)
	  }
	}
}

#THE 2-CLUSTER CASE


#probability that cluster 1 wins as a function of number of inds in cluster 1 - number of inds in cluster 2
data.2clust <- data[which(data$n.clusters==2),]
vote.diff <- data.2clust$votes.clust1 - data.2clust$votes.clust2
clust1wins <- data.2clust$winning.cluster == 1
vote.diff <- sapply(by(vote.diff,INDICES=data.2clust$event.num,max),identity)
clust1wins <- sapply(by(clust1wins,data.2clust$event.num,max),identity)
a<-table(clust1wins,vote.diff)
sums <- colSums(a)
vote.diffs <- sort(unique(vote.diff))
prob.1.wins <- a[2,]/sums

#plot(vote.diffs,prob.1.wins,xlab='votes for cluster 1 - votes for cluster 2',ylab='probability of cluster 1 being chosen',pch=19)

#fit a sigmoid
df<- data.frame(ninds=vote.diffs,prob=prob.1.wins)
params <- list(Asym=1, xmid =0.5, scal=0.5)
model<-nls(prob ~ SSlogis(ninds,Asym,xmid,scal) ,data=df, start = params, trace = F, control = list(maxiter = 500))
newx <- vote.diffs
newy <- predict(model,newdata=data.frame(ninds=newx))
#lines(newx,newy/max(newy),col="red",lwd=2)

#THE N-CLUSTER CASE (this only works for 2 clusters at the moment)

#get data
n.clusters <- 2
foc.cluster<-1
min.inits <- 0
max.inits <- 40
min.angle <- 0*pi/180
data.nclust <- data[which(data$n.clusters==n.clusters & data$success==1 & data$n.inits <= max.inits  & data$n.inits >= min.inits),]
if(foc.cluster==1){vote.diff <- data.nclust$votes.clust1 - data.nclust$votes.clust2 - data.nclust$votes.clust3 - data.nclust$votes.clust4}
if(foc.cluster==2){vote.diff <- data.nclust$votes.clust2 - data.nclust$votes.clust1 - data.nclust$votes.clust3 - data.nclust$votes.clust4}
if(foc.cluster==3){vote.diff <- data.nclust$votes.clust3 - data.nclust$votes.clust1 - data.nclust$votes.clust2 - data.nclust$votes.clust4}
if(foc.cluster==4){vote.diff <- data.nclust$votes.clust4 - data.nclust$votes.clust1 - data.nclust$votes.clust2 - data.nclust$votes.clust3}
foc.clust.wins <- data.nclust$winning.cluster == foc.cluster
vote.diff <- sapply(by(vote.diff,INDICES=data.nclust$event.num,max),identity)
foc.clust.wins <- sapply(by(foc.clust.wins,data.nclust$event.num,max),identity)
a<-table(foc.clust.wins,vote.diff)
sums <- colSums(a)
vote.diffs <- sort(unique(vote.diff))
prob.foc.clust.wins <- a[2,]/sums

#get cluster-level data
clust.level.data <- data.frame(event.num=aggregate(data.nclust$event.num,by=list(data.nclust$event.num,data.nclust$init.cluster),max)$x)
clust.level.data$init.cluster <- aggregate(data.nclust$init.cluster,by=list(data.nclust$event.num,data.nclust$init.cluster),max)$x
clust.level.data$init.direc.x <- aggregate(data.nclust$init.direc.x,by=list(data.nclust$event.num,data.nclust$init.cluster),mean)$x
clust.level.data$init.direc.y <- aggregate(data.nclust$init.direc.y,by=list(data.nclust$event.num,data.nclust$init.cluster),mean)$x
norm <- sqrt(clust.level.data$init.direc.x^2 + clust.level.data$init.direc.y^2)
clust.level.data$init.direc.x <- clust.level.data$init.direc.x/norm
clust.level.data$init.direc.y <- clust.level.data$init.direc.y/norm
clust.level.data$reuslt.direc.x <- aggregate(data.nclust$result.direc.x,by=list(data.nclust$event.num,data.nclust$init.cluster),mean)$x
clust.level.data$result.direc.y <- aggregate(data.nclust$result.direc.y,by=list(data.nclust$event.num,data.nclust$init.cluster),mean)$x
norm <- sqrt(clust.level.data$reuslt.direc.x^2 + clust.level.data$result.direc.y^2)
clust.level.data$result.direc.x <- clust.level.data$reuslt.direc.x/norm
clust.level.data$reuslt.direc.y <- clust.level.data$result.direc.y/norm
clust.level.data$votes.clust1 <- aggregate(data.nclust$votes.clust1,by=list(data.nclust$event.num,data.nclust$init.cluster),max)$x
clust.level.data$votes.clust2 <- aggregate(data.nclust$votes.clust2,by=list(data.nclust$event.num,data.nclust$init.cluster),max)$x
d.theta <- sapply(by(data.frame(x=clust.level.data$init.direc.x,y=clust.level.data$init.direc.y),clust.level.data$event.num,FUN=angle.btwn.vecs),identity)
clust.level.data$d.theta <- d.theta[match(clust.level.data$event.num,names(d.theta))]
clust.level.data$event.clust <- paste(clust.level.data$event.num,'_',clust.level.data$init.cluster,sep='')
clust.level.data$winning.clust <- aggregate(as.numeric(data.nclust$winning.cluster),by=list(data.nclust$event.num,data.nclust$init.cluster),max)$x
clust.level.data$n.inits <- aggregate(data.nclust$n.inits,by=list(data.nclust$event.num,data.nclust$init.cluster),max)$x
clust.level.data$mean.t2 <- aggregate(data.nclust$t2,by=list(data.nclust$event.num,data.nclust$init.cluster),mean)$x
clust.level.data$min.pull.t2 <- aggregate(data.nclust$min.pull.t2,by=list(data.nclust$event.num,data.nclust$init.cluster),min)$x
clust.level.data$max.pull.t3 <- aggregate(data.nclust$max.pull.t3,by=list(data.nclust$event.num,data.nclust$init.cluster),min)$x

clust.level.data <- clust.level.data[which(clust.level.data$d.theta>=min.angle),]

vote.diffs <- bins <- seq(-NB,NB)
prob.foc.clust.wins <- rep(NA,length(bins))
sums <- rep(NA,length(bins))
for(i in 1:length(bins)){
	if(foc.cluster==1){
		curr <- clust.level.data[which((clust.level.data$votes.clust1 - clust.level.data$votes.clust2)==bins[i]),]
	}
	if(foc.cluster==2){
		curr <- clust.level.data[which((clust.level.data$votes.clust2 - clust.level.data$votes.clust1)==bins[i]),]
	}
	wins1 <- length(unique(curr$event.num[which(curr$winning.clust==1)]))
	wins2 <- length(unique(curr$event.num[which(curr$winning.clust==2)]))
	if(foc.cluster==1){
		prob.foc.clust.wins[i] <- wins1/(wins1+wins2)
	}	
	if(foc.cluster==2){
		prob.foc.clust.wins[i] <- wins2/(wins1+wins2)
	}	
	sums[i] <- wins1+wins2	
}

#bootstrap error bars
n.boots <- 1000
boot.probs <- array(NA,dim=c(length(vote.diffs),n.boots))
for(r in 1:n.boots){
	for(i in 1:length(vote.diffs)){
		boot.probs[i,r] <- mean(runif(sums[i]) <= prob.foc.clust.wins[i],na.rm=T)
	}
}
ci025 <- ci975 <- rep(NA,length(vote.diffs))
for(i in 1:length(vote.diffs)){
	ci025[i] <- quantile(boot.probs[i,],0.025,na.rm=T)
	ci975[i] <- quantile(boot.probs[i,],0.975,na.rm=T)
}

#make plot again with error bars
#quartz(width=8,height=7.5)
errbar(vote.diffs,prob.foc.clust.wins,
       ylab='Probability of choosing subgroup 1',
       xlab='Numerical difference (size of subgroup 1 minus subgroup 2)',
       pch=19,
       #xlim=c(-10,10),
       ylim=c(0,1),axes=FALSE,yplus = ci975,yminus=ci025)
abline(h=0.5,lty=2)
#plot(vote.diffs,prob.foc.clust.wins,ylab='probability',xlab='vote difference',pch=19,xlim=c(-10,10),ylim=c(0,1),axes=FALSE,col='#000000')
axis(1,at=c(-21,-16,-12,-8,-4,0,4,8,12,16,21)) #change according to group
axis(2)

#fit a sigmoid
df<- data.frame(ninds=vote.diffs,prob=prob.foc.clust.wins)
params <- list(Asym=1, xmid =0.5, scal=0.5)
model<-nls(prob ~ SSlogis(ninds,Asym,xmid,scal) ,data=df, start = params, trace = F, control = list(maxiter = 500))
newx <- seq(max(-30,-max.inits+1),min(30,max.inits-1),0.05)
newy <- predict(model,newdata=data.frame(ninds=newx))
lines(newx,newy/max(newy),col="red",lwd=2)


#save(list=c('clust.level.data'),file="~/cluster_level_data_gmm_thres35.RData")


