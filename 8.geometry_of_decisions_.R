#GEOMETRY OF DECISIONS
#This code makes plots of direction taken by followers vs angle of disagreement between initiators.
#It also carries out statistical tests to determine whether averaging is occurring at small angles.


library(fields)
library(PerformanceAnalytics)
library(diptest)
library(grDevices)

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
for(i in sort(unique(events$event.num))){
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

#make a histogram of the resultant vector as a function of the difference in angle
theta.bins.y <- seq(-pi,2*pi,pi/30)
theta.bins.x <- seq(0,pi,pi/30)
theta.hist <- matrix(0,nrow=length(theta.bins.x),ncol=length(theta.bins.y)-1)
medians.all <- medians.up <- medians.down <- array(NA,length(theta.bins.x))
for(i in 1:(length(theta.bins.x)-1)){
  th.min <- theta.bins.x[i]
  th.max <- theta.bins.x[i+1]
  curr <- events[which(events$d.theta >= th.min & events$d.theta < th.max),]
  angs <- atan2(curr$result.direc.y.rot,curr$result.direc.x.rot)
  angs[which(angs<= -pi/2)] <- angs[which(angs<= -pi/2)]+2*pi
  histo <- hist(angs,breaks=theta.bins.y,plot=F)
  histo <- histo$counts / sum(histo$counts)
  theta.hist[i,]<-histo
  medians.all[i]<-median(angs,na.rm=T)
  medians.up[i]<-median(angs[which(angs>theta.bins.x[i]/2)])
  medians.down[i]<-median(angs[which(angs<=theta.bins.x[i]/2)])
  
}

new.palette <- colorRampPalette(c("blue3", "white","red"),space="rgb")
new.palette.r <- new.palette(1000)
#make plot of difference in direction of pull vectors vs. probability of choosing each angle
#quartz()
image.plot(theta.bins.x*(180/pi),theta.bins.y*(180/pi),theta.hist,xlab='Difference in angle between pullers',ylab='Direction chosen',ylim=c(-90,180), col=new.palette.r)

abline(h=0,col='black',lty=2,lwd=2)
lines(seq(0,180,1),seq(0,180,1),col='black',lty=2,lwd=2)
lines(theta.bins.x*180/pi,theta.bins.x/2*180/pi,lwd=2,col='black',lty=2)

lines(theta.bins.x[which(theta.bins.x<=130/180*pi)]*180/pi,medians.all[which(theta.bins.x<=130/180*pi)]*180/pi,lwd=4,col='black')
lines(theta.bins.x[which(theta.bins.x<=130/180*pi)]*180/pi,medians.all[which(theta.bins.x<=130/180*pi)]*180/pi,lwd=3,col='white')
#lines(theta.bins.x[which(theta.bins.x<=180/180*pi)]*180/pi,medians.all[which(theta.bins.x<=180/180*pi)]*180/pi,lwd=3,col='white')
lines(theta.bins.x[which(theta.bins.x>=78/180*pi)]*180/pi,medians.up[which(theta.bins.x>=78/180*pi)]*180/pi,lwd=4,col='black')
lines(theta.bins.x[which(theta.bins.x>=78/180*pi)]*180/pi,medians.up[which(theta.bins.x>=78/180*pi)]*180/pi,lwd=3,col='white')
lines(theta.bins.x[which(theta.bins.x>=78/180*pi)]*180/pi,medians.down[which(theta.bins.x>=78/180*pi)]*180/pi,lwd=4,col='black')
lines(theta.bins.x[which(theta.bins.x>=78/180*pi)]*180/pi,medians.down[which(theta.bins.x>=78/180*pi)]*180/pi,lwd=3,col='white')
abline(v=78, col="darkgray")
abline(v=130, col="darkgray")

#test whether they are averaging or choosing based on bimodality of distriubtion at each angle and get standard deviations of upper and lower modes
#library(diptest)
#ps <- array(NA,dim=length(theta.bins.x))
stds.up <- array(NA,dim=length(theta.bins.x))
stds.down <- array(NA,dim=length(theta.bins.x))
for(i in 1:(length(theta.bins.x)-1)){
  th.min <- theta.bins.x[i]
  th.max <- theta.bins.x[i+1]
  curr <- events[which(events$d.theta >= th.min & events$d.theta < th.max),]
  angs <- atan2(curr$result.direc.y.rot,curr$result.direc.x.rot)
  angs[which(angs<= -pi/2)] <- angs[which(angs<= -pi/2)]+2*pi
  #dip.outs<-dip(angs,full.result=FALSE)
  #p <- mean(dip.outs <= replicate(1000,runif(length(angs))))
  #ps[i]<-p
  stds.up[i] <- sd(angs[which(angs > theta.bins.x[i]/2)],na.rm=T)
  stds.down[i] <- sd(angs[which(angs <= theta.bins.x[i]/2)],na.rm=T)
  
}

mean.std <- mean(c(stds.up[which(theta.bins.x*180/pi>=120)],stds.down[which(theta.bins.x*180/pi>=90)]),na.rm=T)

#assume superimposed gaussian distributions with means = angles of opinion and stds = mean.std for each
par(mfrow=c(3,5),mar=c(2,2,1,1))
theta.bins.x <- seq(0,pi,pi/15)
theta.bins.y2 <- seq(-pi,2*pi,pi/30)
kurt.ps <- array(NA,length(theta.bins.x))
bimod.ps <- array(NA,length(theta.bins.x))
dip.ps <- diptest.ps <- array(NA,length(theta.bins.x))
for(i in 1:(length(theta.bins.x)-1)){
  th.min <- theta.bins.x[i]
  th.max <- theta.bins.x[i+1]
  curr <- events[which(events$d.theta >= th.min & events$d.theta < th.max),]
  angs <- atan2(curr$result.direc.y.rot,curr$result.direc.x.rot)
  angs[which(angs<= -pi/2)] <- angs[which(angs<= -pi/2)]+2*pi
  histo.dat <- hist(angs,breaks=theta.bins.y2,plot=F)
  histo.dat$counts <- histo.dat$counts / sum(histo.dat$counts)
  plot(histo.dat$mids*180/pi,histo.dat$counts,type='l',col='red',main=paste(th.min*180/pi,'-',th.max*180/pi,sep=''))
  
  null <- c(rnorm(1000,mean=0,sd=mean.std),rnorm(1000,mean=(theta.bins.x[i]+theta.bins.x[i+1])/2,sd=mean.std))
  if(length(which(null <= -pi/2))>0){
    null[which(null <= -pi/2)] <- null[which(null <= -pi/2)]+2*pi
  }
  histo <- hist(null,plot=F,breaks=theta.bins.y2)
  histo$counts <- histo$counts / sum(histo$counts)
  lines(histo.dat$mids*180/pi,histo$counts)
  
  #get dip statistic
  dip.data <- dip(angs)
  diptest.data <- dip.test(angs)
  
  #run kurtosis test
  kurt.dat <- kurtosis(angs,method='excess')
  skew.dat <- skewness(angs)
  N <- length(angs)
  bimod.dat <- (skew.dat^2+1)/(kurt.dat+((3*(N-1)^2)/((N-2)*(N-3))))
  kurts.null <- array(NA,1000)
  skews.null <- array(NA,1000)
  dips.null <- array(NA,1000)
  for(j in 1:1000){
    null <- c(rnorm(N/2,mean=0,sd=mean.std),rnorm(N/2,mean=(theta.bins.x[i]+theta.bins.x[i+1])/2,sd=mean.std))
    kurts.null[j] <- kurtosis(null)
    skews.null[j] <- skewness(null)
    dips.null[j] <- dip(null)
  }
  kurt.ps[i]<-sum(kurts.null > kurt.dat)/1000
  bimods.null <- (skews.null^2+1)/(kurts.null+((3*(N-1)^2)/((N-2)*(N-3))))
  bimod.ps[i] <- sum(bimods.null > bimod.dat)/1000
  dip.ps[i] <- sum(dips.null >= dip.data)/1000
  diptest.ps[i] <- diptest.data$p.value 
}
#quartz(width=8,height=7.5)

plot(theta.bins.x*180/pi,1-bimod.ps,pch=19,xlab='difference in angle',ylab='p-value of bimodality coefficient compared with null distributions',col=c("black","black")[1+((1-bimod.ps)<0.05)],main='2 pullers',ylim=c(0,1))
points(theta.bins.x*180/pi,diptest.ps,pch=19,col='red')
abline(h=.05,lty=2)



