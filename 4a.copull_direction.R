#CO-PULLING DIRECTION

#This code computes the probability of pulling success as a function of the level of agreement between initiators and the number of
#initiators, fits a GEE to these data, and outputs a surface plot.

#Note: this uses the data frame built in copull_direction2

#Figure:
# Probability of success vs. number of initiators for different levels of agreement (+GEE fit)

library(Hmisc)
library(geepack)
library(colorRamps)
library(fields)


#FUNCS
circ.var <- function(dv){
	dx <- dv[,1]
	dy <- dv[,2]
	cv <- 1 - sqrt((sum(dx)^2) + (sum(dy)^2))/length(dx)
	return(cv)
	
}

#MAIN
load("~/events_with_directions_thresh35.RData")

events <- events[which(events$not.cohesive.event.cluster == 0),]

n.inds <- length(unique(c(events$leader, events$follower)))

#fun <- function(x) { 
#  if (sum(x == "pull") == length(x)) {
#    return(1)
#  } else if (sum(x == "anchor") == length(x)) {
#    return(0)
#  } else {
#    return(NA)
#  }
#}
#success <- sapply(by(events$type,events$event.num,function(x) {fun(x)}),identity)
#events$success <- success[match(events$event.num,names(success))]

events <- events[which(!is.na(events$success_new)),]

agreement <- 1-sapply(by(cbind(events$init.direc.x,events$init.direc.y),events$event.num,circ.var),identity)
events$agreement <- agreement[match(events$event.num,names(agreement))]

###Make a 2D histogram
cv.bins <- seq(0,1,.1)
n.inds.bins <- seq(2,22,1)
p.success <- array(NA,dim=c(length(n.inds.bins),length(cv.bins)-1))
n.events <- array(NA,dim=c(length(n.inds.bins),length(cv.bins)-1))
for(i in 1:length(n.inds.bins)){
	for(j in 1:(length(cv.bins)-1)){
		events.sub <- events[which(events$n.inits==n.inds.bins[i] & events$agreement >= cv.bins[j] & events$agreement < cv.bins[j+1]),]
		succ <- nrow(events.sub[which(events.sub$success==1),])
		fail <- nrow(events.sub[which(events.sub$success==0),])
		p.success[i,j]<-succ / (succ + fail)
		n.events[i,j]<-succ+fail
		if((succ + fail)<3) {p.success[i,j]<-NA}
	}
}


p.events <- n.events / rowSums(n.events)

#get the mean level of agreement of each number of initiators
mean.agreement <- array(NA,dim=c(length(n.inds.bins)))
agreement.upper <- array(NA,dim=c(length(n.inds.bins)))
agreement.lower <- array(NA,dim=c(length(n.inds.bins)))
for(i in 1:length(n.inds.bins)){
	agreements<-events[which(events$n.inits==n.inds.bins[i]),]$agreement
	mean.agreement[i]<-mean(agreements)
	agreement.upper[i]<-quantile(agreements,0.975)
	agreement.lower[i]<-quantile(agreements,0.025)
	
}

########Probability of success vs number of initiators for different levels of agreement

mygee <- geeglm(success ~ agreement*n.inits, data=events, family="binomial", id=n.inits, corstr="independence")

plotcol="black"
zmax <- 1
zmin <- 0

n.inits2 <- seq(2:n.inds)
cv.bins2 <- seq(0,1,0.01)
z <- matrix(NA,nrow=length(n.inits2),ncol=length(cv.bins2))
for (i in 1:length(cv.bins2)) {
	z[,i] <- predict(mygee,data.frame(agreement=rep(cv.bins2[i],length(n.inits2)),n.inits=n.inits2),type='resp')
}

minz <- min(z)
nrz <- nrow(z)
ncz <- ncol(z)
# Create colors
new.palette <- colorRampPalette(c("blue", "white","red"),space="rgb")
jet.colors  <- new.palette(1000)

nbcol <- 100

# CREATE IMAGE PLOT
par(mar=c(5,5,2,6))

image(z,col=jet.colors,axes=FALSE, xlab="Number of initiators", ylab="Agreement",zlim=c(zmin,zmax))
x.labs <- seq(1,n.inds,3)
y.labs <- seq(0,1,0.2)
axis(1,at=seq(0,1,1/(length(x.labs)-1)),labels=x.labs, cex.axis=0.8)
axis(2,at=seq(0,1,1/(length(y.labs)-1)),labels=y.labs, cex.axis=0.8)
box()
contour(z,plot=FALSE,add=TRUE)
image.plot(legend.only=TRUE,col=jet.colors,zlim=c(zmin,zmax),legend.lab="Probability of success")

library(broom)
tidy_m <- tidy(mygee)
#write.csv(tidy_m, "~/Group1_mygee_table1_thres35.csv")

