#COPULLING DIRECTION
#Reads in pull / anchor data (grouped into events) and computes the direction of pulling for each initiator in each event. 

source("~/pull_anchor_funcs.R")

#FUNCTIONS
frac.pulls<-function(types){
	n.pull <- sum(types=='pull')
	n.anch <- sum(types=='anchor')
	return(n.pull / (n.pull + n.anch))
}

is.mixed<-function(types){
	n.pull <- sum(types=='pull')
	n.anch <- sum(types=='anchor')
	if(n.pull==0){return(F)}
	if(n.anch==0){return(F)}
	return(T)
}

min.pull.t2<-function(types.and.t2s){
	if(sum(types.and.t2s[,1])>0){return(min(types.and.t2s[which(types.and.t2s[,1]==1),2],na.rm=T))}
	else{return(NA)}
	
	
}

max.pull.t3<-function(types.and.t3s){
	if(sum(types.and.t3s[,1])>0){return(max(types.and.t3s[which(types.and.t3s[,1]==1),2],na.rm=T))}
	else{return(NA)}
	
	
}

#Load data

load("~/co-pull_co-anchor_thresh35.RData")

events <- mydata[which(mydata$not.cohesive.event.cluster == 0),] # only keep events in which the leaders and followers where within 30 meters from each other

events <- events[which(events$na_prop <= 0.5),] # only keeps events in which more than 50% of the group's tags were working

load("~/xs.Rdata")
load("~/ys.Rdata")


#Make a histogram of fraction of pulls in each event
frac.pulls.by.event<-sapply(by(data=events$type,INDICES=events$event.num,FUN=frac.pulls),identity)
#hist(frac.pulls.by.event,breaks=40,col='gray',ylab='frequency',xlab='fraction of pulls within each event',main='')

#number of mixed events (events where there are some pulls and some anchors)
n.mixed<-length(frac.pulls.by.event) - sum(frac.pulls.by.event==0) - sum(frac.pulls.by.event==1)
frac.mixed <- n.mixed / length(frac.pulls.by.event)

#add a column to events data frame indicating whether or not the event is mixed (containing both pulls and anchors)
mixed.by.event<-aggregate(x=events$type,by=list(events$event.num),FUN=is.mixed)
colnames(mixed.by.event)<-c('event.num','is.mixed')
events$is.mixed<-mixed.by.event$is.mixed[match(x=events$event.num,table=mixed.by.event$event.num)]

#add a column to events data frame indicating whether the event was successful (i.e had at least one pull)
succ.by.event<-aggregate(x=events$type,by=list(events$event.num),FUN=frac.pulls)
colnames(succ.by.event)<-c('event.num','frac.pulls')
succ.by.event$success<-0
succ.by.event$success[which(succ.by.event$frac.pulls>0)]<-1
events$success<-succ.by.event$success[match(x=events$event.num,table=succ.by.event$event.num)]

#add a column to events data frame indicating the number of initiators involved in the event
inits.by.event<-aggregate(x=events$type,by=list(events$event.num),FUN=length)
colnames(inits.by.event)<-c('event.num','n.inits')
events$n.inits<-inits.by.event$n.inits[match(x=events$event.num,table=succ.by.event$event.num)]

#get vector from follower to leader at time t2 of each event, and append to data frame
xi <- xs[cbind(events$leader,events$t2)]
yi <- ys[cbind(events$leader,events$t2)]
xf <- xs[cbind(events$follower,events$t2)]
yf <- ys[cbind(events$follower,events$t2)]
events$pull.direcx.t2 <- (xi - xf) / sqrt((xi-xf)^2 + (yi-yf)^2)
events$pull.direcy.t2 <- (yi - yf) / sqrt((xi-xf)^2 + (yi-yf)^2)
xf3 <- xs[cbind(events$follower,events$t3)]
yf3 <- ys[cbind(events$follower,events$t3)]
events$move.direcx.t3 <- (xf3 - xf) / sqrt((xf3 - xf)^2 + (yf3 - yf)^2)
events$move.direcy.t3 <- (yf3 - yf) / sqrt((xf3 - xf)^2 + (yf3 - yf)^2)

#get min value of t2 for all pulls and anchors within each events (min.t2)
data.by.event<-aggregate(x=events$event.num,by=list(events$event.num),FUN=identity)
colnames(data.by.event)<-c('event.num','min.t2')
data.by.event$min.t2 <- sapply(by(data=events$t2,INDICES=events$event.num,FUN=min),identity)
events$min.t2 <- data.by.event$min.t2[match(x=events$event.num,table=data.by.event$event.num)]

#get max value of t2 for all pulls and anchors within each events (max.t2)
data.by.event<-aggregate(x=events$event.num,by=list(events$event.num),FUN=identity)
colnames(data.by.event)<-c('event.num','max.t2')
data.by.event$max.t2 <- sapply(by(data=events$t2,INDICES=events$event.num,FUN=max),identity)
events$max.t2 <- data.by.event$max.t2[match(x=events$event.num,table=data.by.event$event.num)]

#get min value of t2 for the pulls within each event (min.pull.t2)
events$pull <- 0
events$pull[which(events$type=='pull')]<-1
data.by.event<-aggregate(x=events$event.num,by=list(events$event.num),FUN=identity)
colnames(data.by.event)<-c('event.num','min.pull.t2')
data.by.event$min.pull.t2 <- sapply(by(data=cbind(events$pull,events$t2),INDICES=events$event.num,FUN=min.pull.t2),identity)
events$min.pull.t2 <- data.by.event$min.pull.t2[match(x=events$event.num,table=data.by.event$event.num)]

#get max value of t3 for the pulls within each event (max.pull.t3)
data.by.event$max.pull.t3 <- sapply(by(data=cbind(events$pull,events$t3),INDICES=events$event.num,FUN=max.pull.t3),identity)
events$max.pull.t3 <- data.by.event$max.pull.t3[match(x=events$event.num,table=data.by.event$event.num)]

#get direction of the initiation vector
#this is defined as the vector between the position of the leader at its t2 and the position of the follower at min.pull.t2
#when there is at least one pull in the event. otherwise, it is defined as the vector between the position of the leader at its t2
#and the mean position of follower over the time range min(t2) - max(t2)

	#position of the initiator at its t2
xi.t2<-xs[cbind(events$leader,events$t2)]
yi.t2<-ys[cbind(events$leader,events$t2)]

	#for events that involve at least one pull, use the position of the follower at time min.pull.t2
xf.t2<-xs[cbind(events$follower,events$min.pull.t2)]
yf.t2<-ys[cbind(events$follower,events$min.pull.t2)]

	#for events that involve only anchors, use the mean position of the follower over the time range min.t2 to max.t2
all.anchor.idxs <- which(events$success==0)
for(i in all.anchor.idxs){
	curr.x<-xs[events$follower[i],events$min.t2[i]:events$max.t2[i]]
	curr.y<-ys[events$follower[i],events$min.t2[i]:events$max.t2[i]]
	xf.t2[i]<-mean(curr.x,na.rm=T)
	yf.t2[i]<-mean(curr.y,na.rm=T)
}

dx<-xi.t2 - xf.t2
dy<-yi.t2 - yf.t2
events$init.direc.x<-dx/sqrt(dx^2+dy^2)
events$init.direc.y<-dy/sqrt(dx^2+dy^2)

#get direction of the follower's resultant movement
#this is defined as the vector between the position of the follower at min.pull.t2 and its position at max.pull.t3
#for events that only involve anchors, it is undefined (NA)
xf.max.pull.t3<-xs[cbind(events$follower,events$max.pull.t3)]
yf.max.pull.t3<-ys[cbind(events$follower,events$max.pull.t3)]
dx<-xf.max.pull.t3 - xf.t2
dy<-yf.max.pull.t3 - yf.t2
events$result.direc.x<-dx/sqrt(dx^2+dy^2)
events$result.direc.y<-dy/sqrt(dx^2+dy^2)


# this code adds NAs in the success column if both anchors and pulls are in one event
fun <- function(x) { 
  if (sum(x == "pull") == length(x)) {
    return(1)
  } else if (sum(x == "anchor") == length(x)) {
    return(0)
  } else {
    return(NA)
  }
}
success <- sapply(by(events$type,events$event.num,function(x) {fun(x)}),identity)
events$success_new <- success[match(events$event.num,names(success))]

#save results
save(events, file="~/events_with_directions_thresh35.RData")






