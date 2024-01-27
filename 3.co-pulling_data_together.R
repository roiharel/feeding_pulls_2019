## PUTTING DATA TOGETHER
#Cluster dyadic interactions into "events" involving concurrent initiations and the same follower.
setwd("path")

load("Group1_events_thresholds_30_35_01_01.Rdata")

mydata <- events <- events_
NB <- max(max(events$leader), max(events$follower)) 
first <- min(min(events$leader), min(events$follower)) 


days <- rep(0,NB)
for (i in first:NB) {
	days[i] <- length(unique(mydata$day[mydata$leader==i]))
}

# PUT DATA TOGETHER
mydata <- mydata[order(mydata$day,mydata$t2),]
mydata <- mydata[which(mydata$strength >= 0.1 & mydata$disp >= 0.1),]

mydata_anchor <- mydata[which(mydata$type == "anchor"),]
mydata_pull <- mydata[which(mydata$type == "pull"),]
mydata <- mydata[which(mydata$type %in% c("pull","anchor")),]

mydata_anchor$event.num <- NA
cur_event <- 0
for (i in c(first:NB)) {
	cur_event <- cur_event + 1
	indices <- which(mydata_anchor$follower == i)
	mydata_anchor$event.num[indices][1] <- cur_event
	for (j in 2:nrow(mydata_anchor[indices,])) {
		if (mydata_anchor$t1[indices][j] <= mydata_anchor$t2[indices][j-1]) {
			mydata_anchor$event.num[indices][j] <- cur_event
		} else {
			cur_event <- cur_event + 1
			mydata_anchor$event.num[indices][j] <- cur_event
		}
	}
}


mydata_pull$event.num <- NA
cur_event <- 0
for (i in c(first:NB)) {
	cur_event <- cur_event + 1
	indices <- which(mydata_pull$follower == i)
	mydata_pull$event.num[indices][1] <- cur_event
	for (j in 2:nrow(mydata_pull[indices,])) {
		if (mydata_pull$t1[indices][j] <= mydata_pull$t2[indices][j-1]) {
			mydata_pull$event.num[indices][j] <- cur_event
		} else {
			cur_event <- cur_event + 1
			mydata_pull$event.num[indices][j] <- cur_event
		}
	}
}


# ALL

tmp <- mydata$leader[mydata$type=="anchor"]
mydata$leader[mydata$type=="anchor"] <- mydata$follower[mydata$type=="anchor"]
mydata$follower[mydata$type=="anchor"] <- tmp
mydata$event.num <- NA
cur_event <- 0
for (i in c(first:NB)) {
	cur_event <- cur_event + 1
	indices <- which(mydata$follower == i)
	mydata$event.num[indices][1] <- cur_event
	for (j in 2:nrow(mydata[indices,])) {
		if (mydata$t1[indices][j] <= mydata$t2[indices][j-1]) {
			mydata$event.num[indices][j] <- cur_event
		} else {
			cur_event <- cur_event + 1
			mydata$event.num[indices][j] <- cur_event
		}
	}
}




# PUT DATA TOGETHER: ORDERINGS
mydata_pull$Order <- NA
mydata_pull$Order.norm <- NA
mydata_pull$time <- NA
mydata_pull <- mydata_pull[order(mydata_pull$day,mydata_pull$t1),]
indices <- which(table(mydata_pull$event.num[which(mydata_pull$disp>=0.1 & mydata_pull$strength>=0.1)]) > 1)

for (i in indices) {
	mydata_pull$Order[which(mydata_pull$event.num == i & mydata_pull$disp>=0.1 & mydata_pull$strength>=0.1)] <- c(1:sum(mydata_pull$event.num[which(mydata_pull$disp>=0.1 & mydata_pull$strength>=0.1)] == i))
	mydata_pull$Order.norm[which(mydata_pull$event.num == i & mydata_pull$disp>=0.1 & mydata_pull$strength>=0.1)] <- (mydata_pull$Order[which(mydata_pull$event.num == i & mydata_pull$disp>=0.1 & mydata_pull$strength>=0.1)]-1)/sum(which(mydata_pull$event.num == i & mydata_pull$disp>=0.1 & mydata_pull$strength>=0.1))
	mydata_pull$time[which(mydata_pull$event.num == i & mydata_pull$disp>=0.1 & mydata_pull$strength>=0.1)] <- mydata_pull$t1[which(mydata_pull$event.num == i & mydata_pull$disp>=0.1 & mydata_pull$strength>=0.1)]-min(mydata_pull$t1[which(mydata_pull$event.num == i & mydata_pull$disp>=0.1 & mydata_pull$strength>=0.1)])
}



mydata_pull2 <- mydata_pull[which(!is.na(mydata_pull$Order)),]
order_of_pulling <- rep(NA,NB)
for (i in first:NB) {
	order_of_pulling[i] <- mean(sapply(by(mydata_pull2$Order.norm[mydata_pull2$leader==i], mydata_pull2$event.num[mydata_pull2$leader==i],min),identity))
}

mydata_anchor$Order <- NA
mydata_anchor$Order.norm <- NA
mydata_anchor$time <- NA
mydata_anchor <- mydata_anchor[order(mydata_anchor$day,mydata_anchor$t1),]
indices <- which(table(mydata_anchor$event.num[which(mydata_anchor$disp>=0.1 & mydata_anchor$strength>=0.1)]) > 1)

for (i in indices) {
	mydata_anchor$Order[which(mydata_anchor$event.num == i & mydata_anchor$disp>=0.1 & mydata_anchor$strength>=0.1)] <- c(1:sum(mydata_anchor$event.num[which(mydata_anchor$disp>=0.1 & mydata_anchor$strength>=0.1)] == i))
	mydata_anchor$Order.norm[which(mydata_anchor$event.num == i & mydata_anchor$disp>=0.1 & mydata_anchor$strength>=0.1)] <- (mydata_anchor$Order[which(mydata_anchor$event.num == i & mydata_anchor$disp>=0.1 & mydata_anchor$strength>=0.1)]-1)/sum(which(mydata_anchor$event.num == i & mydata_anchor$disp>=0.1 & mydata_anchor$strength>=0.1))
	mydata_anchor$time[which(mydata_anchor$event.num == i & mydata_anchor$disp>=0.1 & mydata_anchor$strength>=0.1)] <- mydata_anchor$t1[which(mydata_anchor$event.num == i & mydata_anchor$disp>=0.1 & mydata_anchor$strength>=0.1)]-min(mydata_anchor$t1[which(mydata_anchor$event.num == i & mydata_anchor$disp>=0.1 & mydata_anchor$strength>=0.1)])
}

 mydata_anchor2 <- mydata_anchor[which(!is.na(mydata_anchor$Order)),]
order_of_anchor <- rep(NA,NB)
for (i in first:NB) {
	order_of_anchor[i] <- mean(sapply(by(mydata_anchor2$Order.norm[mydata_anchor2$leader==i], mydata_anchor2$event.num[mydata_anchor2$leader==i],max),identity))
}


save(mydata_pull,file="~/co-pulled_thresh35.Rdata")
save(mydata_anchor,file="~/co-anchored_thresh35.Rdata")
save(mydata,file="~/co-pull_co-anchor_thresh35.Rdata")
