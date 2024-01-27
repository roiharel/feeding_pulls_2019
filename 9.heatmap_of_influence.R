library(lattice)
library(openxlsx)

load("~/events_with_directions_thresh35.RData")

events_ <- events[which(events$not.cohesive.event.cluster == 0),]

#events_$leader <-events_$leader-1 #that's because indiv. W1261 from Group 2 had no data and we had to exclude it from the analysis otherwise it gives a white column in the heatmap
#events_$follower <-events_$follower-1 


foltab <- matrix(NA, max(events_$leader), max(events_$leader))
load("~/xs.Rdata")

#xs <-xs[-1,]

rownames(foltab) <- rownames(xs)
colnames(foltab) <- rownames(xs)

l<- matrix(NA, max(events_$leader), max(events_$leader))

for(i in 1:max(events_$leader)){
  for(j in 1:max(events_$leader)){
    if(i != j){
  l[i,j] <- length(which(events_$leader == i & events_$follower == j))
    }
  }
}

for(i in 1:max(events_$leader)){
  for(j in 1:max(events_$leader)){
    if(i != j){
      foltab[i,j] <- (l[i,j] - l[j,i])/(l[i,j] + l[j,i])
    }
  }
}

rSc <- rowSums(foltab, na.rm=T)

rSc <- rSc[order(rSc)]
idscnts <- names(rSc)

foltab_ <- foltab[idscnts,idscnts]

new.palette <- colorRampPalette(c("blue", "white","red"),space="rgb")
levelplot(foltab_[ncol(foltab_):1,ncol(foltab_):1],cuts=49,col.regions=new.palette(50), xlab="influencer", ylab="influenced", colorkey=list(col=new.palette(50)))

#plot against dominance rank
load("~/dominance_ranks_Group1.Rdata")# load ids of individuals in descending dominance rank
ids_ <- ids_[which(ids_ %in% idscnts)]

foltab_ <- foltab[rev(ids_),rev(ids_)]

new.palette <- colorRampPalette(c("blue", "white","red"),space="rgb")
levelplot(foltab_[ncol(foltab_):1,ncol(foltab_):1],cuts=49,col.regions=new.palette(50), xlab="Influencer", ylab="Influenced", colorkey=list(col=new.palette(50)), scales=list(x=list(rot=90)))


###############statistical testing######################
inf_dom <- data.frame(ID = ids_, dom= 1:length(ids_), infl= match(rev(idscnts), ids_))

sex <- read.xlsx("~/VGF_habituated_groups.xlsx", sheet=1)
cap <- read.xlsx("~/VGF_capture.xlsx", sheet=1)
sex$id <- NA
sex$id <- cap$Original.Ring.number[match(sex$Ind, cap$Colour.Bands)]
sex$id[is.na(sex$id)] <- cap$Original.Ring.number[match(sex$Ind[is.na(sex$id)], cap$Wing.Tag)]

inf_dom$sex <- sex$Sex[match(inf_dom$ID, sex$id)]

rm(list=setdiff(ls(), "inf_dom"))

females <- mean(abs(inf_dom$dom[inf_dom$sex %in% c("F")]- inf_dom$infl[inf_dom$sex %in% c("F")])) 

#all
all <- mean(abs(inf_dom$dom- inf_dom$infl))

dom_r <- list()
inf_r <- list()
mean_r <- c()
  
for (i in 1:1000){
dom_r[[i]] <- sample(1:length(inf_dom$ID), length(inf_dom$ID))
inf_r[[i]] <- sample(1:length(inf_dom$ID), length(inf_dom$ID))
mean_r <-c(mean_r, mean(abs(dom_r[[i]]- inf_r[[i]])))
}

mean(mean_r)



#only males
males <- mean(abs(inf_dom$dom[inf_dom$sex %in% c("M")]- inf_dom$infl[inf_dom$sex %in% c("M")])) 

inf_dom_m <- inf_dom[inf_dom$sex %in% c("M"),]
dom_r_m <- list()
inf_r_m <- list()
mean_r_m <- c()

for (i in 1:1000){
  dom_r_m[[i]] <- sample(1:length(inf_dom_m$ID), length(inf_dom_m$ID))
  inf_r_m[[i]] <- sample(1:length(inf_dom_m$ID), length(inf_dom_m$ID))
  mean_r_m <-c(mean_r_m, mean(abs(dom_r_m[[i]]- inf_r_m[[i]])))
}


mean(mean_r_m)

par(mfrow=c(1,3))

#boxplot(mean_r_m, ylim=c(2,12), main="All")
#points(1, all, col="red", cex=3, pch=8)

#only females
females <- mean(abs(inf_dom$dom[inf_dom$sex %in% c("F")]- inf_dom$infl[inf_dom$sex %in% c("F")])) 

inf_dom_f <- inf_dom[inf_dom$sex %in% c("F"),]
dom_r_f <- list()
inf_r_f <- list()
mean_r_f <- c()

for (i in 1:1000){
  dom_r_f[[i]] <- sample(1:length(inf_dom_f$ID), length(inf_dom_f$ID))
  inf_r_f[[i]] <- sample(1:length(inf_dom_f$ID), length(inf_dom_f$ID))
  mean_r_f <-c(mean_r_f, mean(abs(dom_r_f[[i]]- inf_r_f[[i]])))
}

mean(mean_r_f)


#par(mfrow=c(1,3))
#boxplot(mean_r, ylim=c(1,12), main="All")
#points(1, all, col="red", cex=3, pch=8)

#boxplot(mean_r_m, ylim=c(1,12), main="Males")
#points(1, males, col="red", cex=3, pch=8)

#boxplot(mean_r_f, ylim=c(1,12), main="Females")
#points(1, females, col="red", cex=3, pch=8)



#more tests
nMales <- length(which(inf_dom$sex == "M"))

inf_dom$inf0 <- 0
inf_dom$inf0[which(inf_dom$infl < length(which(inf_dom$sex =="M")))] <- 1
inf_dom$sexbin <- 0
inf_dom$sexbin[inf_dom$sex == "M"] <-1

mean(abs(inf_dom$inf0- inf_dom$sexbin))

inf0_r <- list()
mean_t_r <- c()
for (i in 1:1000){
  inf0_r[[i]] <- sample(c(rep(1, nMales),rep(0,length(which(inf_dom$sex == "F")))),length(inf_dom$ID))
  mean_t_r <- c(mean_t_r, mean(abs(inf0_r[[i]]- inf_dom$sexbin)))
}

mean(mean_t_r)


#par(mfrow=c(1,3))
#boxplot(mean_t_r, ylim=c(0,1), main="Test 2")
#points(1, mean(abs(inf_dom$inf0- inf_dom$sexbin)), col="red", cex=3, pch=8)

#final plot
par(mfrow=c(1,2))
plot(c(1,2,3), c(mean(mean_r), mean(mean_r_m), mean(mean_r_f)), xlim=c(0.5,3.5), ylim=c(0,13), pch=16,
     ylab="", xaxt="n", xlab="")
segments(1,quantile(mean_r,0.05),1,quantile(mean_r,1))
segments(2,quantile(mean_r_m,0.05),2,quantile(mean_r_m,1))
segments(3,quantile(mean_r_f,0.05),3,quantile(mean_r_f,1))
points(1, all, col="red", cex=1, pch=8)
points(2, males, col="red", cex=1, pch=8)
points(3, females, col="red", cex=1, pch=8)
axis(1,at=1:3, labels=c("All", "Males", "Females"))

p1 <- sum(all >= mean_r) / length(mean_r)
p2 <- sum(males >= mean_r_m) / length(mean_r_m)
p3 <- sum(females >= mean_r_f) / length(mean_r_f)

text(1.3,all,paste("p=", round(p1, digits=3)), cex=0.5)
text(2.3,males, paste("p=", round(p2,digits = 3)), cex=0.5)
text(3.3,females,paste("p=", round(p3,digits = 3)), cex=0.5)

text(c(1:3), rep(7,3), c("(i)", "(ii)", "(iii)"))


plot(1, mean(mean_t_r), ylim=c(0,1), pch=16,
     ylab="", xaxt="n", xlab="", xlim=c(0.8, 1.2))
segments(1,quantile(mean_t_r,0.05),1,quantile(mean_t_r,1))
points(1, mean(abs(inf_dom$inf0- inf_dom$sexbin)), col="red", cex=1, pch=8)

p4 <- sum(mean(abs(inf_dom$inf0- inf_dom$sexbin)) >= mean_t_r) / length(mean_t_r)

text(1.1,mean(abs(inf_dom$inf0- inf_dom$sexbin)),paste("p=",p4), cex=0.7)
text(1,0.9, "(iv)")



