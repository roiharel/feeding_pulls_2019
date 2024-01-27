library(reshape); library(lubridate); library(dplyr)

load("~/events_with_directions_thresh35.RData")
load("~/times.Rdata")

times2 <- times
minute(times2) <- 0
second(times2) <- 0

load("~/xs.Rdata")
#xs <-xs[-1,]#that's because indiv. W1261 from Group2 had no data and we have to exclude it from the analysis 

load("~/nas_prop.Rdata")
na <- which(nas_prop$na_prop > 0.5)
#I have to exclude cases when less than 50% of tags where collecting data
xs[,na] <- NA

hours_tracked <- c()

for (i in 1: nrow(xs)){
  a <- which(!is.na(xs[i,]))
  timesa <- times2[a]
  hours_tracked[i] <- length(unique(timesa))
}

load("~/dominance_ranks_wts.Rdata")

events_ <- events[which(events$not.cohesive.event.cluster == 0),]

succY <- cast(events_[which(events_$type == "pull"),], leader ~ type)
succN <- cast(events_[which(events_$type == "anchor"),], leader ~ type)

succY <- succY[order(succY$leader),]
succN <- succN[order(succN$leader),]

succYN <- cbind(succY, succN)

succYN$leader.1 <- rownames(xs)

rank <- data.frame(ids = ids_, rank = c(1:length(ids_)))

succYN$rank <- rank$rank[match(succYN$leader.1, rank$ids)]

tableY <- as.table(succYN$pull)
tableN <- as.table(succYN$anchor)
names(tableY) <- 1:nrow(succYN)
names(tableN) <- 1:nrow(succYN)

a <- rbind(tableY, tableN, hours_tracked)
a_ <- a[c(1,2),]
a_[1,] <- a[1,]/a[3,]
a_[2,] <- a[2,]/a[3,]

b<-barplot(a_/(nrow(xs)^2), xlab="Dominance rank", ylab="Rates of initiations per hour / sqrt(group size)", col=c("black", "darkgrey"),
           border=NA, ylim=c(0,max(a_/(nrow(xs)^2))+0.01), cex.names=0.8)


legend(x=nrow(succYN)-2,y=max(a_/(nrow(xs)^2)),
       legend = c("Pulls", "Anchors"),
       fill = c("black", "darkgrey"), cex=.8, bty="n")


############################## test
library(openxlsx)

succYN <- succYN[order(succYN$rank),]
succYN$pull_rate <- a_[1,]

succYN$new_rank <- c(1:nrow(succYN))

succYN <- succYN[rev(order(succYN$pull_rate)),]
succYN$pull_rate_rank <- c(1:nrow(succYN))

sex <- read.xlsx("~/VGF_habituated_groups.xlsx", sheet=1) #Group 1 is sheet 1 and Group 2 sheet 4
cap <- read.xlsx("~/VGF_capture.xlsx", sheet=1)
sex$id <- NA
sex$id <- cap$Original.Ring.number[match(sex$Ind, cap$Colour.Bands)]
sex$id[is.na(sex$id)] <- cap$Original.Ring.number[match(sex$Ind[is.na(sex$id)], cap$Wing.Tag)]

succYN$sex <- sex$Sex[match(succYN$leader.1, sex$id)]

#all
all <- mean(abs(succYN$new_rank- succYN$pull_rate_rank))

dom_r <- list()
rate_r <- list()
mean_r <- c()

for (i in 1:1000){
  dom_r[[i]] <- sample(1:length(succYN$leader.1), length(succYN$leader.1))
  rate_r[[i]] <- sample(1:length(succYN$leader.1), length(succYN$leader.1))
  mean_r <-c(mean_r, mean(abs(dom_r[[i]]- rate_r[[i]])))
}

mean(mean_r)

#only males
males <- mean(abs(succYN$new_rank[succYN$sex %in% c("M")]- succYN$pull_rate_rank[succYN$sex %in% c("M")])) 

succYN_m <- succYN[succYN$sex %in% c("M"),]
dom_r_m <- list()
rate_r_m <- list()
mean_r_m <- c()


for (i in 1:1000){
  dom_r_m[[i]] <- sample(1:length(succYN_m$leader.1), length(succYN_m$leader.1))
  rate_r_m[[i]] <- sample(1:length(succYN_m$leader.1), length(succYN_m$leader.1))
  mean_r_m <-c(mean_r_m, mean(abs(dom_r_m[[i]]- rate_r_m[[i]])))
}


mean(mean_r_m)

par(mfrow=c(1,3))

#only females
females <- mean(abs(succYN$new_rank[succYN$sex %in% c("F")]- succYN$pull_rate_rank[succYN$sex %in% c("F")])) 

succYN_f <- succYN[succYN$sex %in% c("F"),]
dom_r_f <- list()
rate_r_f <- list()
mean_r_f <- c()

for (i in 1:1000){
  dom_r_f[[i]] <- sample(1:length(succYN_f$leader.1), length(succYN_f$leader.1))
  rate_r_f[[i]] <- sample(1:length(succYN_f$leader.1), length(succYN_f$leader.1))
  mean_r_f <-c(mean_r_f, mean(abs(dom_r_f[[i]]- rate_r_f[[i]])))
}

mean(mean_r_f)

par(mfrow=c(1,1))
plot(c(1,2,3), c(mean(mean_r), mean(mean_r_m), mean(mean_r_f)), xlim=c(0.5,3.5), ylim=c(0,6), pch=16,
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



#more tests
nMales <- length(which(succYN$sex == "M"))

succYN$inf0 <- 0
succYN$inf0[which(succYN$pull_rate_rank < length(which(succYN$sex =="M")))] <- 1
succYN$sexbin <- 0
succYN$sexbin[succYN$sex == "M"] <-1

mean(abs(succYN$inf0- succYN$sexbin))

inf0_r <- list()
mean_t_r <- c()
for (i in 1:1000){
  inf0_r[[i]] <- sample(c(rep(1, nMales),rep(0,length(which(succYN$sex == "F")))),length(succYN$leader.1))
  mean_t_r <- c(mean_t_r, mean(abs(inf0_r[[i]]- succYN$sexbin)))
}

mean(mean_t_r)

plot(1, mean(mean_t_r), ylim=c(0,1), pch=16,
     ylab="", xaxt="n", xlab="", xlim=c(0.8, 1.2))
segments(1,quantile(mean_t_r,0.05),1,quantile(mean_t_r,1))
points(1, mean(abs(succYN$inf0- succYN$sexbin)), col="red", cex=1, pch=8)

p4 <- sum(mean(abs(succYN$inf0- succYN$sexbin)) >= mean_t_r) / length(mean_t_r)

text(1.1,mean(abs(succYN$inf0- succYN$sexbin)),paste("p=",p4), cex=0.7)
text(1,0.9, "(iv)")







