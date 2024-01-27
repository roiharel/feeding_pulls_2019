elo.scores <- function(winners,losers,n.inds=NULL,sigmoid.param=1/100,K=200,init.score=0,n.rands=1000,return.trajectories=FALSE){
	if(is.null(n.inds)){
		n.inds <- max(c(unique(winners),unique(losers)))	
	}
	
	T <- length(winners)
	
	if(return.trajectories){
		all.scores <- array(0,c(n.inds,T+1,n.rands))
	} else{
		all.scores <- array(0,c(n.inds,n.rands))
	}
	
	if (length(K) == 1) {
		K <- rep(K,T)
	}
	
	for(r in 1:n.rands){
		ord <- sample(1:T,T,replace=F)
		winners.perm <- winners[ord]
		losers.perm <- losers[ord]
		scores<-array(NA,c(n.inds,T+1))
		scores[,1]<-init.score
		
		for(i in 1:T){
			
			scores[,i+1] <- scores[,i]
			
			winner <- winners.perm[i]
			loser <- losers.perm[i]
			p<-1/(1+exp(-sigmoid.param*(scores[winner,i]-scores[loser,i]))) #prob that winner wins
			
			if(scores[winner,i] >= scores[loser,i]){
				scores[winner,i+1] <- scores[winner,i] + (1-p)*K[i]
				scores[loser,i+1] <- scores[loser,i] - (1-p)*K[i]
			}
			else{
				scores[winner,i+1] <- scores[winner,i] + p*K[i]
				scores[loser,i+1] <- scores[loser,i] - p*K[i]
			}
		}
		
		if(return.trajectories){
			all.scores[,,r]<-scores
		} else{
			all.scores[,r]<-scores[,T+1]
		}
		
		
	}
	
	invisible(all.scores)	
}




plot_matrix_age_sex_class <- function(mat, orders=NULL, remove=NULL, zlims=NULL, ids_filename="~/Dropbox/baboons_shared/data/IDs.csv",colors='redblue') {

	require(fields)
	
	ids <- read.csv(ids_filename,stringsAsFactors=FALSE)
	colours <- data.frame(SEX=c("M","M","M","F","F"),AGE=c("A","SA","J","A","SA"),COLOUR=c("black","blue","yellow","red","green"),COLOUR_RGB=c("0,0,255","102,204,255","102,102,102","255,0,0","255,204,102"),COLOUR_HEX=c("#0000FF","#66CCFF","#666666","#FF0000","#FFCC66"))
	for (i in 1:nrow(ids)) { ids$colours[i] <- as.character(colours$COLOUR_HEX[which(colours$SEX==ids$Sex[i] & colours$AGE==ids$Age[i])]) }

	if(colors=='redblue'){
		cr <- colorRampPalette(c("blue","white","red"),interpolate="linear")
	}
	if(colors=='greenred'){
		cr <- colorRampPalette(c('red','black','green'))
	}
	
	if (is.null(zlims)) {
		zlims <- c(min(mat,na.rm=T),max(mat,na.rm=T))
	}
	
	if (is.null(orders)) {
		orders <- c(1:ncol(mat))
	}
	
	if (!is.null(remove)) {
		remove <- remove[order(remove,decreasing=TRUE)]
		for (i in 1:length(remove)) {
			orders <- orders[-which(orders == remove[i])]
			orders[orders > remove[i]] <- orders[orders > remove[i]]-1
		}
		mat <- mat[-remove,-remove]
		ids <- ids[-remove,]
	}
	
	n.inds <- ncol(mat)

	ids$Days2 <- ids$Days
	ids$Days2[ids$Days2>=5] <- 14

	quartz(width=8,height=7)
	par(mar=c(5,5,2,7),las=2,mgp=c(3,2,0),xpd=NA)
	image(mat[orders,orders],col=cr(64),axes=FALSE,zlim=zlims)
	axis(1,axTicks(1),at=c(0:(n.inds-1))/(n.inds-1),labels=FALSE)
	axis(2,axTicks(1),at=c(0:(n.inds-1))/(n.inds-1),labels=FALSE)
	text(x=c(0:(n.inds-1))/(n.inds-1), y=-0.12, labels=ids$Collar[orders],col = gray.colors(14,start=0,end=0.5)[15-ids$Days2[orders]], srt=90)
	text(y=c(0:(n.inds-1))/(n.inds-1), x=-0.12, labels=ids$Collar[orders],col = gray.colors(14,start=0,end=0.5)[15-ids$Days2[orders]])
	points(c(0:(n.inds-1))/(n.inds-1),rep(-1.3/(n.inds-1),n.inds),pch=21,cex=2,bg=ids$colours[orders],col = gray.colors(14,start=0,end=0.5)[15-ids$Days2[orders]])
	points(rep(-1.3/(n.inds-1),n.inds),c(0:(n.inds-1))/(n.inds-1),pch=21,cex=2,bg=ids$colours[orders],col = gray.colors(14,start=0,end=0.5)[15-ids$Days2[orders]])
	box()
	image.plot(mat, legend.only=TRUE,col=cr(64),zlim=zlims)
	
}


plot_vector_age_sex_class<-function(x,y,remove=NULL,ids_filename="~/Dropbox/baboons_shared/data/IDs.csv",xlab='',ylab='',xlim=c(min(x,na.rm=T)-diff(range(x,na.rm=T))*0.1,max(x,na.rm=T)+diff(range(x,na.rm=T))*0.1),ylim=c(min(y,na.rm=T)-diff(range(y,na.rm=T))*0.1,max(y,na.rm=T)+diff(range(y,na.rm=T))*0.1), CIs_lower=NULL, CIs_upper=NULL, rank_x=FALSE, rank_y=FALSE, width=8,height=7.5, add=FALSE){
	ids <- read.csv(ids_filename,stringsAsFactors=FALSE)
	colours <- data.frame(SEX=c("M","M","M","F","F"),AGE=c("A","SA","J","A","SA"),COLOUR=c("black","blue","yellow","red","green"),COLOUR_RGB=c("0,0,255","102,204,255","102,102,102","255,0,0","255,204,102"),COLOUR_HEX=c("#0000FFBB","#66CCFFBB","#666666BB","#FF0000BB","#FFCC66BB"))
	for (i in 1:nrow(ids)) { ids$colours[i] <- as.character(colours$COLOUR_HEX[which(colours$SEX==ids$Sex[i] & colours$AGE==ids$Age[i])]) }
	
	if(!is.null(remove)){
		colours <- colours[-remove]
		x <- x[-remove]
		y <- y[-remove]
		ids <- ids[-remove,]
		CIs_lower <- CIs_lower[-remove]
		CIs_upper <- CIs_upper[-remove]
	}
	
	ids$Days2 <- ids$Days
	ids$Days2[ids$Days2>=5] <- 14

	if (rank_x==TRUE) {
		x <- rank(x)
	}
	
	if (rank_y==TRUE) {
		y <- rank(y)
		if (!is.null(CIs_lower) & !is.null(CIs_upper)) {
			CIs_lower <- rank(CIs_lower)
			CIs_upper <- rank(CIs_upper)
		}
	}
	if (!is.null(CIs_lower) & !is.null(CIs_upper)) {
		ylim=c(min(CIs_lower,na.rm=T)-diff(range(y,na.rm=T))*0.1,max(CIs_upper,na.rm=T)+diff(range(y,na.rm=T))*0.1)
	}
	
	if (add == FALSE) {
		quartz(width=width,height=height)
	}
	par(bty='n')
	plot(NULL,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim)
	if (!is.null(CIs_lower) & !is.null(CIs_upper)) {
		arrows(x,CIs_lower,x,CIs_upper,angle=90,code=3,col=ids$colours,lwd=2,length=0.15)
	}
	points(x,y,col=ids$colours,pch=19,cex=6)
	points(x,y,col=gray.colors(14,start=0,end=0.5)[15-ids$Days2],pch=21,cex=6,lwd=1.5)
	text(cbind(x,y),labels=ids$Collar,col=gray.colors(14,start=0,end=0.5)[15-ids$Days2])
}

circ.var <- function(dv){
	dx <- dv[,1]
	dy <- dv[,2]
	cv <- 1 - sqrt((sum(dx)^2) + (sum(dy)^2))/length(dx)
	return(cv)
	
}

#get pull directionality network
get.pull.network<-function(data,N=26,remove=16,n.rands=1){
	E <- matrix(NA,nrow=N,ncol=N)
	winners <- losers <- vector()
	for(i in 1:N){
		for(j in 1:N){
			E[i,j]<-length(which(data$leader==i & data$follower==j & data$type=='pull'))
			winners <- c(winners,rep(i,E[i,j]))
			losers <- c(losers,rep(j,E[i,j]))
		}
	}
	D <- (E - t(E)) / (E + t(E))
	
	elos <- elo.scores(winners,losers,n.rands=n.rands)
	elo.order <- order(rowMeans(elos))
	elo.ranks <- rank(rowMeans(elos))
	plot_matrix_age_sex_class(D,orders=elo.order,remove=remove)
	return(list(E,D,elo.ranks))
	
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

prob.map.2d <- function(x,y,response,xbins,ybins){
	probs <- tots <- array(NA,dim=c(length(xbins)-1,length(ybins)-1))
	for(i in 1:(length(xbins)-1)){
		for(j in 1:(length(ybins)-1)){
			curr <- response[which((x >= xbins[i]) & (x < xbins[i+1]) & (y >= ybins[j]) & (y < ybins[j+1]))]
			probs[i,j] <- mean(curr,na.rm=T)
			tots[i,j] <-  length(curr)
		}
	}
	print(probs)
	out <- list(probs=probs,tots=tots,xbins=xbins,ybins=ybins)
	return(out)
}

i<-1
j<-1


distance <- function(N1,E1,N2,E2) {
	dist <- sqrt((E1 - E2)^2+(N1 - N2)^2)
	return(dist)
}

hist.2d <- function(x,y,xbins,ybins,output_freqs=F,output_plot=T,xlab='',ylab=''){
	nx <- length(xbins)
	ny <- length(ybins)
	histo <- array(0,dim=c(nx-1,ny-1))
	for(i in 1:(nx-1)){
		for(j in 1:(ny-1)){
			histo[i,j] <- length(which((x >= xbins[i]) & (x < xbins[i+1]) & (y >= ybins[j]) & (y < ybins[j+1])))
		}
	}
	if(!output_freqs){
		histo <- histo / sum(histo,na.rm=T)
	}
	
	if(output_plot){
		image.plot(histo,x=xbins[1:(nx-1)],y=ybins[1:(ny-1)],xlab=xlab,ylab=ylab)
	}
	invisible(histo)
	
}

mantel.test <- function(A,B,method='pearson',n.rands=1000){
	N <- dim(A)[1]
	cor.data <- cor(as.vector(A),as.vector(B),method=method,use='pairwise.complete.obs')
	null.corrs <- rep(NA,n.rands)
	for(r in 1:n.rands){
		swap.order <- sample(N)
		A.rand <- A[swap.order,swap.order]
		null.corrs[r] <- cor(as.vector(A.rand),as.vector(B),method=method,use='pairwise.complete.obs')
	}
	p<-mean(null.corrs >= cor.data)
	out <- list(cor=cor.data,p=p)
	return(out)
}

#perform spatial discretization on a path (vectors x and y, which must be the same length), using a step length r
#returns new x and y coordinates (approximately r apart) and the times at which they occurred
#t0 is the initial timestep
spatial.discretization <- function(x,y,r,t0){
	i<-1
	spat.x <- spat.y <- spat.t <- c()
	ref.idx <- ref.x <- ref.y <- NA
	while(i <= length(x)){

		if(is.na(ref.x) | is.na(ref.y)){
			ref.idx <- i
			ref.x <- x[ref.idx]
			ref.y <- y[ref.idx]
			i <- i + 1
			next
		}
		
		#get current location
		x.curr <- x[i]
		y.curr <- y[i]
		
		#if current location is missing, just set ref.idx to i and continue
		if(is.na(x.curr) | is.na(y.curr)){
			ref.idx <- i
			ref.x <- x[i]
			ref.y <- y[i]
			i <- i + 1
			next
		}
		
		sq.dist <- (ref.x - x[i])^2 +(ref.y - y[i])^2
		
		if(sq.dist >= r^2){
			spat.x <- c(spat.x,x.curr)
			spat.y <- c(spat.y,y.curr)
			spat.t <- c(spat.t,i + t0 - 1)
			
			ref.x <- x.curr
			ref.y <- y.curr
			
		}
		else{
			i <- i + 1
		}
		
	}
	
	out <- list(x=spat.x,y=spat.y,t=spat.t)
	invisible(out)
	
}

correlation.randomization <- function(x,y,n.rands=1000,method='pearson',two.tailed=F){
	cor.dat <- cor(x,y,method=method)
	cor.null <- rep(NA,n.rands)
	for(i in 1:n.rands){
		cor.null[i] <- cor(x,y[sample(length(y))],method=method)
	}
	if(!two.tailed){
		p <- mean(cor.dat <= cor.null)
	}
	if(two.tailed){
		p <- mean(abs(cor.dat) >= cor.null)*2
	}
	out <- list(cor.dat=cor.dat,p=p)
	invisible(out)
	
}


















