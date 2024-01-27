load("~/times.Rdata")
load("~/x_axis.Rdata")
load("~/y_axis.Rdata") 
ids <- rownames(xs)

'%!in%' <- function(x,y)!('%in%'(x,y))

m <- c()
b <- list()
for (i in 1: length(ids)){
  
  a <- rle(is.na(xs[i,]))
  b[[i]] <- rep(a$values*a$lengths,a$lengths)
  
  xs[i, which(b[[i]] == 1)] <-  (xs[i, which(b[[i]] == 1)+ 1] + xs[i, which(b[[i]] == 1)-1])/2
  ys[i, which(b[[i]] == 1)] <-  (ys[i, which(b[[i]] == 1)+ 1] + ys[i, which(b[[i]] == 1)-1])/2
  m <- c(m, length(which(b[[i]] == 1)))
}

sum(m)*100/ sum(!is.na(xs))
# xx% of the overall dataset is comprised of interpolated missing values

#save(xs, file="~/xs.Rdata")
#save(ys, file="~/ys.Rdata")
