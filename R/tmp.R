ord.data <- matrix(c(3,2,3,1,2,
                     6,7,5,8,7,
                     8,5,7,6,4,
                     1,3,2,3,5,
                     4,1,1,2,1,
                     5,8,6,7,8,
                     2,4,4,5,3,
                     7,6,8,4,6),nrow=5)

Height.Blood <- matrix(c(2,1,2,1,4,0,1,1,3),nrow=3)
colnames(Height.Blood) <- c("High","Mid","Short")
rownames(Height.Blood) <- c("Hi-pre","Normal","Low-pre")



source("R/make_dominance.R")
d.data <- make.dominance(ord.data)
dat.tmp <- d.data

#dat <- Height.Blood
#dat.tmp <- dat




# e <- matrix(c(1,-3,3,-1,
#               -2,0,3,-1,
#               -3,-1,3,1,
#               0,-1,1,0,
#               -1,-2,3,0,
#               1,1,-1,-1,
#               1,1,1,-3,
#               -1,-3,1,3),ncol=4,byrow=T)
#
# dat.tmp <- e
# dat.tmp <- dat
# #
# dat.tmp <- Height.Blood
#
MRA <- function(dat.tmp){
  dat.tmp <- make.dominance(ord.data)
  nc <- ncol(dat.tmp)   # n
  nr <- nrow(dat.tmp)   # N
  mg.row <- rep(nc * (nc-1),nrow(dat.tmp))
  mg.col <- rep(nr * (nc-1),ncol(dat.tmp))
  mg.t   <- nr * nc * (nc-1)
  u <- rnorm(nc)
  tmp <- 0
  for( ppp in 1:1000){
    #algorithm
    v <- u %*% t(dat.tmp)
    v <- v / mg.row
    av <- (mg.row %*% t(v))/mg.t
    #      v <- v - (av * rep(1,nr))
    gy <- max(abs(v))
    v <- v/gy

    u <- v %*% dat.tmp
    u <- u / mg.col
    av <- (mg.col %*% t(u))/mg.t
    #      u <- u - (av *rep(1,nc))
    gx <- max(abs(u))
    u <- u/gx
  }

  eta2 <- gx * gy  # Correlation ratio
  return(list(eta2,u,v))
}

nc <- ncol(dat.tmp)
nr <- nrow(dat.tmp)
Hn <- t(dat.tmp)%*%dat.tmp/(nr*nc*(nc-1)^2)
eigen(Hn)$values[1]
MRA(dat.tmp)^2
MRA(Height.Blood)
