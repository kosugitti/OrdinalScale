#' DS for multiple choice
#'
#' @param dat A matrix or data frame
#' @param eps A criteria whether it reached conversion. Default value is 10e-5
#' @param iter.max A muximum numbers of iterations
#' @export

DS_multi <- function(dat,eps=10e-5,iter.max=1000){
  dat.dm <- make.dummy(dat)
  ds <- DualScaling(dat.dm,eps=eps,iter.max=iter.max)

  # Scoring Table
  ScoringTable <- array(NA,dim=c(NROW(dat),NCOL(dat),ds$ndims))
  opt <- apply(dat,2,max)
  opt <- cumsum(opt)-opt[1]
  for(d in 1:ds$ndims){
    for(i in 1:NCOL(dat)){
      for(j in 1:NROW(dat)){
        ScoringTable[j,i,d] = ds$NormedCol[opt[i]+dat[j,i],d]
      }
    }
  }

  ## Informations with each Dimensions
  SSj <- matrix(ncol=NCOL(dat),nrow=ds$ndims)
  RTj <- matrix(ncol=NCOL(dat),nrow=ds$ndims)
  RTj2 <- RTj
  alpha <- rep(NA,ds$ndims)
  CorrTable <- array(dim=c(NCOL(dat),NCOL(dat),ds$ndims))
  for(d in 1:ds$ndims){
    SSj[d,] <- apply(ScoringTable[,,d],2,function(x) sum(x^2))
    sumXit <- apply(ScoringTable[,,d],1,sum)
    RTj[d,] <- cor(ScoringTable[,,d],sumXit)
    RTj2[d,] <- cor(ScoringTable[,,d],sumXit)^2
    alpha[d] <- NCOL(dat)/(NCOL(dat)-1) * ((sum(RTj2[d,])-1)/sum(RTj2[d,]))
    CorrTable[,,] <- cor(ScoringTable[,,d])
  }

  ## for Output
  res = rbind(ds$result,alpha)
  dimName <- paste0("Dim", 1:ds$ndims)
  itemName <- paste0("Item",1:NCOL(dat))
  rownames(SSj) <- dimName
  rownames(RTj) <- dimName
  rownames(RTj2) <- dimName
  colnames(SSj) <- itemName
  colnames(RTj) <- itemName
  colnames(RTj2) <- itemName


  return(list(
    result = res,
    SSj = t(SSj),
    RTj = t(RTj),
    RTj2 = t(RTj2),
    ndims = ds$ndims,
    singularValue = ds$singularValue,
    eigenValue = ds$eigenValue,
    delta = ds$delta,
    cumulativeDelta = ds$cumulativeDelta,
    NormedCol = ds$NormedCol,
    NormedRow = ds$NormedRow,
    ProjectedCol = ds$ProjectedCol,
    ProjectedRow = ds$ProjectedRow
  ))
}


