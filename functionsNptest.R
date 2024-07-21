
## basis expansion
SobBasis <- function(x, d,n) {
  Phi<-matrix(0,ncol=d+1,nrow=n)
  n <- length(x)
  x.scale0 <- x - mean(x)
  x.scale <- (x.scale0 - ((n+1)/n) * min(x.scale0))/
    (((n+1)/n)*(max(x.scale0) - min(x.scale0)))
  H <- matrix(0, nrow = n, ncol = d + 1) #eigen basis for H
  H[,1] <- x.scale
  for(i in 1:(d/2)) {
    H[,2*i]<- 2^(1/2)*cos(2*i*pi*(x.scale))
    H[,2*i+1]<- 2^(1/2)*sin(2*i*pi*(x.scale))
  }
  
  for( i in 1:n) {
    x0.scale <- x[i] - mean(x)
    x0.scale <- (x0.scale - ((n+1)/n) * min(x.scale0))/
      (((n+1)/n) * (max(x.scale0) - min(x.scale0)))
    
    phi <- numeric(d+1)
    phi[1] <- x0.scale
    phi[1 + seq(1, d, 2)] <- 2^(1/2)*cos(2*(1:(d/2))*pi*(x0.scale))
    phi[2 + seq(1, d, 2)] <- 2^(1/2)*sin(2*(1:(d/2))*pi*(x0.scale))
    
    Phi[i,]<-phi
  }
  return(Phi - apply(H, 2, mean))
}

## Variance of h(x)
var.h<-function(d,U){
  V<-matrix(0,ncol=d,nrow=d)
  for(i in 1:d) {
    for(j in 1:d) {
      V[i,j]<- mean(U[,i]*U[,j])
    }
  }
  return(V)
}


## cauchy combination test
CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  if(sum(is.na(pvals)) > 0){
    stop("Cannot have NAs in the p-values!")
  }
  
  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }
  
  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  
  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }
  
  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }
  
  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}

