

source("functionsNptest.R")

#'  Test conditional average treatment effect equal to average treatment effect
#'  
#' @param Y a vector. The outcome variable. The length of Y equals to the number of observations.
#' @param X Covariate (1 dimensional for now)
#' @param A Treatment indicator
#' @param boot.samp Number of bootstrap samples
#' @param d Number of basis
#' @return A list containing necessary information to compute p-value of the test
#' \item{test.stats}{Test statistics for different lambda values}
#' \item{p.values}{P-values for the test for each value of lambda}
#' \item {null.dstr.lambdas} {Multiplier bootstrap sample for approximating the distribution
#'                             of the test statistic for different lambda values.}
#' 
#' @export
#' 
testnullcate<-function(Y,W,A,boot.samp,lambdas,d=100) {
  
  n<-length(Y)
  null.dist.matrix<-matrix(NA,nrow=boot.samp,ncol=length(lambdas))
  lamb.pval<-numeric()
  pvals.aggreg<-numeric()
  rr<-matrix(NA,ncol=4,nrow=length(lambdas))
  test.sup<-numeric()
  
  ## scale the covariate to be between (0,1)
  x<-W
  x.scale0 <- x - mean(x)
  x.scale <- (x.scale0 - ((n+1)/n) * min(x.scale0))/
    (((n+1)/n)*(max(x.scale0) - min(x.scale0)))
  W<-x
  
  for(k in 1:length(lambdas)){
    
    lambda=lambdas[k]
    
    fit_or <- glm (Y ~ ., data = data.frame(A,W))
    fit_ps<- glm(A ~ .,data=data.frame(W),family =binomial())
    # probability of receiving treatment
    g1W <- predict(fit_ps,newdata=data.frame(W),type ="response")
    g1W<-0.5
    # predict on data setting A=1
    Qbar1W <- predict(fit_or,newdata=data.frame(W,A=1))
    # predict on data setting A=0
    Qbar0W <- predict(fit_or ,newdata=data.frame(W,A=0))
    
    U<-SobolevBasis(W,d,n)
    U<-U[,-1] # not include intercept
    #A=1
    psiu1<-matrix((1/g1W)*(Y-Qbar1W)-mean(Qbar1W)+Qbar1W,ncol=1)
    S.theta1<-(1/n)*(t(psiu1)%*%U)
    #A=0
    psiu0<-matrix((1/(1-g1W))*(Y-Qbar0W)-mean(Qbar0W)+Qbar0W,ncol=1)
    S.theta0<-(1/n)*(t(psiu0)%*%U)
    
    S.theta<-S.theta1-S.theta0
    Pen<-diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))
    V<-var.h(d,U)
    h.hat<-S.theta%*%solve(Pen+lambda*V)
    rkhs.norm<-(h.hat)%*%Pen%*%t(h.hat)
    V.hat<-h.hat%*%V%*%t(h.hat)
    stat<-as.numeric(S.theta%*%t(h.hat))
    
    #A=1
    psiu1<-matrix((1/g1W)*(Y-Qbar1W)-mean(Qbar1W)+Qbar1W,ncol=1)
    S.theta1<-(1/n)*(t(psiu1)%*%U)
    #A=0
    psiu0<-matrix((1/(1-g1W))*(Y-Qbar0W)-mean(Qbar0W)+Qbar0W,ncol=1)
    S.theta0<-(1/n)*(t(psiu0)%*%U)
    S.theta<-S.theta1-S.theta0
    stat<-as.numeric(S.theta%*%t(h.hat))
    
    rr[k,]<-c(lambda,stat,rkhs.norm,V.hat)
    
    ### approximating asymptotic distribution under the null 
    
    for(i in 1:boot.samp) {
      epsilon<-rnorm(n)
      #A=1
      psiu1<-matrix(epsilon*((1/g1W)*(Y-Qbar1W)-mean(Qbar1W)+Qbar1W)-mean(epsilon)*Qbar1W,ncol=1)
      S.theta1<-(1/n)*(t(psiu1)%*%U)
      #A=0
      psiu0<-matrix(epsilon*(1/(1-g1W)*(Y-Qbar0W)-mean(Qbar0W)+Qbar0W)-mean(epsilon)*Qbar0W,ncol=1)
      S.theta0<-(1/n)*(t(psiu0)%*%U)
      S.theta<-S.theta1-S.theta0
      
      psiu<-psiu1-psiu0
      h.hat<-S.theta%*%solve(Pen+lambda*V)
  
      #A=1
      psiu1<-matrix(epsilon*((1/g1W)*(Y-Qbar1W)-mean(Qbar1W)+Qbar1W)-mean(epsilon)*Qbar1W,ncol=1)
      S.theta1<-(1/n)*(t(psiu1)%*%U)
      #A=0
      psiu0<-matrix(epsilon*(1/(1-g1W)*(Y-Qbar0W)-mean(Qbar0W)+Qbar0W)-mean(epsilon)*Qbar0W,ncol=1)
      S.theta0<-(1/n)*(t(psiu0)%*%U)
      
      S.theta<-S.theta1-S.theta0
      h.hat<-S.theta%*%solve(Pen+lambda*V)
      test.sup[i]<-S.theta%*%t(h.hat)
      
    }
    
    null.dist.matrix[,k]<- test.sup
    p.val<-mean(test.sup > stat,na.rm=T)
    lamb.pval[k]<-p.val
  }
  return(list("test.stats"=rr[,2] ,"p.values"=lamb.pval,"null.dstr.lambdas"=null.dist.matrix))
}
