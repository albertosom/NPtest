

source("functionsNptest.R")

#'  Test for conditional mean equal to the unconditional mean
#'  
#' @param Y a vector. The outcome variable. The length of Y equals to the number of observations.
#' @param X covariate (1 dimensional for now)
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
CondMean<-function(Y,X,boot.samp,lambdas,d=100) {
  
  n<-length(Y)
  null.dist.matrix<-matrix(NA,nrow=boot.samp,ncol=length(lambdas))
  lamb.pval<-numeric()
  pvals.aggreg<-numeric()
  rr<-matrix(NA,ncol=4,nrow=length(lambdas))
  test.sup<-numeric()
  
  ## scale the covariate to be between (0,1)
  x<-X
  x.scale0 <- x - mean(x)
  x.scale <- (x.scale0 - ((n+1)/n) * min(x.scale0))/
    (((n+1)/n)*(max(x.scale0) - min(x.scale0)))
  X<-x
  
    for(k in 1:length(lambdas)){
      
      # lambda1<-0.001
      lambda=lambdas[k]
      
      U<-SobolevBasis(X,d,n)
      U<-U[,-1] # not include intercept
      psiu<-matrix(Y-mean(Y),ncol=1)
      S.theta<-(1/n)*(t(psiu)%*%U)
      Pen<-diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))
      fit1<-loess(Y~X,control=loess.control(surface="direct"))
      out.reg1<-predict(fit1, newdata = data.frame(X=X))
      out.reg<-predict(fit1, newdata = data.frame(X=X))
      V<-var.h(d,U)
      h.hat<-S.theta%*%solve(Pen+lambda*V)
      rkhs.norm<-(h.hat)%*%diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))%*%t(h.hat)
      V.hat<-h.hat%*%V%*%t(h.hat)
      
      U<-SobolevBasis(X,d,n)
      U<-U[,-1] # not include intercept
      psiu<-matrix(Y-mean(Y),ncol=1)
      S.theta<-(1/n)*(t(psiu)%*%U)
      stat<-as.numeric(S.theta%*%t(h.hat))
      rr[k,]<-c(lambda,stat,rkhs.norm,V.hat)
      
      ### approximating asymptotic distribution under the null 
    
        for ( i in 1:boot.samp) {
          
        epsilon<-rnorm(n)
        psiu<-matrix(epsilon*(Y-mean(Y))-mean(epsilon)*out.reg,ncol=1)
        S.theta<-(1/n)*(t(psiu)%*%U)
        h.hat<-S.theta%*%solve(Pen+lambda*V)
        psiu<-matrix(epsilon*(Y-mean(Y))-mean(epsilon)*out.reg1,ncol=1)
        S.theta<-(1/n)*(t(psiu)%*%U)
        test.sup[i]<-S.theta%*%t(h.hat)
        
        }
       
      null.dist.matrix[,k]<- test.sup
      p.val<-mean(test.sup > abs(stat),na.rm=T)
      lamb.pval[k]<-p.val
  
    }
    return(list("test.stats"=rr[,2] , "p.values"=lamb.pval,"null.dstr.lambdas"=null.dist.matrix))
  }
  