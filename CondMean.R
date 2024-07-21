
## conditional mean test 
source("functionsNptest.R")
library(MASS)
library(parallel)
library(mvtnorm)

CondMean<-function(Y,X,boot.samp,lambdas,lambda1=0.001,d=100,ncores=detectCores() - 1) {
  
  n<-length(Y)
  paths<-rmvnorm(boot.samp, sigma=diag(1,nrow=n))
  null.dist.matrix<-matrix(NA,nrow=boot.samp,ncol=length(lambdas))
  lamb.pval<-numeric()
  pvals.aggreg<-numeric()
  rr<-matrix(NA,ncol=4,nrow=length(lambdas))
  
    for(k in 1:length(lambdas)){
      
      # lambda1<-0.001
      lambda=lambdas[k]
      
      U<-SobBasis(X,d,n)
      U<-U[,-1]
      psiu<-matrix(Y-mean(Y),ncol=1)
      S.theta<-(1/n)*(t(psiu)%*%U)
      Pen<-diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))
      fit1<-loess(Y~X,control=loess.control(surface="direct"))
      out.reg1<-predict(fit1, newdata = data.frame(X=X))
      out.reg<-predict(fit1, newdata = data.frame(X=X))
      V<-var.h(d,U)
      h.hat<-(1/lambda1)*S.theta%*%solve(Pen+lambda*V)
      rkhs.norm<-(h.hat)%*%diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))%*%t(h.hat)
      V.hat<-h.hat%*%V%*%t(h.hat)
      
      U<-SobBasis(X,d,n)
      U<-U[,-1]
      psiu<-matrix(Y-mean(Y),ncol=1)
      S.theta<-(1/n)*(t(psiu)%*%U)
      stat<-as.numeric(S.theta%*%t(h.hat))
      rr[k,]<-c(lambda,stat,rkhs.norm,V.hat)
      
      ### approximating asymptotic distribution under the null 
      
      stats <- mclapply(1:boot.samp, function(x0) {
        
        epsilon<-paths[x0,]
        psiu<-matrix(epsilon*(Y-mean(Y))-mean(epsilon)*out.reg,ncol=1)
        U<-SobBasis(X,d,n)
        U<-U[,-1]
        S.theta<-(1/n)*(t(psiu)%*%U)
        Pen<-diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))
        
        h.hat<-(1/lambda1)*S.theta%*%solve(Pen+lambda*V)
        U<-SobBasis(X,d,n)
        U<-U[,-1]
        psiu<-matrix(epsilon*(Y-mean(Y))-mean(epsilon)*out.reg1,ncol=1)
        S.theta<-(1/n)*(t(psiu)%*%U)
        test.sup<-S.theta%*%t(h.hat)
        
        return(test.sup)
      } ,mc.cores = ncores)
      
      null.dist.matrix[,k]<- unlist(stats)
      p.val<-mean(abs(unlist(stats)) > abs(stat),na.rm=T)
      lamb.pval[k]<-p.val
    }
    
    ## aggregating the test statistic
    
    for(b in 1:boot.samp ){
      mu<-apply(null.dist.matrix[-b,],2,mean)
      Ts.b<-null.dist.matrix[b,]
      sigma<- sqrt(apply(null.dist.matrix[-b,],2,var))
      pvals.aggreg[b]<-max(abs(( mu- Ts.b)/sigma))
    }
    
    Q.b<-pvals.aggreg
    
    mu<-apply(null.dist.matrix,2,mean)
    sigma<- sqrt(apply(null.dist.matrix,2,var))
    Q.0<- max(abs(mu-rr[,2])/sigma)
    
    aggreg.pvalue<- (1/(boot.samp+1))*(1+ sum(ifelse(Q.b>Q.0,T,F)))
    aggreg.Pvalue<-aggreg.pvalue     
    
    ## cauchy combination test
    cct.pvalue<-CCT(lamb.pval,weights=rep(1/length(lambdas),length(lambdas)))
    p.value.cauchy<-cct.pvalue

    return(list("p-value:aggregate"=aggreg.Pvalue,"p-value:cauchy"=p.value.cauchy))
  }
  