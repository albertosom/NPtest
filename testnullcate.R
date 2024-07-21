
## testing treatment effect heterogeneity 
library(parallel)
library(mvtnorm)
source("functionsNptest.R")

testnullcate<-function(Y,W,A,boot.samp,lambdas,lambda1=0.001,d=100,ncores=detectCores() - 1) {
  
  n<-length(Y)
  paths<-rmvnorm(boot.samp, sigma=diag(1,nrow=n))
  null.dist.matrix<-matrix(NA,nrow=boot.samp,ncol=length(lambdas))
  lamb.pval<-numeric()
  pvals.aggreg<-numeric()
  rr<-matrix(NA,ncol=4,nrow=length(lambdas))
  
  for(k in 1:length(lambdas)){
    
    # lambda1<-0.001
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
    
    U<-SobBasis(W,d,n)
    U<-U[,-1]
    #A=1
    psiu1<-matrix((1/g1W)*(Y-Qbar1W)-mean(Qbar1W)+Qbar1W,ncol=1)
    S.theta1<-(1/n)*(t(psiu1)%*%U)
    #A=0
    psiu0<-matrix((1/(1-g1W))*(Y-Qbar0W)-mean(Qbar0W)+Qbar0W,ncol=1)
    S.theta0<-(1/n)*(t(psiu0)%*%U)
    
    S.theta<-S.theta1-S.theta0
    Pen<-diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))
    V<-var.h(d,U)
    h.hat<-(1/lambda1)*S.theta%*%solve(Pen+lambda*V)
    rkhs.norm<-(h.hat)%*%diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))%*%t(h.hat)
    V.hat<-h.hat%*%V%*%t(h.hat)
    stat<-as.numeric(S.theta%*%t(h.hat))
    
    U<-SobBasis(W,d,n)
    U<-U[,-1]
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
    
    stats <- mclapply(1:boot.samp, function(x0) {
      
      epsilon<-paths[x0,]
    
      U<-SobBasis(W,d,n)
      U<-U[,-1]
      #A=1
      psiu1<-matrix(epsilon*((1/g1W)*(Y-Qbar1W)-mean(Qbar1W)+Qbar1W)-mean(epsilon)*Qbar1W,ncol=1)
      S.theta1<-(1/n)*(t(psiu1)%*%U)
      #A=0
      psiu0<-matrix(epsilon*(1/(1-g1W)*(Y-Qbar0W)-mean(Qbar0W)+Qbar0W)-mean(epsilon)*Qbar0W,ncol=1)
      S.theta0<-(1/n)*(t(psiu0)%*%U)
      S.theta<-S.theta1-S.theta0
      
      psiu<-psiu1-psiu0
      Pen<-diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))
      h.hat<-(1/lambda1)*S.theta%*%solve(Pen+lambda*V)
      
      U<-SobBasis(W,d,n)
      U<-U[,-1]
      #A=1
      psiu1<-matrix(epsilon*((1/g1W)*(Y-Qbar1W)-mean(Qbar1W)+Qbar1W)-mean(epsilon)*Qbar1W,ncol=1)
      S.theta1<-(1/n)*(t(psiu1)%*%U)
      #A=0
      psiu0<-matrix(epsilon*(1/(1-g1W)*(Y-Qbar0W)-mean(Qbar0W)+Qbar0W)-mean(epsilon)*Qbar0W,ncol=1)
      S.theta0<-(1/n)*(t(psiu0)%*%U)
      
      S.theta<-S.theta1-S.theta0
      Pen<-diag(c(1/(2*pi*seq(1,d,length.out=d))^(-4)))
      h.hat<-(1/lambda1)*S.theta%*%solve(Pen+lambda*V)
      
      test.sup<-S.theta%*%t(h.hat)
      
      return(test.sup)
    } ,mc.cores = ncores)
    
    null.dist.matrix[,k]<- unlist(stats)
    p.val<-mean(abs(unlist(stats)) > abs(stat),na.rm=T)
    lamb.pval[k]<-p.val
  }
  
  ## aggregating the test statistic approach
  
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
