
library(np)
source("functionsNptest.R")
source("CondMean.R")
source("testnullcate.R")
library(MASS)
library(parallel)
library(mvtnorm)


## Example 1 testing effect of explanatory variable on outcome

n=500
A<-rbinom(n,size=1,prob=0.5)
W<-rnorm(n,0,1)
epsilon<-rnorm(n,0,1)
Y<- sin(W*sign(W)) + epsilon

## our non-parametric test 
lambdas<-seq(0,5e5,length.out=5)
rr<-CondMean(Y,W,boot.samp=500,lambdas,d=100,ncores=detectCores() - 1) 
  


## Example 3 Testing null conditional average treatment effect 

n=50
A<-rbinom(n,size=1,prob=0.5)
W<-rnorm(n,0,1)
epsilon<-rnorm(n,0,1)
Y<- 0.5 + 0.5*W + 10*A + 1.6*W*A + epsilon
#Y<-0.4 + 0.05*W + 0.8*A +epsilon

## approach by Racine (2006)
bw <- npregbw(formula=Y~A+W,regtype="ll",bwmethod="cv.aic")
npsigtest(bws=bw,index=c(1))

lambdas<-seq(0,5e5,length.out=5)
rr<-testnullcate(Y,W,A,boot.samp=500,lambdas,d=100,ncores=detectCores() - 1)

