
library(np)
source("functionsNptest.R")
library(MASS)
library(parallel)
library(mvtnorm)


## Example 1 testing effect of explanatory variable on outcome

n=50
A<-rbinom(n,size=1,prob=0.5)
W<-rnorm(n,0,1)
epsilon<-rnorm(n,0,1)
Y<- sin(W*sign(W)) + epsilon

summary(lm(Y~W))

## approach by Racine (2006)
bw <- npregbw(formula=Y~W,regtype="ll",bwmethod="cv.aic")
npsigtest(bws=bw,index=c(1))

## our non-parametric test 
lambdas<-seq(0,5e5,length.out=5)
rr<-CondMean(Y,W,boot.samp=500,lambdas,lambda1=0.001,d=100,ncores=detectCores() - 1) 
  

library(NHANES)

colnames(NHANES)

Chol<-NHANES$TotChol
syst<-NHANES$BPSysAve
df<-data.frame(y=Chol,x=syst)
df<-na.omit(df)
plot(y=df$y,x=df$x)
summary(lm(df$y~df$x))

df<-df[sample(1:nrow(df),size=500),]

## approach by Racine (2006)
bw <- npregbw(formula=df$y~df$x,regtype="ll",bwmethod="cv.aic")
npsigtest(bws=bw,index=c(1))

rr<-CondMean(df$y,df$x,boot.samp=800,lambdas=lambdas,lambda1=0.001,d=100,ncores=detectCores() - 1) 


## Testing null conditional average treatment effect 

n=50
A<-rbinom(n,size=1,prob=0.5)
W<-rnorm(n,0,1)
#W<-runif(n,-1,1)
epsilon<-rnorm(n,0,1)
beta=-0.5
Y<- beta*A*(1+W^2) + 0.4 + W + epsilon
#Y<-0.4 + 0.05*W + 0.8*A

## approach by Racine (2006)
bw <- npregbw(formula=Y~A+W,regtype="ll",bwmethod="cv.aic")
npsigtest(bws=bw,index=c(1))

lambdas<-seq(0,5e5,length.out=5)
rr<-testnullcate(Y,W,A,boot.samp,lambdas,lambda1=0.001,d=100,ncores=detectCores() - 1)

