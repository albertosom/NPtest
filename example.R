
## simulation to compare with Racine (2006)
install.packages("np")
library(np)
n=50
A<-rbinom(n,size=1,prob=0.5)
W<-rnorm(n,0,1)
epsilon<-rnorm(n,0,1)
beta=0
Y<- - beta*A*(1+W^2)+ sin(W*sign(W)) + epsilon

bw <- npregbw(formula=Y~W,regtype="ll",bwmethod="cv.aic")

npsigtest(bws=bw,index=c(1))

lambdas<-seq(0,5e5,length.out=5)
rr<-CondMean(Y,W,boot.samp=500,lambdas,lambda1=0.001,d=100,ncores=detectCores() - 1) 
  
