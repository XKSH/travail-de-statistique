#Q1 and Q2
betapara<-matrix(c(0.2,0.4,0.5,1,2,2,0.5,2,0.5,1,2,6),nrow=6,ncol=2)#initialize the parameters
theta_hatm<-matrix(0,nrow=6,ncol=2)#initialize the moments estimators
#MNS calculate the moments
MNS<-function(x){
  xmean<-mean(x)
  #sample variance
  xvar<-var(x)
  #compute moment estimator alpha and beta
  alpham=xmean*(xmean*(1-xmean)/xvar-1)
  betam=(1-xmean)*(xmean*(1-xmean)/xvar-1)
  theta_hatm<-c(alpham,betam)
  return(theta_hatm)
}
#calculate the order 1 derivative
DL<-function(a,x){
  dalpha<-digamma(sum(a))-digamma(a[1])+mean(log(x))
  dbeta<-digamma(sum(a))-digamma(a[2])+mean(log(1-x))
  dl<-c(dalpha,dbeta)
  return(dl)
}
#calculate the order 2 derivative
DDL<-function(a){
  ddalpha<-trigamma(sum(a))-trigamma(a[1])
  ddbeta<-trigamma(sum(a))-trigamma(a[2])
  ddab<-trigamma(sum(a))
  ddl<-matrix(c(ddalpha,ddab,ddab,ddbeta),nrow=2)
  return(ddl)
}
ns<-6 #define the number of the groups to be compared
mns<-matrix(rep(0,2*ns),nrow=2,ncol=ns*6)#initialize the vector of moment estimator
mlens<-matrix(rep(0,2*ns),nrow=2,ncol=ns*6)#initialize the vector of moment estimator
tn<-50*c(1:ns)#use 50 as one unit of the number of the samples
for(num in 1:6){
  for (i in 1:ns){
    set.seed(50)
    x<-rbeta(tn[i], betapara[num,1],betapara[num,2] )#create the simulated data
    mns[,ns*num-ns+i]<-MNS(x)
    theta_hat<-MNS(x) #initialize the original point
    ddL<-DDL(theta_hat)
    dL<-DL(theta_hat,x)
    n<-0
    epi=1 #initialize the episilon
    #update the theta_hat
    while(epi>10^(-8)){
      #record the old theta_hat
      old<-theta_hat
      #update theta_hat
      theta_hat<-theta_hat-solve(ddL)%*%dL
      #update dL
      dL<-DL(theta_hat,x)
      #updata ddL
      ddL<-DDL(theta_hat)
      epi<-sum((theta_hat-old)^2)
      n<-n+1
    }
    mlens[,ns*num-ns+i]=theta_hat
  }
}
#Q3:the graphs for the different numbers of sample,for the simplicity,just use moment estimator as the original point
#the function calculate the moment estimators
MNS<-function(x){
  xmean<-mean(x)
  #sample variance
  xvar<-var(x)
  #compute moment estimator alpha and beta
  alpham=xmean*(xmean*(1-xmean)/xvar-1)
  betam=(1-xmean)*(xmean*(1-xmean)/xvar-1)
  thetachapeaum<-c(alpham,betam)
  return(thetachapeaum)
}
ns<-15 #define the number of the groups to be compared
mns<-matrix(rep(0,2*ns),nrow=2,ncol=ns)#initialize the vector of moment estimator
mlens<-matrix(rep(0,2*ns),nrow=2,ncol=ns)#initialize the vector of moment estimator
tn<-20*c(1:ns)#use 50 as one unit of the number of the samples
#the iteration of MLEs
for (i in 1:ns){
  set.seed(50)
  x<-rbeta(tn[i], 6, 9)
  mns[,i]<-MNS(x)
  thetachapeau<-MNS(x) #initialize the original point
  ddL<-DDL(thetachapeau)
  dL<-DL(thetachapeau,x)
  n<-0
  epi=1
  while(epi>10^(-8)){
    #record the old thetachapeau
    old<-thetachapeau
    #update thetachapeau
    thetachapeau<-thetachapeau-solve(ddL)%*%dL
    #update dL
    dL<-DL(thetachapeau,x)
    #updata ddL
    ddL<-DDL(thetachapeau)
    epi<-sum((thetachapeau-old)^2)
    n<-n+1
  }
  mlens[,i]=thetachapeau
}
#the plot precedure
plot(tn, mns[1,], type="b",
     pch=6, lty=1, col="red", ylim=c(1, 3),
     main="Moment estimator vs. MLE by N-R",
     xlab="number of samples", ylab="alpha")
lines(tn, mlens[1,], type="b",
      pch=8, lty=2, col="blue")
abline(h=2)
legend("topright", inset=.05, title="Estimator Type", c("ME","MLE"),
       lty=c(1, 2), pch=c(15, 17), col=c("red", "blue"))
plot(tn, mns[2,], type="b",
     pch=6, lty=1, col="red", ylim=c(5, 10),
     main="Moment estimator vs. MLE by N-R",
     xlab="number of samples", ylab="beta")
lines(tn, mlens[2,], type="b",
      pch=8, lty=2, col="blue")
abline(h=6)
legend("topright", inset=.05, title="Estimator Type", c("ME","MLE"),
       lty=c(1, 2), pch=c(15, 17), col=c("red", "blue"))
#pvalue and pvalueb is the pvalue according to chisq distribution and bootstrap respectively
#define the dL ddL for alpha-beta,use alpha-beta and beta as the parameters to  estimate
DLw<-function(a,x){
  dab<-digamma(sum(a))-digamma(sum(a)+a[2])-mean(log(x))
  dbeta<-digamma(sum(a))+digamma(a[2])-2*digamma(sum(a)+a[2])-mean(log(1-x))-mean(log(x))
  dl<-c(dab,dbeta)
  return(dl)
}
DDLw<-function(a){
  ddamb<-trigamma(sum(a))-trigamma(sum(a)+a[2])
  ddbeta<-trigamma(sum(a))+trigamma(a[2])-4*trigamma(sum(a)+a[2])
  ddab<-trigamma(sum(a))-2*trigamma(sum(a)+a[2])
  ddl<-matrix(c(ddamb,ddab,ddab,ddbeta),nrow=2)
  return(ddl)
}
ns<-15 #define the number of the groups to be compared
wald<-rep(0,ns)#initialize the vector of wald test statistics
pvalue<-rep(0,ns)
mlensw<-matrix(rep(0,2*ns),nrow=2,ncol=ns)#initialize mle
tn<-20*c(1:ns)#use 20 as one unit of the number of the samples
for (i in 1:ns){
  set.seed(50)
  x<-rbeta(tn[i],0.5, 0.5) #change parameter here
  #initialize the original point
  th<-MNS(x)
  th[1]<-th[1]-th[2]
  thw<-th
  ddL<-DDLw(thw)
  dL<-DLw(thw,x)
  n<-0
  epi=1
  while(epi>10^(-8)){
    #record the old thw
    old<-thw
    #update thw
    thw<-thw-solve(ddL)%*%dL
    #update dL
    dL<-DLw(thw,x)
    #updata ddL
    ddL<-DDLw(thw)
    epi<-sum((thw-old)^2)
    n<-n+1
  }
  mlensw[,i]=thw
  iddL<-solve(tn[i]*ddL)
  #calculate the wald satistic
  wald[i]<-mlensw[1,i]^2/iddL[1,1]
  #calculate the pvalue
  pvalue[i]<-2*pchisq(wald[i],1)*(pchisq(wald[i],1)<0.5)+2*(1-pchisq(wald[i],1))*(pchisq(wald[i],1)>0.5)
}
# parametric bootstrap 
mlensb<-c(0,0)#initialize mle
B=1000
pvalueb<-rep(0,ns)
for (i in 1:ns){
  theta1<-rep(0,B)
  for(j in 1:B){
    set.seed(30+10*j)
    xb<-rbeta(tn[i],mlensw[2,i] ,mlensw[2,i])#resampling from sampling distribution obeying the null hypothesis
    th<-MNS(xb)
    th[1]<-th[1]-th[2]
    thw<-th
    ddL<-DDLw(thw)
    dL<-DLw(thw,xb)
    epi=1
    while(epi>10^(-8)){
      #record the old thw
      old<-thw
      #update thw
      thw<-thw-solve(ddL,tol=1e-20)%*%dL
      #update dL
      dL<-DLw(thw,xb)
      #updata ddL
      ddL<-DDLw(thw)
      epi<-sum((thw-old)^2)
    }
    mlensb=thw
    iddL<-solve(tn[i]*ddL,tol=1e-20)
    theta1[j]<-mlensb[1]
  }
  p<-sum(abs(theta1)>abs(mlensw[1,i]))/B#calculate the p value
  pvalueb[i]<-p
}