# the functions below are used to generate the data within a dataset following VAR(2)
library(aricode)
library(dplyr)
library(MASS)
library(ggplot2)
library(lattice)
library(caret)
library(MixSim)
library(sandwich)
library(urca)
library(lmtest)
library(vars)
library(mclust)
library(svMisc)
library(iterators)
library(foreach)
library(doParallel)
library(NbClust)
library(portes)
library(parallel)
library(portes)


# generate ILDs
var.sim=function(phi,sigma,y,time,p)## y is intercept
{
  phi=as.vector(phi)
  dim(phi)=c(k,k,p)
  res=varima.sim(phi = phi, theta = NULL,  n =time, sigma = sigma)
  return(res)
  
}
# A VAR(p)-process is stable
varroot <- function (phi, k, p) 
{
  
  # all eigenvalues of the companion matrix A have modulus less than 1
  
  companion =matrix(0, nrow = k * p, ncol = k * p)
  companion[1:k, 1:(k * p)] = phi
  if (p > 1) {
    j <- 0
    for (i in (k + 1):(k * p)) {
      j = j + 1
      companion[i, j] = 1
    }
  }
  roots = eigen(companion)$values
  roots = Mod(roots)
  return(all(roots < 1))
  # returns TRUE if VAR(p) is stationary
}       

### set the coefficients and innovation matrix
var.con=function(k,phi1.low,phi1.high,phi2.low,phi2.high,phi3.low, phi3.high,dist,num_c,ef,p){
  # phi1.low: the lower bound of reference clusters AR effect in lag-1
  # phi1.high: the upper bound of reference clusters AR effect in lag-1
  # phi2.low: the lower bound of reference clusters CR effect in lag-1
  # phi2.high: the upper bound of reference clusters CR effect in lag-1
  # phi3.low: the lower bound of reference clusters AR and CR effect in lag-2
  # phi3.high: the upper bound of reference clusters AR and CR effect in lag-2
  # dist: distance between two coefficients in phi matrix
  # num_c: number of clusters
  # ef: number of effective coefficients. the effective coefficients distribute in lag-1 and lag-2 evenly
  a=Sys.time()
  while(TRUE){
    phi.lag1=diag(runif(k, phi1.low, phi1.high))
    
    for (i in 1:k){
      for (j in 1:k){
        if(i!=j){
          phi.lag1[i,j]=sample(runif(1,phi2.low,phi2.high))
        }
      }
    }
    
    phi.lag2=diag(runif(k, phi3.low, phi3.high))
    for (i in 1:k){
      for (j in 1:k){
        if(i!=j){
          phi.lag2[i,j]=sample(runif(1,phi3.low,phi3.high))
        }
      }}
    
    phi=cbind(phi.lag1,phi.lag2)
    ## check the stationary
    
    if(varroot(phi,k,p)){
      break
    }}
  
  sigma=1/2*diag(k)+0.5
  
  phi.1=as.vector(phi)
  print(Sys.time()-a)
  ## creat phi.2
  ### set a distance from 0.05-0.25
  ### equally assign the negative and positive sign to each entry
  ### check the stationary of phi.2
  
  
  while(TRUE){
    while(TRUE){
      phi.f1=c()
      phi.p=phi.1[1:(k*k)]
      phi.f1=rbind(phi.f1,phi.p)
      for(i in 2:num_c){
        
        coef=phi.p
        # cluster differences occur in lag-1 
        random=sample(1:(k*k),ef/2,replace=F)
        sign=sample(c(-1,1),size=ef/2,replace=T,prob=c(0.5,0.5))
        add=dist*sign
        coef[random]=coef[random]+add
        phi.2=matrix(coef,ncol=k)
        
        phi.2=as.vector(phi.2)
        
        phi.f1=rbind(phi.f1,phi.2)}
      
      if(length(unique(as.vector(dist(phi.f1,method="euclidean"))))==1){
        break
      }}
    print(Sys.time()-a)
    ## cluster differences occur in lag-2
    while(TRUE){
      phi.f2=c()
      phi.p=phi.1[(k*k+1):(k*k*2)]
      phi.f2=rbind(phi.f2,phi.p)
      
      for(i in 2:num_c){
        coef=phi.p
        
        random=sample(1:(k*k),ef/2,replace=F)
        sign=sample(c(-1,1),size=ef/2,replace=T,prob=c(0.5,0.5))
        add=dist*sign
        coef[random]=coef[random]+add
        phi.2=matrix(coef,ncol=k)
        
        phi.2=as.vector(phi.2)
        
        phi.f2=rbind(phi.f2,phi.2)}
      
      if(length(unique(as.vector(dist(phi.f2,method="euclidean"))))==1){
        # check the whether the distance between any pair of clusters equally
        break
      }}
    error=0
    phi.f=cbind(phi.f1,phi.f2)
    for(test in 2:num_c){
      if(varroot(phi.f[test,],k,p)==FALSE){
        error=1+error
      }}
    if(error==0){
      break}}
  
  print(Sys.time()-a)
  return(list(phi.f,sigma))} # return phi matrix for each cluster and innovation matrix

sim.data.var2=function(nob,percent,k,dist,ncluster,time,ef,p){
  # nob: sample size
  # percent: relative size of clusters
  # ncluster: number of cluster
  # time: length of ILDs
  
  ###create #n cluster
 
  mem.t=c()
  par=var.con(k,0.5,0.7,-0.4,0.4,-0.2,0.2,dist,ncluster,ef,p)
  slope=par[[1]]
  ser=list()
  j=0
  for (class in 1:nrow(par[[1]])){
    nob.1=nob*percent[class]
    
    phi=matrix(par[[1]][class,],nrow=k)
    
    sigma=par[[2]]
    
    
    for(i in (1+j):(j+nob.1)){
      
      ser[[i]]=var.sim(phi,sigma,y,time,p)}
    j=i
    mem.t=c(mem.t,rep(class,nob.1))}
  return(list(ser,mem.t,slope))
}