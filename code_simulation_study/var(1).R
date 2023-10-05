# the functions below are used to generate the data within a dataset following VAR(1)
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
### create VAR with lag 1
k=4 # number of ILDs
y0=rep(0,k) # intercept

# generate ILDs
var.sim=function(phi,sigma,y,time,p)## y is intercept
{# time-length of ILDs sigma-innovation matrix
  phi=as.vector(phi)
  dim(phi)=c(k,k,p)
  res=varima.sim(phi = phi, theta = NULL,  n =time, sigma = sigma) 
  return(res)
  
}
### set the coefficients and innovation matrix for reference and other clusters
var.con=function(k,phi1.low,phi1.high,phi2.low,phi2.high,dist,num_c,ef){
  # phi1.low: the lower bound of reference clusters AR effect
  # phi1.high: the upper bound of reference clusters AR effect
  # phi2.low: the lower bound of reference clusters CR effect
  # phi2.high: the upper bound of reference clusters CR effect
  # dist: distance between two coefficients in phi matrix
  # num_c: number of clusters
  # ef: number of effective coefficients
  while(TRUE){
    # set AR and CR effect for reference group
    phi=diag(runif(k, phi1.low, phi1.high))
    
    for (i in 1:k){
      for (j in 1:k){
        if(i!=j){
          phi[i,j]=sample(runif(1,phi2.low,phi2.high))
        }
      }
    }
    ## check the stationary
    roots = eigen(phi)$values
    roots =Mod(roots)
    if(all(roots < 1)==1){
      break
    }}
  
  sigma=1/2*diag(k)+0.5
  
  phi.1=as.vector(phi)
  ## creat phi.2 (other cluster except reference cluster)
  ### set a distance from 0.05-0.25
  ### equally assign the negative and positive sign to each entry
  ### check the stationary of phi.2
  while(TRUE){
    phi.f=c()
    phi.f=rbind(phi.f,phi.1)
    for(i in 2:num_c){
      while(TRUE){
        coef=as.vector(phi)
        
        random=sample(1:(k*k),ef,replace=F)
        sign=sample(c(-1,1),size=ef,replace=T,prob=c(0.5,0.5))
        add=dist*sign
        coef[random]=coef[random]+add
        phi.2=matrix(coef,ncol=k)
        roots = eigen(phi.2)$values
        roots =Mod(roots)
        if(all(roots < 1)==1){
          break
        }}
      
      phi.2=as.vector(phi.2)
      
      phi.f=rbind(phi.f,phi.2)}
    
    if(length(unique(as.vector(dist(phi.f,method="euclidean"))))==1){
      # check whether the totall distance equal between any pair of clusters
      break
    }}
  
  return(list(phi.f,sigma))} # phi.f: phi matrix, sigma:innovation matrix

### simulate data with lag 1
sim.data.var1=function(nob,percent,k,dist,ncluster,time,ef){
  # nob: sample size
  # percent: relative size of clusters
  # ncluster: number of cluster
  # time: length of ILDs

  ###create #n cluster
  
  mem.t=c()
  par=var.con(k,0.5,0.7,-0.4,0.4,dist,ncluster,ef)
  slope=par[[1]]
  ser=list()
  j=0
  for (class in 1:nrow(par[[1]])){
    nob.1=nob*percent[class]
    
    phi=matrix(par[[1]][class,],nrow=k)
    
    sigma=par[[2]]
    
    
    for(i in (1+j):(j+nob.1)){
      
      ser[[i]]=var.sim(phi,sigma,y,time)}
    j=i
    mem.t=c(mem.t,rep(class,nob.1))}
  return(list(ser,mem.t,slope))
} # ser: ILDs within a sample. mem.t: membership of each individual within a cluster
# slope: the phi matrix for each cluster