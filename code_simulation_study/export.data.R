
# function used for derive the coefficents from the raw ILDs for each individual
export.data=function(nob,percent,k,p,repet,dist,ncluster,ef) #mix_num only used in mixed.R
  # nob: sample size
  # percent: relative size of each cluster within a dataset
  # k: number of Intensive longitudinal variables
  # dist: distance between two cluster
  # ncluster: number of cluster
  # ef: effective coefficients
  # repet: number of replicates for each condition
  # p: temporal order for data generated model 
  #     p=1 when proportion of data following var(2)=0% 
  #     p=2 when proportion of data following var(2)=10%,40%,60%,90%,100%
{ 
  #set.seed(1)
  a=Sys.time()
  phi.t=c()
  phi.ind=c()
  phi.ind1=c()
  phi.ind2=c()
  phi.ind3=c()
  # r=c()
  # r1=c()
  # d=c()
  # d1=c()
  # optimal_order=c()
  cl1<-makeCluster(10)
  registerDoParallel(cl1) 
  clusterEvalQ(cl1, library(MASS))
  clusterEvalQ(cl1, library(portes))
  clusterEvalQ(cl1, library(BigVAR))
  clusterExport(cl1,list('sim.data.var1','var.con','var.sim'))
  clusterExport(cl1,c('repet','k','nob','percent','dist','ncluster','time','ef','y','p'))
  
  # generate data for each replicate under each condition
  data=foreach(1:repet) %dopar% {
    data= sim.data.var1(nob,percent,k,dist,ncluster,time,ef)
    return(data)}
  
  stopCluster(cl1)
  closeAllConnections() 
  print(Sys.time()-a)
  for (times in 1:repet){
    
    ### get AR and CR coefficients from VAR(p) model
    ser=data[[times]][[1]] # ILDs
    mem.t=data[[times]][[2]] # membership
    phi=data[[times]][[3]] # data generated model
    model_per_pers = matrix(0, nob, (k*k)*1) # LO method
    model_per_pers1 = matrix(0, nob, (k*k)*2) # HO method
    
    
    # get AR and CR effect for each person whithin each dataset
    for(i in 1:nob){
      # LO  fit VAR(1)
      v <- VAR(ser[[i]], type = "none", #  none means no intercept was estimated due to the intercept set to 0
               season = NULL,p=1)
      terms = length(v$varresult[[1]]$coefficients) # number of coefficients
      coeff_pers = matrix(0, k, k*1)
      
      
      ###HO fit VAR(2)
      v2 <- VAR(ser[[i]], type = "none", 
                season = NULL,p=2)
      terms2 = length(v2$varresult[[1]]$coefficients)
      coeff_pers2 = matrix(0, k, k*2)
      for(vas in 1:k){
        
        coeff_pers[vas,1:(terms) ] = v$varresult[[vas]]$coefficients[1:terms]
        
        coeff_pers2[vas,1:(terms2) ] = v2$varresult[[vas]]$coefficients[1:terms2]
      }
 
      
      model_per_pers[i, ] = as.vector(t(regrcoeff_pers))
      
      model_per_pers2[i, ] = as.vector(t(regrcoeff_pers2))
     
    }
    
    
    
    ### calculate the overlap
    coef_over=cbind(model_per_pers,mem.t)
    coef_over2=cbind(model_per_pers2,mem.t)
    
    ### calculate the overlap
    phi.ind=rbind(phi.ind,coef_over)
    phi.ind2=rbind(phi.ind2,coef_over2)
    
    phi.t=rbind(phi.t,phi)
  }
  print(Sys.time()-a)
  
  result=list(phi.ind,phi.ind2,phi.t)
  
  return(result)} 
# phi.ind coefficients for VAR(1)
# phi.ind2 coefficients for VAR(2)
# phi.t coefficients for data generated model

## test
#dd=export.data(randomstart=10,nob=60,time=50,percent=percent,k=k,p=2,repet=10,dist=0.15,ncluster=4,ef=4,mix_num=1)


