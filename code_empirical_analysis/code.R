# The code is used to do model-based clustering analysis for empirical 
# data about Positive affect and negative affect.
# And the code to create figures and tables for description the clustering solution 
# and evaluate the clustering solution
library(haven)
library(dplyr)
library(Amelia)
library(mclust) # 
library(vars)
library(factoextra)
library(cluster)
library(TSstudio)
library(gridExtra)
library(grid)
load(data)
names(data)
# calculate the pa and na from their corresponding 6 items
data=
  data%>%
  rowwise()%>%
  mutate(
    pad=mean(mad_diary_5,mad_diary_11,mad_diary_13,na.rm=T),
    nad=mean(mad_diary_6,mad_diary_14,mad_diary_16,na.rm=T),
    paa=mean(mad_diary_7,mad_diary_9,mad_diary_15,na.rm=T),
    naa=mean(mad_diary_8,mad_diary_10,mad_diary_12,na.rm=T),
    pa=mean(pad,paa,na.rm=T),
    na=mean(naa,nad.na.rm=T)
  )%>%
  mutate(event=case_when(
    morning==1~1,
    afternoon==1~2,
    evening==1~3
  ))%>%
  dplyr::select(id,time,pa,na,meting,Real_period,event)
# find id number of measurment larger than 81
ds_id=
  data%>%
  dplyr::mutate(pam =sum(!is.na(pa)), # calculate the number of measurement which is missing for each person
                nam=sum(!is.na(na)))%>%
  group_by(id)%>%
  summarise(spa=sum(pam),sna=sum(nam))%>%
  
  dplyr::filter(spa>=81& sna>=81)

# select the data
ds=data%>%
  filter(id %in% ds_id$id)

# find the person who finish over 30 days
l9=ds%>%
  filter(time>90)%>%
  dplyr::select(id)
l9=unique(l9$id)  


ds=ds%>%
  group_by(id,event)%>%
  mutate(n=1:n())%>%
  mutate(tm=time+n-1)%>%
  dplyr::select(-n,-time)%>%
  rename(time="tm")
# create new event between previous night and next day morning event
rep_id=data.frame(id=rep(unique(ds$id),each=124),event=rep(c(1,2,3,4),31*length(unique(ds$id))),day=rep(seq(1,31,1),each=4)) %>%
  group_by(id)%>%
  mutate(time=1:n())%>%
  ungroup()%>%
  dplyr::select(id,event,time,day)


# deal with the person who finish 31 days measurements
ds3=dplyr::left_join(rep_id,ds,by=c("id","time","event"))%>%
  group_by(id)%>% 
  arrange(time,.by_group = TRUE)%>%
  filter(!((!id %in%l9) &time>120))

# change NAN to NA


ds$pa[is.nan(ds$pa)]<-NA
ds$na[is.nan(ds$na)]<-NA
ds=as.data.frame(ds)



# used amelia to imputate data
set.seed(123457)
# bds <- matrix(c(5, 0, 100,6, 0, 100,7, 0, 100,8, 0, 100), nrow = 4, ncol = 3,byrow = T)
ds_imp=amelia(x=ds,cs="id",ts="time",m = 10,p2s = 2,lags=c("pa","na"),leads=c("pa","na"),idvars = "day")
#bounds = bds)

summary(ds_imp)
#save(ds_imp, file = "ds_imp.RData")

#save(mem, file = "mem.RData")

#--------------clustering----------------------------------
# extracting coefficients of VAR models
# 
ds1=(ds_imp[[1]]$imp1+ds_imp[[1]]$imp2+ds_imp[[1]]$imp3+ds_imp[[1]]$imp4+ds_imp[[1]]$imp5+
       ds_imp[[1]]$imp6+ds_imp[[1]]$imp7+ds_imp[[1]]$imp8+ds_imp[[1]]$imp9+ds_imp[[1]]$imp10)/10
  # center by each individuals mean of pa and na
ds1_m=
  ds1%>%
  group_by(id)%>%
  mutate(m_pa=mean(pa),m_na=mean(na))%>%
  mutate(
    pa=pa-m_pa,
    na=na-m_na)%>%
  dplyr::select(id,time,pa,na)


#---------find optimal p for real data-----------------
nob=length(unique(ds1$id))

totp=c()
for(i in 1:nob){
  id=unique(ds1$id)[i]
  ser=ds1_m[ds1_m$id==id,]
  a<- VAR(ser[,-c(1:2)], type = "none", 
          season = NULL,lag.max = 3,ic="AIC")
  totp=c(totp,a$p)
}

table(totp)
#--------fit-----------------------
## LO for VAR(1) HO for VAR(3)
  # get the data for each individual
coef_f=list()
coef_f2=list()
for (m in 1:10){
  ds1=ds_imp[[1]][[m]]

nob=length(unique(ds1$id))
k=ncol(ds1_m)-2
ds2=list()
ds1_m=
  ds1%>%
  group_by(id)%>%
  mutate(m_pa=mean(pa),m_na=mean(na))%>%
  mutate(
    pa=pa-m_pa,
    na=na-m_na)%>%
  dplyr::select(id,time,pa,na)
for(i in 1:nob){
  id=unique(ds1$id)[i]
  ds2[[i]]=ds1_m[ds1_m$id==id,]
}

### get the VAR model's coefficients 
ser=ds2
model_per_pers = matrix(0, nob, (k)*k*1)
model_per_pers2 = matrix(0, nob, (k*3)*k)

for(i in 1:nob){
 # fit VAR(1) model
  v <- VAR(ser[[i]][,-c(1,2)], type = "none", 
           season = NULL,p=1)
  terms = length(v$varresult[[1]]$coefficients)
  coeff_pers = matrix(0, k, k*1)

  ###fit VAR(3) model
  v2 <- VAR(ser[[i]][,-c(1,2)], type = "none", 
            season = NULL,p=3)
  terms2 = length(v2$varresult[[1]]$coefficients)
  coeff_pers2 = matrix(0, k, k*3)
  for(vas in 1:k){
    # combined coefficients for each time serise variables
    coeff_pers[vas,1:(terms) ] = v$varresult[[vas]]$coefficients[1:(terms)]
    
    coeff_pers2[vas,1:(terms2) ] = v2$varresult[[vas]]$coefficients[1:(terms2)]
  }
  # combined cofficients for each person whithin in a sample
  model_per_pers[i, ] = as.vector(t(coeff_pers))
  
  model_per_pers2[i, ] = as.vector(t(coeff_pers2))
 
}
coef_f[[m]]=model_per_pers # combined coefficients for each imputation
coef_f2[[m]]=model_per_pers2
}


# calculate the mean of coefficent across imputations
model_per_pers=Reduce("+", coef_f)/m
model_per_pers2=Reduce("+", coef_f2)/m
#-----clustering---------------------
#-LO method--------------------
#--3 cluster---------------------
p=1
repet=100
mod.r=c()


# gmm
for (time in 1:repet){
  
  m=tryCatch(Mclust(model_per_pers, G=1:8,initialization = list(hcPairs =hcRandomPairs(model_per_pers))),
               error=function(e){return(NA)})
 
  mod.r=rbind(mod.r,m)
}

bic=rep(NA,1*repet)
for(run in 1:(1*repet)){
  bic[run]=mod.r[run,]$bic}
member=mod.r[which(bic==max(bic)),]$classification

table(member)

#member=m$classification

#--kmeans----------------

hc_dist=function(x) dist(x, method = 'euclidean')
hcl <- function(x, k) list(cluster=cutree(hclust(hc_dist(x), method = "average"),k=k))
a=clusGap(model_per_pers, FUN = hcl, K.max =8, B = 500)  
op.n=which(a$Tab[,3]==max(a$Tab[,3]))
h=hclust(dist(model_per_pers, method = 'euclidean'),method = 'average')
member_c=cutree(h, k =  op.n)
ini=c()
for(c in 1: op.n){
  mean=colMeans(model_per_pers[member_c==c,])
  ini=rbind(ini,mean)
}
m.k=kmeans(model_per_pers,centers=ini)

#m=kmeans(model_per_pers,op.n,nstart=500)


member.k=m.k$cluster

#----------HO method---------------------
#----------p=3--------------------------------
# gmm 3 cluster
p=3
repet=100
mod2=c()
for (time in 1:repet){
  m.r2=tryCatch(Mclust(model_per_pers2, G=1:8,initialization = list(hcPairs =hcRandomPairs(model_per_pers2))),
                       
                error=function(e){return(NA)})
  mod2=rbind(mod2,m.r2)
  
}

bic=rep(NA,1*repet)
for(run in 1:(1*repet)){
  bic[run]=mod2[run,]$bic}
member2=mod2[which(bic==max(bic)),]$classification

table(member2)

#---------kmeans----------------------------------
a=clusGap(model_per_pers2, FUN = hcl, K.max = 8, B = 500)  
op.n=which(a$Tab[,3]==max(a$Tab[,3]))
h=hclust(dist(model_per_pers2, method = 'manhattan'),method = 'average')
member2=cutree(h, k =op.n)
ini=c()
for(c in 1: op.n){
  mean=colMeans(model_per_pers2[member==c,])
  ini=rbind(ini,mean)
}
m.k2=kmeans(model_per_pers2,centers=ini)

member2.k=m.k2$cluster
#----summarized and evaluated the clustering results------------------
#-------coefficients by LO method with GMM --------------------------
mem=cbind(id=unique(ds1$id),member=member)%>%
  as.data.frame()
co_mean=cbind(mem,model_per_pers)%>%
  group_by(member)%>%
  summarize(across(where(is.numeric),~mean(.x,na.rm=T)))

cn=rep(c("PA(t-1)","NA(t-1)"),2)
cn=factor(cn,levels=c("PA(t-1)","NA(t-1)"))
rn=rep(c("PA(t)","NA(t)"),each=2)
a=data.frame(Value=t(co_mean[1,-c(1,2)]),cn=cn,rn=rn)%>%
  ggplot(aes(x=cn,y=rn))+
  geom_tile(aes(fill = Value))+ 
  scale_fill_gradient2(low = "black", high = "black")+ 
  labs(x = "",
       y = "",
       title="Cluster 1")+
  theme(axis.text.x=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12,face="bold"))+
  theme_bw()

b=data.frame(Value=t(co_mean[2,-c(1,2)]),cn=cn,rn=rn)%>%
  ggplot(aes(x=cn,y=rn))+
  geom_tile(aes(fill = Value))+ 
  scale_fill_gradient2(low = "black", high = "black")+ 
  labs(x = "",
       y = "",
       title="Cluster 2")+
  theme(axis.text.x=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12,face="bold"))+
  theme_bw()

c=data.frame(Value=t(co_mean[3,-c(1,2)]),cn=cn,rn=rn)%>%
  ggplot(aes(x=cn,y=rn))+
  geom_tile(aes(fill = Value))+ 
  scale_fill_gradient2(low = "black", high = "black")+ 
  labs(x = "",
       y = "",
       title="Cluster 3")+
  theme(axis.text.x=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12,face="bold"))+
  theme_bw()

grid.arrange(arrangeGrob(a,b,c,ncol=1),bottom=textGrob("(a)"))



#------coefficients by HO method with GMM---------------------
mem2=cbind(id=unique(ds1$id),member=member2)%>%
  as.data.frame()
co_mean2=cbind(mem2,model_per_pers2)%>%
  group_by(member)%>%
  summarize(across(where(is.numeric),~mean(.x,na.rm=T)))

matrix(co_mean2[1,-c(1,2)],nrow=2,ncol=6,byrow=T)
matrix(co_mean2[2,-c(1,2)],nrow=2,ncol=6,byrow=T)
matrix(co_mean2[3,-c(1,2)],nrow=2,ncol=6,byrow=T)

cn=rep(c("PA(t-1)","NA(t-1)","PA(t-2)","NA(t-2)","PA(t-3)","NA(t-3)"),2)
cn=factor(cn,levels=c("PA(t-1)","NA(t-1)","PA(t-2)","NA(t-2)","PA(t-3)","NA(t-3)"))
rn=rep(c("PA(t)","NA(t)"),each=6)
d=data.frame(Value=t(co_mean2[1,-c(1,2)]),cn=factor(cn),rn=rn)%>%
  ggplot(aes(x=cn,y=rn))+
  geom_tile(aes(fill = Value))+ 
  scale_fill_gradient2(low = "black", high = "black")+ 
  labs(x = "",
       y = "",
       title="Cluster 1")+
  theme(axis.text.x=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12,face="bold"))+
  theme_bw()

e=data.frame(Value=t(co_mean2[2,-c(1,2)]),cn=factor(cn),rn=rn)%>%
  ggplot(aes(x=cn,y=rn))+
  geom_tile(aes(fill = Value))+ 
  scale_fill_gradient2(low = "black", high = "black")+ 
  labs(x = "",
       y = "",
       title="Cluster 2")+
  theme(axis.text.x=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12,face="bold"))+
  theme_bw()

f=data.frame(Value=t(co_mean2[3,-c(1,2)]),cn=factor(cn),rn=rn)%>%
  ggplot(aes(x=cn,y=rn))+
  geom_tile(aes(fill = Value))+ 
  scale_fill_gradient2(low = "black", high = "black")+ 
  labs(x = "",
       y = "",
       title="Cluster 3")+
  theme(axis.text.x=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12,face="bold"))+
  theme_bw()

grid.arrange(arrangeGrob(a,b,c,ncol=1))



grid.arrange(arrangeGrob(d,e,f,ncol=1),bottom=textGrob("(b)"))




#------------evaluation for LO with GMM----------------------
library(ggplot2)
mem=cbind(id=unique(ds1$id),member=member)%>%
  as.data.frame()
dt_ex <- load(dt_ex)

# using big-five traits, happyness depression anxiety stress
pl=mem%>%
  left_join(dt_ex,by="id")%>%
  dplyr::select(id,member,neo_extraversion,neo_openness,neo_agreeable,neo_conscient,neo_neurot,happyindex,dassdep,dassanx,dassstress,age,gender)%>%
  group_by(member)%>%
  summarize(
    neo_extraversion_s=sd(neo_extraversion,na.rm=T),
    neo_openness_s=sd(neo_openness,na.rm=T),
    neo_agreeable_s=sd(neo_agreeable,na.rm=T),
    neo_conscient_s=sd(neo_conscient,na.rm=T),
    neo_neurot_s=sd(neo_neurot,na.rm=T),
    neo_extraversion=mean(neo_extraversion,na.rm=T),
    neo_openness=mean(neo_openness,na.rm=T),
    neo_agreeable=mean(neo_agreeable,na.rm=T),
    neo_conscient=mean(neo_conscient,na.rm=T),
    neo_neurot=mean(neo_neurot,na.rm=T),
    
    happyindex=sd(happyindex,na.rm=T),
    dassdep=sd(dassdep,na.rm=T),
    dassanx=sd(dassanx,na.rm=T),
    dassstress=sd(dassstress,na.rm=T),
    agep=sd(age,na.rm=T),
    genderp=sd(gender,na.rm=T),
    
    
    # happyindex=mean(happyindex,na.rm=T),
    # dassdep=mean(dassdep,na.rm=T),
    # dassanx=mean(dassanx,na.rm=T),
    # dassstress=mean(dassstress,na.rm=T),
    # agep=mean(age,na.rm=T),
    # genderp=mean(gender,na.rm=T)
  )



pp1=pl[,c(1,7,2)]
pp2=pl[,c(1,8,3)]
pp3=pl[,c(1,9,4)]
pp4=pl[,c(1,10,5)]
pp5=pl[,c(1,11,6)]
names(pp1)=names(pp2)=names(pp3)=names(pp4)=names(pp5)=c("member","neo","neo_sd")
pl1=rbind(pp1,pp2,pp3,pp4,pp5)
pl1=cbind(pl1,type=rep(c("Extraversion","Openness","Agreeableness","Conscientiousness","Neuroticism"),each=length(unique(mem$member))))
pl1$member=as.factor(pl1$member)
names(pl1)[1]="Cluster"

# plot for mean of big-five trait for three clusters
pl1%>%
  ggplot(aes(x=type,y=neo,group=Cluster,color=Cluster))+
  geom_point()+
  geom_line(aes(linetype=Cluster,color = Cluster), size=1)+
  # geom_errorbar(aes(ymin=neo-1.96*neo_sd, ymax=neo+1.96*neo_sd,color=Cluster),
  #               width =0.005,linetype = 1,size=1)
  labs(
    y = "Score",
    x="Personality")+
  
  geom_point(size=4,aes(shape=Cluster,color=Cluster))+
  scale_shape_manual(values=c(15, 16, 17))+
  scale_color_manual(values=c("black","#5A5A5A","lightgrey"))+
  theme(
    panel.background = element_rect(fill = "white"),
    legend.key  = element_rect(fill = "white"),
    axis.line.x = element_line(colour = "black", size = 1),
    axis.line.y = element_line(colour = "black", size = 1),
    axis.text.x=element_text(size=12,face="bold"),
    axis.text.y=element_text(size=12,face="bold"),
    axis.title.y = element_text(size = 12,face="bold"),
    axis.title.x = element_text(size = 12,face="bold")
  )  



pp11=pl[,c(1,7)]
pp12=pl[,c(1,8)]
pp13=pl[,c(1,9)]
pp14=pl[,c(1,10)]
# plot for mean of happyness depression anxiety and stress for three clusters
names(pp11)=names(pp12)=names(pp13)=names(pp14)=c("member","emo")
pl2=rbind(pp11,pp12,pp13,pp14)
pl2=cbind(pl2,type=rep(c("happy","depression","anxiety","stress"),each=(length(unique(mem$member)))))
pl2$member=as.factor(pl2$member)
pl2%>%
  ggplot(aes(x=type,y=emo,group=member,color=member))+
  geom_point()+
  geom_line()

###------calculate the the gender age for each cluster---------------------

mem%>%
  left_join(data,by="id")%>%
  group_by(member)%>%
  summarize(
            age=mean(age_diary,na.rm=T),
            agesd=sd(age_diary,na.rm=T),
            gender=mean(sex,na.rm=T))


#---------mean and sd for trait level of NA and PA and correlation between trait level NA and PA-----------------


data_or=mem%>%
  left_join(data,by="id")%>%
  rowwise()%>%
  mutate(
    pad=(mad_diary_5+mad_diary_11+mad_diary_13)/3,
    nad=(mad_diary_6+mad_diary_14+mad_diary_16)/3,
    paa=(mad_diary_7+mad_diary_9+mad_diary_15)/3,
    naa=(mad_diary_8+mad_diary_10+mad_diary_12)/3)%>%
  mutate(pa=(pad+paa)/2,
         na=(nad+naa)/2)%>%
  dplyr::select(id,member,pa,na)

data_or%>%
  dplyr::group_by(member)%>%
summarize(pa_sd=sd(pa,na.rm = T),pa=mean(pa,na.rm = T),na_sd=sd(na,na.rm = T),na=mean(na,na.rm = T))

ds_anova=ds%>%
  dplyr::group_by(id)%>%
  
  dplyr::summarize(pa=mean(pa,na.rm=T),na=mean(na,na.rm=T))%>%
  dplyr::left_join(mem,by="id")


cor.test(ds_anova$pa,ds_anova$na)
te=cbind(ds_anova[ds_anova$member==1,"pa"],ds_anova[ds_anova$member==1,"na"])
cor.test(te$pa,te$na)

te=cbind(ds_anova[ds_anova$member==2,"pa"],ds_anova[ds_anova$member==2,"na"])
cor.test(te$pa,te$na)

te=cbind(ds_anova[ds_anova$member==3,"pa"],ds_anova[ds_anova$member==3,"na"])
cor.test(te$pa,te$na)



cor.test(ds_anova$pa,ds_anova$na)

ds_anova$member=as.factor(ds_anova$member)
anova(lm(pa~member,data=ds_anova))
anova(lm(na~member,data=ds_anova))


TukeyHSD(aov(pa~member,data=ds_anova))
TukeyHSD(aov(na~member,data=ds_anova))




