library(mvtnorm)
library(truncnorm)
library(ROCR)
library(polspline)
library(spatstat)
library(MCMCpack)
library(MASS)
library(parallel)
library(Matrix)
library(sparseMVN)
library(fields)
library(RcppTN)
library(devtools)
library(Rcpp)
library(RcppArmadillo)
library(inline)
library(tmvtnorm)
n.cores=44
nearest_neighbor=10
n.knots=10
n.sample=40000
index=1:44
b=0.3
rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

source("Base_model_functions.R")
PRED=function(data.pred, data.train, alpha){
  k=1
  dec.mat=array(list(NULL), c(npat, 3))
  logpden=apply(llist, 1, FUN=function(x){LogPredDen(data.train, data.pred[1,,drop=F], x, delta,omega=omega.hat)})
  pden=exp(logpden-mean(logpden))
  
  ##### P
  PHAT=data.pred[[1,3]]
  priprob=apply(llist, 1, FUN=function(x){kk=sum(x==1); PHAT^kk*(1-PHAT)^(n.new-kk)})
  pden=priprob*pden
  
  risk=sapply(1:nrow(llist) , FUN=function(loop){
    d=llist[loop,]
    n.fp = apply(llist, 1,FUN=function(y){sum(d==1 & y==0)})
    n.fn = apply(llist, 1,FUN=function(y){sum(d==0 & y==1)})
    risk=sum((pden)*((1-alpha)*n.fn+alpha*n.fp))
  })
  
  dec.mr=as.numeric(llist[which.min(risk),])
  dec.mat[[1,1]]=data.pred[[1,4]]
  dec.mat[[1,2]]=dec.mr
  dec.mat[[1,3]]=pden
  return(dec.mat)
}
source("functions.R")
sourceCpp('functions.cpp')

alpha=0.5
m=4 
delta=m-1 
p=4
npat=1
n.new=1
llist=matrix(0:1,nrow=2)

#### Using the motivating data set to estimate region-specific MRI distributions:
load("msi_fillrecenter.RData") # motivating data set
fillin=fillrecenter
rm(fillrecenter)
index=1:46
for (i in 1:46)
{
  fillin[[i]][,"ADC"]=fillin[[i]][,"ADC"]/100
  fillin[[i]][,c("KTRANS","KEP","AUGC")]=log(fillin[[i]][,c("KTRANS","KEP","AUGC")])
}
########## Generate data for jags: 
data=cbind(as.data.frame(fillin[[1]][,c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone")]),1)
colnames(data)=c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone","subject")
for (k in c(index[-1]))
{
  temp=cbind(as.data.frame(fillin[[k]][,c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone")]),k)
  colnames(temp)=c("ADC","KTRANS","KEP","AUGC","cancer","case","x","y","zone","subject")
  data=rbind(data,temp)
}
dataNC=data[data[,"cancer"]==0,]
dataC=data[data[,"cancer"]==1,]
DATA=data
rm(data)
####### average pcancer per region:
ppz=nrow(DATA[(DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ"),])/nrow(DATA[(DATA[,"zone"]=="PZ"),])
pcg=nrow(DATA[(DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG"),])/nrow(DATA[(DATA[,"zone"]=="CG"),])

meancpz=colMeans(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
meanccg=colMeans(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
meanncpz=colMeans(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
meannccg=colMeans(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
meannc=colMeans(DATA[DATA[,"cancer"]==0,c("ADC","KTRANS","KEP","AUGC")])
meanc=colMeans(DATA[DATA[,"cancer"]==1,c("ADC","KTRANS","KEP","AUGC")])
#means=list(meannc,meanc)
means=list(meancpz,meanccg,meanncpz,meannccg)
covcpz=var(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
covccg=var(DATA[((DATA[,"cancer"]==1)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
covncpz=var(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="PZ")),c("ADC","KTRANS","KEP","AUGC")])
covnccg=var(DATA[((DATA[,"cancer"]==0)&(DATA[,"zone"]=="CG")),c("ADC","KTRANS","KEP","AUGC")])
covnc=var(DATA[(DATA[,"cancer"]==0),c("ADC","KTRANS","KEP","AUGC")])
covc=var(DATA[(DATA[,"cancer"]==1),c("ADC","KTRANS","KEP","AUGC")])
#covs=list(covnc,covc)
covs=list(covcpz,covccg,covncpz,covnccg)
#save(means,covs,file="~/Google Drive/mpMRI/17.1210/means_covs_region.RData")
Sigmas=list()
Sigmas[[1]]=diag(c(5,0.2,0.25,0.2))
Sigmas[[2]]=diag(c(10,0.4,0.5,0.4))
omega.delta=diag(rep(1,4))
################# 
tausq=0.2
the1=c(1,5,20) #sigmasq
the2=c(0.5,2,10) # phi: inverse of range parameters
the3=c(0.5,1.5) #nu

setting=expand.grid(x=the1,y=the2,z=the3)
### Sigmasq: 
Sigma=Sigmas[[2]]
tausq=1
pcg=0.2
ppz=pcg+0.2
index=1:44
rm(fillin)

######## code for one setting:
set=1

auc_mregion=tp_mregion=auc_msmooth=tp_msmooth=auc_sse=tp_sse=auc_ssenngp=tp_ssenngp=auc_sselr=tp_sselr=auc_ssecar=tp_ssecar=auc_full=tp_full=numeric()
aucmregion=tpmregion=aucmsmooth=tpmsmooth=aucsse=tpsse=aucssenngp=tpssenngp=aucsselr=tpsselr=aucssecar=tpssecar=aucfull=tpfull=numeric()

for (S in 1:100)
{
sigmasq=setting[set,1]
phi=setting[set,2]
nu=setting[set,3]
load("fillin_reduced.RData") # simulate from reduced-size images
fillnew=list()
new=1
select=sample(1:46,size=44,replace=TRUE)
for (i in select)
{
  temp=fillin[[i]]
  zone=1*(temp[,"zone"]=="PZ")+0*(temp[,"zone"]=="CG")
  coords=temp[,c("x","y")]
  D <- as.matrix(dist(coords))
  #R <- exp(-phi*D)
  R <- materncovcpp(D,sigmasq,phi,nu)
  qs0 <- rmvn(1, rep(0,nrow(temp)), R)+sapply(1:nrow(temp),FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
  qs <- qs0 + rnorm(n=nrow(temp),mean=0,sd=sqrt(tausq))
  while (sum(qs>0)==0) {
    qs0 <- rmvn(1, rep(0,nrow(temp)), R)+sapply(1:nrow(temp),FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
    qs <- qs0 + rnorm(n=nrow(temp),mean=0,sd=sqrt(tausq))
  }
  ps <- sapply(1:nrow(temp),FUN=function(x){pnorm(q=qs[x],mean=0,sd=1)})
  cancer=1*(qs>0)+0*(qs<0)
  #print(i)
  temp[,"cancer"]=cancer
  temp=cbind(temp,qs0,qs,ps)
  B=mvrnorm(n=1,mu=rep(0,4),Sigma=Sigma)
  for (j in 1:nrow(temp))
  {
    indicator=1*((cancer[j]==1)&(zone[j]==1))+2*((cancer[j]==1)&(zone[j]==0))+3*((cancer[j]==0)&(zone[j]==1))+4*((cancer[j]==0)&(zone[j]==0))
    A=mvrnorm(n=1,mu=means[[indicator]],Sigma=covs[[indicator]])
    temp[j,c("ADC","KTRANS","KEP","AUGC")]=A+B
  }
  fillnew[[new]]=temp
  new=new+1
}

fillrecenter=fillin=fillnew
for(i in index)
{
  temp=fillin[[i]]
  temp[,"area"]=1*(temp[,"cancer"]==1)+2*((temp[,"cancer"]==0)&(temp[,"zone"]=="PZ"))+3*((temp[,"cancer"]==0)&(temp[,"zone"]=="CG"))
  temp[,"parea"]=1*(temp[,"cancer"]==1)+2*((temp[,"cancer"]==0)&(temp[,"Z"]==2))+3*((temp[,"cancer"]==0)&(temp[,"Z"]==1))
  temp[,"region"]=1*((temp[,"cancer"]==1)&(temp[,"zone"]=="PZ"))+2*((temp[,"cancer"]==1)&(temp[,"zone"]=="CG"))+3*((temp[,"cancer"]==0)&(temp[,"zone"]=="PZ"))+4*((temp[,"cancer"]==0)&(temp[,"zone"]=="CG"))
  temp[,"pregion"]=1*((temp[,"cancer"]==1)&(temp[,"Z"]==2))+2*((temp[,"cancer"]==1)&(temp[,"Z"]==1))+3*((temp[,"cancer"]==0)&(temp[,"Z"]==2))+4*((temp[,"cancer"]==0)&(temp[,"Z"]==1))
  fillrecenter[[i]]=temp
  #print(i)
}
fillin=fillnew=fillrecenter
data=fillnew[[1]]
data[,"subject"]=1
for (k in 2:44)
{
  temp=fillnew[[k]]
  temp[,"subject"]=k
  data=rbind(data,temp)
}
dataNC=DATANC=data[data[,"cancer"]==0,]
dataC=DATAC=data[data[,"cancer"]==1,]
DATA=data
rm(data)

rdata=list()
ldata=list()
ydata=list()
adata=list()
region=list()
zone=list()
TZ=list()
Z=list()
ZP=list()

for (K in 1:44)
{
  rdata[[K]]=cbind(fillrecenter[[K]][,c("ADC","KTRANS","KEP","AUGC")])
  #rdata[[K]][,c("KTRANS","KEP","AUGC")]=log(rdata[[K]][,c("KTRANS","KEP","AUGC")])
  ldata[[K]]=as.matrix(cbind(fillrecenter[[K]][,c("|x|","y","zone")]))
  colnames(ldata[[K]])=c("loc.x","loc.y","loc")
  ydata[[K]]=fillrecenter[[K]][,c("cancer")]
  adata[[K]]=fillrecenter[[K]][,c("area")]
  region[[K]]=fillrecenter[[K]][,c("region")]
  zone[[K]]=fillrecenter[[K]][,c("zone")]
  TZ[[K]]=fillrecenter[[K]][,c("TZ")]
  Z[[K]]=fillrecenter[[K]][,c("Z")]
  ZP[[K]]=fillrecenter[[K]][,c("ZP")]
}
####################################################
######## Step 0: Smoothing Bandwidth Tuning ########
####################################################
whole=DATA[,-1]
wholePZ=DATA[DATA$zone=="PZ",-1]
wholeCG=DATA[DATA$zone=="CG",-1]

auc_mregion=tp_mregion=numeric()
auc_msmooth=tp_msmooth=numeric()
##############################################
######## Step 1: Mregion and Msmooth: ########
##############################################
trainC=DATAC[(DATAC[,"subject"]%in%c(1:34)),]
trainNC=DATANC[(DATANC[,"subject"]%in%c(1:34)),]
trainPZ=DATA[((DATA[,"subject"]%in%c(1:34))&(DATA$zone=="PZ")),]
trainCG=DATA[((DATA[,"subject"]%in%c(1:34))&(DATA$zone=="CG")),]
train=DATA[(DATA[,"subject"]%in%c(1:34)),]

test=DATA[(DATA[,"subject"]%in%c(35:44)),]
testNC=DATANC[(DATANC[,"subject"]%in%c(35:44)),]

### PZ:2, CG:1
train[,"ZONE"]=ifelse(train[,"zone"]=="PZ",2,1)
fitp=polymars(train[,"cancer"], as.factor(train[,"ZONE"]), classify = T,maxsize=4)
W=list()
for (qq in 1:44)
{
  temp=whole[whole[,"subject"]==qq,]
  p0=predict(fitp,as.factor(temp[,"Z"]),classify=F)[,1]
  OW=ifelse(temp[,"Z"]==2,1,2)
  p1=predict(fitp,as.factor(OW),classify=F)[,1]
  pPZ=ifelse(temp[,"Z"]==2,p0,p1)
  pCG=ifelse(temp[,"Z"]==1,p0,p1)
  W[[qq]]=pPZ*ZP[[qq]]+pCG*(1-ZP[[qq]])
}
############## Training Data:
mydata.train=data.frame()
mydata.train[1:34,1]=1:34
mydata.train$par=rdata[1:34]
mydata.train$paxis=W[1:34] # predicted p
mydata.train$cancer=ydata[1:34]
mydata.train$area=adata[1:34]
mydata.train$region=region[1:34]
mydata.train$predzone=Z[1:34]  #### dat[[i,7]]: zone for test voxel! in logpredden
mydata.train$predzonep=ZP[1:34]

phat=sum(unlist(mydata.train[,4])==1)/length(unlist(mydata.train[,4]))
dec.mat=array(list(NULL), c(npat, 3))

ParmFit=DataFit(mydata.train)
Nreg=omega.hat=NULL
for(z in 1:4)
{omega.hat[z]=list(ParmFit$sigma.hat[[z]]*delta/(length(ParmFit$r[[z]])-1))}

bpre=btru=list()
auc=tp=auc_s=tp_s=numeric()
ss=1
for (j in 35:44)
{
  bpre[[ss]]=numeric()
  btru[[ss]]=numeric()
  simu=function(P)
  {
    mydata.test=data.frame()
    mydata.test[1,1]=j
    mydata.test$par=list(rdata[[j]][P,])
    mydata.test$paxis=list(W[[j]][P])
    mydata.test$cancer=list(ydata[[j]][P])
    mydata.test$area=list(adata[[j]][P])
    mydata.test$region=list(region[[j]][P])
    mydata.test$predzone=list(Z[[j]][P])
    mydata.test$predzonep=list(ZP[[j]][P])
    a=PRED(data.pred=mydata.test,data.train=mydata.train, alpha)
    return(a)
  }
  results=mclapply(1:nrow(rdata[[j]]),simu,mc.cores = n.cores)
  for (i in 1:length(results))
  {
    bpre[[ss]][i]=as.numeric(results[[i]][1,3][[1]][2]/sum(results[[i]][1,3][[1]]))
    btru[[ss]][i]=as.numeric(results[[i]][1,1][[1]])
  }
  bpred <- prediction(bpre[[ss]], btru[[ss]])
  bperf<- performance(bpred,"tpr","fpr")
  bau <- performance(bpred,"auc")
  auc[ss] <- unlist(slot(bau, "y.values"))
  tp[ss]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  
  Xcoord=fillin[[j]][,"x"];Ycoord=fillin[[j]][,"y"];
  XX <- ppp(Xcoord, Ycoord, c(-1,1),c(-1,1),marks=bpre[[ss]])
  #b <- bw.smoothppp(X)  # least-squares CV to select a smoothing bandwidth for spatial smoothing of marks.
  Spre=as.numeric(Smooth(XX, sigma=b,at="points",edge=TRUE, diggle=FALSE))
  #### Msmooth:
  bpred <- prediction(Spre, btru[[ss]])
  bperf<- performance(bpred,"tpr","fpr")
  bau <- performance(bpred,"auc")
  auc_s[ss] <- unlist(slot(bau, "y.values"))
  tp_s[ss]=max(bperf@y.values[[1]][bperf@x.values[[1]]<=0.1])
  ss=ss+1
}
aucmregion[S]=mean(auc)
tpmregion[S]=mean(tp)
#print(auc_mregion)
aucmsmooth[S]=mean(auc_s)
tpmsmooth[S]=mean(tp_s)
#print(auc_msmooth)

########################################################
###################### SSE + SP : ######################
########################################################
whole=DATA
####### Estimating delta's and Sigma:
ind=1:34
trainC=DATAC[(DATAC[,"subject"]%in%c(1:34)),]
trainNC=DATANC[(DATANC[,"subject"]%in%c(1:34)),]
trainPZ=DATA[((DATA[,"subject"]%in%c(1:34))&(DATA$zone=="PZ")),]
trainCG=DATA[((DATA[,"subject"]%in%c(1:34))&(DATA$zone=="CG")),]
train=DATA[(DATA[,"subject"]%in%c(1:34)),]

### PZ:2, CG:1
train[,"ZONE"]=ifelse(train[,"zone"]=="PZ",2,1)
fitp=polymars(train[,"cancer"], as.factor(train[,"ZONE"]), classify = T,maxsize=4)
W=list()
for (qq in index)
{
  temp=whole[whole[,"subject"]==qq,]
  p0=predict(fitp,as.factor(temp[,"TZ"]),classify=F)[,1]
  W[[qq]]=p0
}

############## Training Data:
mydata.train=data.frame()
mydata.train[1:34,1]=1:34
mydata.train$par=rdata[1:34]
mydata.train$paxis=W[1:34] # predicted p
mydata.train$cancer=ydata[1:34]
mydata.train$area=adata[1:34]
mydata.train$region=region[1:34]
mydata.train$predzone=Z[1:34]  #### dat[[i,7]]: zone for test voxel! in logpredden
mydata.train$predzonep=ZP[1:34]

phat=sum(unlist(mydata.train[,4])==1)/length(unlist(mydata.train[,4]))
dec.mat=array(list(NULL), c(npat, 3))

ParmFit=DataFit(mydata.train)
Nreg=omega.hat=NULL
for(z in 1:4)
{omega.hat[z]=list(ParmFit$sigma.hat[[z]]*delta/(length(ParmFit$r[[z]])-1))}

ni1=nrow(trainC[trainC[,"zone"]=="PZ",])
ni2=nrow(trainC[trainC[,"zone"]=="CG",])
ni3=nrow(trainNC[trainNC[,"zone"]=="PZ",])
ni4=nrow(trainNC[trainNC[,"zone"]=="CG",])

###### Initial values:
mus=gammas=list()
s=1
mus=means
gammas=covs
Sigmasq=diag(c(10,0.4,0.5,0.4))

y=lapply(1:44,FUN=function(x){as.matrix(fillnew[[x]][,c("ADC","KTRANS","KEP","AUGC")])})
#### for LR:
knots=t(sapply(1:44,FUN=function(x){sample(1:nrow(fillnew[[x]]),size=n.knots,replace=FALSE)}))
ci=Cistar=covmat=invCistar=list()
for (x in 1:44)
{
  coordsi=fillnew[[x]][,c("x","y")]
  Di=as.matrix(dist(coordsi))
  Diknots=Di[knots[x,],knots[x,]]
  Cistar[[x]]=materncovcpp(Diknots,sigmasq,phi,nu) #
  invCistar[[x]]=solvecpp(Cistar[[x]])
  Direst=Di[knots[x,],]
  ci[[x]]= materncovcpp(Direst,sigmasq,phi,nu) #
  covmat[[x]]=t(ci[[x]])%*%solvecpp(Cistar[[x]])%*%ci[[x]]+diag(tausq,nrow(fillnew[[x]]))
}
indicators=1*((train[,"cancer"]==1)&(train[,"zone"]=="PZ"))+2*((train[,"cancer"]==1)&(train[,"zone"]=="CG"))+3*((train[,"cancer"]==0)&(train[,"zone"]=="PZ"))+4*((train[,"cancer"]==0)&(train[,"zone"]=="CG"))

ss=1
for (slice in 35:44)
{
  u=kk=k=j=slice
  test=DATA[DATA$subject==k,]
  testNC=DATANC[DATANC[,"subject"]==k,]
  load("nngpdata_reduced_201-46.RData")
  nIndx=nngpdata[[select[slice]]]$nIndx
  nnIndx=nngpdata[[select[slice]]]$nnIndx
  nnIndxLU=nngpdata[[select[slice]]]$nnIndxLU
  uIndx=nngpdata[[select[slice]]]$uIndx
  uIndxLU=nngpdata[[select[slice]]]$uIndxLU
  uiIndx=nngpdata[[select[slice]]]$uiIndx
  CIndx=nngpdata[[select[slice]]]$CIndx
  c.update=nngpdata[[select[slice]]]$c
  d=nngpdata[[select[slice]]]$d
  D=nngpdata[[select[slice]]]$D
  n=nngpdata[[select[slice]]]$n
  rm(nngpdata)
  
  bcf_updates=updateBF_materncpp(nearest_neighbor,d, D, nIndx, nnIndxLU, CIndx, n, sigmasq, tausq, phi, nu)
  c.update=bcf_updates$c
  B.update=bcf_updates$B
  C.update=bcf_updates$C
  FF.update=bcf_updates$FF
  
  teind=slice
  ind=1:34
  
  test=DATA[DATA[,"subject"]==slice,]
  ####### Priors:
  
  cancertest_sse=qupdate_sse=matrix(0,n.sample,nrow(test))
  cancertest_ssenngp=wupdate_ssenngp=kupdate_ssenngp=matrix(0,n.sample,nrow(test))
  cancertest_sselr=wupdate_sselr=matrix(0,n.sample,nrow(test))
  cancertest_ssecar=qupdate_ssecar=matrix(0,n.sample,nrow(test))
  cancertest_full=qupdate_full=matrix(0,n.sample,nrow(test))
  
  ###### Random Initial for LR: qupdate follows the nngp, qupdate +nu is the probit of p(cancer)
  wupdate_sselr[1,]=sapply(1:n,FUN=function(x){rnorm(1,0,sqrt(sigmasq+tausq))})
  cancertest_sselr[1,]=sapply(1:n,FUN=function(x){rbinom(n=1,size=1,prob=pnorm(0,0,1))})
  
  ###### if we dont update theta, then this part is constant:
  mecoef=matrix(NA,n,n-1)
  va=numeric()
  for (ii in 1:n)
  {
    cijtilde=covmat[[slice]][ii,-ii]
    ci_j=ci[[slice]][,-ii]
    mecoef[ii,]=cijtilde%*%(diag(1/tausq,(n-1))-(1/(tausq)^2)*t(ci_j)%*%solvecpp(Cistar[[slice]]+ci_j%*%t(ci_j)/tausq)%*%ci_j)
    va[ii]=cijtilde%*%(diag(1/tausq,(n-1))-(1/(tausq)^2)*t(ci_j)%*%solvecpp(Cistar[[slice]]+ci_j%*%t(ci_j)/tausq)%*%ci_j)%*%t(t(cijtilde))
  }
  #### CAR:
  test=DATA[DATA$subject==slice,]
  testNC=DATANC[DATANC[,"subject"]==slice,]
  temp=fillnew[[slice]]
  zone=1*(temp[,"zone"]=="PZ")+0*(temp[,"zone"]=="CG")
  coords=temp[,c("x","y")]
  Dinv <- 1/as.matrix(dist(coords))
  for (x in 1:nrow(Dinv))
  {
    Dinv[x,x]=0
  }
  Dinv <- Dinv/rowSums(Dinv)
  cancertest_ssecar=qupdate_ssecar=matrix(0,n.sample,nrow(test))
  qupdate_ssecar[1,]=sapply(1:n,FUN=function(x){rnorm(1,0,sqrt(sigmasq+tausq))})
  cancertest_ssecar[1,]=sapply(1:n,FUN=function(x){rbinom(n=1,size=1,prob=pnorm(0,0,1))})
  
  ###### Full model:
  coords=test[,c("x","y")]
  DD <- as.matrix(dist(coords))
  #RR <- exp(-phi*D)
  Cifull <- materncovcpp(DD,sigmasq,phi,nu) + diag(tausq,n)
  ###### if we dont update theta, then this part is constant:
  meanfull=matrix(NA,n,n-1)
  varfull=numeric()
  for (ii in 1:n)
  {
    Cj_j=Cifull[ii,-ii]
    C_j_j=Cifull[-ii,-ii]
    invtemp=solvecpp(C_j_j)
    meanfull[ii,]=Cj_j%*%invtemp
    varfull[ii]=Cj_j%*%invtemp%*%t(t(Cj_j))
  }
  ###### Random Initial: qupdate follows the nngp, qupdate +nu is the probit of p(cancer)
  ###### 
  Delta=lapply(1:34,FUN=function(x){rep(0,4)})
  zone=1*(fillnew[[slice]][,"zone"]=="PZ")+0*(fillnew[[slice]][,"zone"]=="CG")
  lambda=sapply(1:n,FUN=function(x){qnorm(p=zone[x]*ppz+(1-zone[x])*pcg,mean=0,sd=1)})
  pcancers=zone*ppz+(1-zone)*pcg
  
  ##### Update q:
  cancertest_ssenngp[1,]=sapply(1:n,FUN=function(x){rbinom(n=1,size=1,prob=pnorm(0,0,1))})
  wupdate_ssenngp[1,]=rep(0,n)
  kijresults = nngpkijupdatecpp(n, lambda, rep(0,n), cancertest_ssenngp[1,], tausq)
  kupdate_ssenngp=kijresults$kij
  
  wupdate_sselr[1,]=sapply(1:n,FUN=function(x){rnorm(1,0,sqrt(sigmasq+tausq))})
  cancertest_sselr[1,]=sapply(1:n,FUN=function(x){rbinom(n=1,size=1,prob=pnorm(0,0,1))})
  qupdate_ssecar[1,]=sapply(1:n,FUN=function(x){rnorm(1,0,sqrt(sigmasq+tausq))})
  cancertest_ssecar[1,]=sapply(1:n,FUN=function(x){rbinom(n=1,size=1,prob=pnorm(0,0,1))})
  qupdate_full[1,]=sapply(1:n,FUN=function(x){rnorm(1,0,sqrt(sigmasq+tausq))})
  cancertest_full[1,]=sapply(1:n,FUN=function(x){rbinom(n=1,size=1,prob=pnorm(0,0,1))})
  
  auc_sse_temp=auc_ssenngp_temp=auc_sselr_temp=auc_ssecar_temp=auc_full_temp=numeric()
  ####### MCMC Iteration:
  for (s in 2:n.sample)
  {
    ##### update Gamma's:
    zo=1*(fillnew[[ind[1]]][,"zone"]=="PZ")+0*(fillnew[[ind[1]]][,"zone"]=="CG")
    can=fillnew[[ind[1]]][,"cancer"]
    indicator=1*((can==1)&(zo==1))+2*((can==1)&(zo==0))+3*((can==0)&(zo==1))+4*((can==0)&(zo==0))
    Y=y[[ind[1]]]-t(sapply(1:nrow(fillnew[[ind[1]]]),FUN=function(x){mus[[indicator[x]]]}))-matrix(rep(Delta[[ind[1]]],nrow(fillnew[[ind[1]]])),ncol=4,byrow=T)
    for (x in ind[-1])
    {
      zo=1*(fillnew[[x]][,"zone"]=="PZ")+0*(fillnew[[x]][,"zone"]=="CG")
      can=fillnew[[x]][,"cancer"]
      indicator=1*((can==1)&(zo==1))+2*((can==1)&(zo==0))+3*((can==0)&(zo==1))+4*((can==0)&(zo==0))
      tem=y[[x]]-t(sapply(1:nrow(fillnew[[x]]),FUN=function(xx){mus[[indicator[xx]]]}))-matrix(rep(Delta[[x]],nrow(fillnew[[x]])),ncol=4,byrow=T)
      Y=rbind(Y,tem)
    }
    S1=t(Y[indicators==1,])%*%Y[indicators==1,]
    S2=t(Y[indicators==2,])%*%Y[indicators==2,]
    S3=t(Y[indicators==3,])%*%Y[indicators==3,]
    S4=t(Y[indicators==4,])%*%Y[indicators==4,]
    gammas[[1]]=riwish(v=ni1+delta,S=S1+omega.hat[[1]])
    gammas[[2]]=riwish(v=ni2+delta,S=S2+omega.hat[[2]])
    gammas[[3]]=riwish(v=ni3+delta,S=S3+omega.hat[[3]])
    gammas[[4]]=riwish(v=ni4+delta,S=S4+omega.hat[[4]])
    
    ##### Update mu's:
    Y=y[[ind[1]]]-matrix(rep(Delta[[ind[1]]],nrow(fillnew[[ind[1]]])),ncol=4,byrow=T)
    for (x in ind[-1])
    {
      tem=y[[x]]-matrix(rep(Delta[[x]],nrow(fillnew[[x]])),ncol=4,byrow=T)
      Y=rbind(Y,tem)
    }
    mus[[1]]=mvrnorm(n=1,mu=colSums(Y[indicators==1,]/ni1),Sigma=gammas[[1]]/ni1)
    mus[[2]]=mvrnorm(n=1,mu=colSums(Y[indicators==2,]/ni2),Sigma=gammas[[2]]/ni2)
    mus[[3]]=mvrnorm(n=1,mu=colSums(Y[indicators==3,]/ni3),Sigma=gammas[[3]]/ni3)
    mus[[4]]=mvrnorm(n=1,mu=colSums(Y[indicators==4,]/ni4),Sigma=gammas[[4]]/ni4)
    
    ########## Update Delta's:
    for (ii in ind)
    {
      zo=1*(fillnew[[ii]][,"zone"]=="PZ")+0*(fillnew[[ii]][,"zone"]=="CG")
      can=fillnew[[ii]][,"cancer"]
      indicator=1*((can==1)&(zo==1))+2*((can==1)&(zo==0))+3*((can==0)&(zo==1))+4*((can==0)&(zo==0))
      invgamsq=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])})
      V_deltai=solvecpp( matrix(rowSums(invgamsq),4,4) + solvecpp(Sigmasq))
      leftover=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])%*%(y[[ii]][x,]-mus[[indicator[x]]])})
      Mu_deltai=V_deltai%*%rowSums(leftover)
      Delta[[ii]]=mvrnorm(n=1,mu=Mu_deltai,Sigma=V_deltai)
    }
    # test Delta:
    zo=1*(fillnew[[slice]][,"zone"]=="PZ")+0*(fillnew[[slice]][,"zone"]=="CG")
    ####### SSE:
    can=cancertest_sse[s-1,]
    indicator=1*((can==1)&(zo==1))+2*((can==1)&(zo==0))+3*((can==0)&(zo==1))+4*((can==0)&(zo==0))
    invgamsq=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])})
    V_deltai=solvecpp( matrix(rowSums(invgamsq),4,4) + solvecpp(Sigmasq))
    leftover=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])%*%(y[[slice]][x,]-mus[[indicator[x]]])})
    Mu_deltai=V_deltai%*%rowSums(leftover)
    Delta_sse=mvrnorm(n=1,mu=Mu_deltai,Sigma=V_deltai)
    ####### SSE + NNGP:
    can=cancertest_ssenngp[s-1,]
    indicator=1*((can==1)&(zo==1))+2*((can==1)&(zo==0))+3*((can==0)&(zo==1))+4*((can==0)&(zo==0))
    invgamsq=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])})
    V_deltai=solvecpp( matrix(rowSums(invgamsq),4,4) + solvecpp(Sigmasq))
    leftover=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])%*%(y[[slice]][x,]-mus[[indicator[x]]])})
    Mu_deltai=V_deltai%*%rowSums(leftover)
    Delta_ssenngp=mvrnorm(n=1,mu=Mu_deltai,Sigma=V_deltai)
    ####### SSE + Low Rank:
    can=cancertest_sselr[s-1,]
    indicator=1*((can==1)&(zo==1))+2*((can==1)&(zo==0))+3*((can==0)&(zo==1))+4*((can==0)&(zo==0))
    invgamsq=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])})
    V_deltai=solvecpp( matrix(rowSums(invgamsq),4,4) + solvecpp(Sigmasq))
    leftover=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])%*%(y[[slice]][x,]-mus[[indicator[x]]])})
    Mu_deltai=V_deltai%*%rowSums(leftover)
    Delta_sselr=mvrnorm(n=1,mu=Mu_deltai,Sigma=V_deltai)
    ####### SSE + CAR:
    can=cancertest_ssecar[s-1,]
    indicator=1*((can==1)&(zo==1))+2*((can==1)&(zo==0))+3*((can==0)&(zo==1))+4*((can==0)&(zo==0))
    invgamsq=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])})
    V_deltai=solvecpp( matrix(rowSums(invgamsq),4,4) + solvecpp(Sigmasq))
    leftover=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])%*%(y[[slice]][x,]-mus[[indicator[x]]])})
    Mu_deltai=V_deltai%*%rowSums(leftover)
    Delta_ssecar=mvrnorm(n=1,mu=Mu_deltai,Sigma=V_deltai)
    ####### SSE + Full Spatial Model:
    can=cancertest_full[s-1,]
    indicator=1*((can==1)&(zo==1))+2*((can==1)&(zo==0))+3*((can==0)&(zo==1))+4*((can==0)&(zo==0))
    invgamsq=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])})
    V_deltai=solvecpp( matrix(rowSums(invgamsq),4,4) + solvecpp(Sigmasq))
    leftover=sapply(1:length(indicator),FUN=function(x){solvecpp(gammas[[indicator[x]]])%*%(y[[slice]][x,]-mus[[indicator[x]]])})
    Mu_deltai=V_deltai%*%rowSums(leftover)
    Delta_full=mvrnorm(n=1,mu=Mu_deltai,Sigma=V_deltai)
    ##### Update Sigmasq:
    S_Delta=sapply(ind,FUN=function(x){t(t(Delta[[x]]))%*%t(Delta[[x]])})
    Sigmasq=riwish(v=34+delta,S=matrix(rowSums(S_Delta),4,4)+omega.delta[[1]])
    
    ##### Update q with NNGP:
    ##### Update q:
    wijresults=nngpwupdatecpp(sigmasq, tausq, lambda, FF.update, cancertest_ssenngp[s-1,], wupdate_ssenngp, d, n, B.update,
                                                          nIndx, nnIndx, nnIndxLU, uIndx, uIndxLU, uiIndx, kupdate_ssenngp)
    wupdate_ssenngp=wijresults$wij
    ###### Update kij:
    kijresults = nngpkijupdatecpp(n, lambda, wupdate_ssenngp, cancertest_ssenngp[s-1,], tausq)
    kupdate_ssenngp=kijresults$kij
    #qiupdateresult_ssenngp=qiupdatecpp(sigmasq,tausq,lambda,FF.update,cancertest_ssenngp[s-1,],qupdate_ssenngp[s-1,],
    #                                   d,n,B.update,nIndx,nnIndx,nnIndxLU,
    #                                   uIndx,uIndxLU,uiIndx)
    #qupdate_ssenngp[s,]=qiupdateresult_ssenngp$q
    
    ##### Update q with Low Rank:
    wiupdateresult=qiupdate_lrcpp(sigmasq,tausq,lambda,cancertest_sselr[s-1,],wupdate_sselr[s-1,],
                                  n,mecoef,va)
    wupdate_sselr[s,]=wiupdateresult$w
    
    ##### Update q with CAR:
    qiupdateresult_ssecar=qicarcpp(sigmasq,lambda,cancertest_ssecar[s-1,],qupdate_ssecar[s-1,],Dinv)
    qupdate_ssecar[s,]=qiupdateresult_ssecar$q
    
    ##### Update q with Full model:
    qiupdateresult_full=qiupdate_fullcpp(q0=lambda,cancer=cancertest_full[s-1,],w=qupdate_full[s-1,],n=n,meanfull=meanfull,varfull=varfull)
    qupdate_full[s,]=qiupdateresult_full$w
    
    ######## Update auc's:
    btru=test[,"cancer"]
    ###### SSE Only:
    cancertest_sse[s,] = sapply(1:n,FUN=function(x){rbinomcpp(zone[x],y[[slice]][x,],means,Delta_sse,covs,pcancers[x])}) #Rcpp 2.986 s -> 0.073s 
    bpre_sse=colMeans(cancertest_sse[-1,])
    bpred_sse <- prediction(bpre_sse, btru)
    bperf_sse<- performance(bpred_sse,"tpr","fpr")
    bau_sse <- performance(bpred_sse,"auc")
    auc_sse_temp[s] <- unlist(slot(bau_sse, "y.values"))
    ###### SSE + NNGP:
    cancertest_ssenngp[s,] = sapply(1:n,FUN=function(x){rbinomcpp(zone[x],y[[slice]][x,],mus,Delta_ssenngp,gammas,kijresults$pstar[x])}) #Rcpp 2.986 s -> 0.073s 
    #cancertest_ssenngp[s,] = sapply(1:n,FUN=function(x){rbinomcpp(zone[x],y[[slice]][x,],mus,Delta_ssenngp,gammas,qiupdateresult_ssenngp$pstar[x])}) #Rcpp 2.986 s -> 0.073s 
    bpre_ssenngp=colMeans(cancertest_ssenngp[1:s,])
    bpred_ssenngp <- prediction(bpre_ssenngp, btru)
    bperf_ssenngp<- performance(bpred_ssenngp,"tpr","fpr")
    bau_ssenngp <- performance(bpred_ssenngp,"auc")
    auc_ssenngp_temp[s] <- unlist(slot(bau_ssenngp, "y.values"))
    ###### SSE + Low Rank:
    cancertest_sselr[s,] = sapply(1:n,FUN=function(x){rbinomcpp(zone[x],y[[slice]][x,],mus,Delta_sselr,gammas,wiupdateresult$pstar[x])}) #Rcpp 2.986 s -> 0.073s 
    bpre_sselr=colMeans(cancertest_sselr[1:s,])
    bpred_sselr <- prediction(bpre_sselr, btru)
    bperf_sselr<- performance(bpred_sselr,"tpr","fpr")
    bau_sselr <- performance(bpred_sselr,"auc")
    auc_sselr_temp[s] <- unlist(slot(bau_sselr, "y.values"))
    ###### SSE + CAR:
    cancertest_ssecar[s,] = sapply(1:n,FUN=function(x){rbinomcpp(zone[x],y[[slice]][x,],mus,Delta_ssecar,gammas,qiupdateresult_ssecar$pstar[x])}) #Rcpp 2.986 s -> 0.073s 
    bpre_ssecar=colMeans(cancertest_ssecar[1:s,])
    bpred_ssecar <- prediction(bpre_ssecar, btru)
    bperf_ssecar<- performance(bpred_ssecar,"tpr","fpr")
    bau_ssecar <- performance(bpred_ssecar,"auc")
    auc_ssecar_temp[s] <- unlist(slot(bau_ssecar, "y.values"))
    ###### SSE + Full:
    cancertest_full[s,] = sapply(1:n,FUN=function(x){rbinomcpp(zone[x],y[[slice]][x,],mus,Delta_full,gammas,qiupdateresult_full$pstar[x])}) #Rcpp 2.986 s -> 0.073s 
    bpre_full=colMeans(cancertest_full[1:s,])
    bpred_full <- prediction(bpre_full, btru)
    bperf_full<- performance(bpred_full,"tpr","fpr")
    bau_full <- performance(bpred_full,"auc")
    auc_full_temp[s] <- unlist(slot(bau_full, "y.values"))
    #print(set)
    #print(slice)
    #print(s)
  }
  btru=test[,"cancer"]
  ###### SSE:
  bpre_sse=colMeans(cancertest_sse[201:s,])
  bpred_sse <- prediction(bpre_sse, btru)
  bperf_sse<- performance(bpred_sse,"tpr","fpr")
  bau_sse <- performance(bpred_sse,"auc")
  tp_sse[ss]=max(bperf_sse@y.values[[1]][bperf_sse@x.values[[1]]<=0.1])
  auc_sse[ss]=auc_sse_temp[s]
  ###### SSE + NNGP:
  bpre_ssenngp=colMeans(cancertest_ssenngp[201:s,])
  bpred_ssenngp <- prediction(bpre_ssenngp, btru)
  bper_ssenngpf<- performance(bpred_ssenngp,"tpr","fpr")
  bau_ssenngp <- performance(bpred_ssenngp,"auc")
  tp_ssenngp[ss]=max(bperf_ssenngp@y.values[[1]][bperf_ssenngp@x.values[[1]]<=0.1])
  auc_ssenngp[ss]=auc_ssenngp_temp[s]
  ###### SSE + Low Rank:
  bpre_sselr=colMeans(cancertest_sselr[201:s,])
  bpred_sselr <- prediction(bpre_sselr, btru)
  bperf_sselr<- performance(bpred_sselr,"tpr","fpr")
  bau_sselr <- performance(bpred_sselr,"auc")
  tp_sselr[ss]=max(bperf_sselr@y.values[[1]][bperf_sselr@x.values[[1]]<=0.1])
  auc_sselr[ss]=auc_sselr_temp[s]
  ###### SSE + CAR:
  bpre_ssecar=colMeans(cancertest_ssecar[201:s,])
  bpred_ssecar <- prediction(bpre_ssecar, btru)
  bperf_ssecar<- performance(bpred_ssecar,"tpr","fpr")
  bau_ssecar <- performance(bpred_ssecar,"auc")
  tp_ssecar[ss]=max(bperf_ssecar@y.values[[1]][bperf_ssecar@x.values[[1]]<=0.1])
  auc_ssecar[ss]=auc_ssecar_temp[s]
  ###### SSE + Full:
  bpre_full=colMeans(cancertest_full[201:s,])
  bpred_full <- prediction(bpre_full, btru)
  bperf_full<- performance(bpred_full,"tpr","fpr")
  bau_full <- performance(bpred_full,"auc")
  tp_full[ss]=max(bperf_full@y.values[[1]][bperf_full@x.values[[1]]<=0.1])
  auc_full[ss]=auc_full_temp[s]
  #print(auc_sprem[[set]][ss])
  ss=ss+1
}
tpsse[S]=mean(tp_sse)
aucsse[S]=mean(auc_sse)
tpssenngp[S]=mean(tp_ssenngp)
aucssenngp[S]=mean(auc_ssenngp)
tpsselr[S]=mean(tp_sselr)
aucsselr[S]=mean(auc_sselr)
tpssecar[S]=mean(tp_ssecar)
aucssecar[S]=mean(auc_ssecar)
tpfull[S]=mean(tp_full)
aucfull[S]=mean(auc_full)
}
