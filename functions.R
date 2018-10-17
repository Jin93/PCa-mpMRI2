####### Define functions:
dist2=function(a1,a2,b1,b2)
{
  sqrt((a1-b1)^2+(a2-b2)^2)
}
getNNIndx=function(i, m){ # i=1:n
  
  if(i == 1){
    iNNIndx = 0;#this should never be accessed
    iNN = 0;
    return(list(iNNIndx=iNNIndx,iNN=iNN));
  }else if(i <= m){
    iNNIndx = (i)/2*(i-1);
    iNN = i-1;
    return(list(iNNIndx=iNNIndx,iNN=iNN));
  }else{
    iNNIndx = (m)/2*(m-1)+(i-m)*m;
    iNN = m;
    return(list(iNNIndx=iNNIndx,iNN=iNN));
  } 
}

# Allocated for the nearest neighbor index vector (note, first location has no neighbors):
#nnDist=d
mkNNIndxTree0=function(n, m, coords,nIndx) #output: nnIndx, nnDist, nnIndxLU
{
  nnIndx=nnDist=rep(NA,nIndx)
  nnIndxLU=rep(NA,2*n)
  nnDist=rep(NA,nIndx)
  for (i in 2:m)
  {
    #nnIndx[((i-1)*(i-2)/2+1):(i*(i-1)/2)]=i
    for (j in 1:(i-1))
    {
      nnDist[((i-1)*(i-2)/2+1):(i*(i-1)/2)][j]=dist2(coords[i,1],coords[i,2],coords[j,1],coords[j,2])
      nnIndx[((i-1)*(i-2)/2+1):(i*(i-1)/2)][j]=j
    }
    nnIndxLU[i] = (i-1)/2*(i-2); # the location BEFORE the starting location!!
    nnIndxLU[n+i] = i-1
  }
  ###### i>m: nearest neighbor:
  for (i in (m+1):n)
  {
    nnIndx[(m*(m-1)/2+1+(i-m-1)*m):(m*(m-1)/2+(i-m)*m)]=i
    ### get the row indicator in coords of the n.neighbors:
    dis=sapply(1:(i-1),FUN=function(x){dist2(coords[i,1],coords[i,2],coords[x,1],coords[x,2])})
    indx=order(dis)[1:m]
    nnDist[(m*(m-1)/2+1+(i-m-1)*m):(m*(m-1)/2+(i-m)*m)]=dis[indx]
    nnIndx[(m*(m-1)/2+1+(i-m-1)*m):(m*(m-1)/2+(i-m)*m)]=indx
    nnIndxLU[i] = m*(m-1)/2+(i-m-1)*m # the location BEFORE the starting location!!
    nnIndxLU[n+i] = m
    #print(i)#40 seconds, can do parallel
  }
  return(list(nnIndx=nnIndx,nnDist=nnDist,nnIndxLU=nnIndxLU))
}
  
  
##### Make U Index:
mkUIndx=function(n, m, nIndx,nnIndx,nnIndxLU)
{
  uIndx=rep(NA,nIndx)
  uIndxLU=rep(NA,2*n) 
  ### 1st column: the location BEFORE the starting location!!
  ### 2nd column: how many locations have i as one of the neighbors~
  l=0
  for (i in 1:n)
  {
    uIndxLU[i] = l
    h=1
    for (j in 2:n)
    {
      jNNIndx=nnIndxLU[j]
      jNN=nnIndxLU[j+n]
      for (k in 1:jNN)
      {
        if (nnIndx[jNNIndx+k]== i)
        {
          uIndx[l+h]=j
          h=h+1
        }
      }
    }
    l=l+h-1
    uIndxLU[n+i]=h-1 #important!!!
  }
  return(list(uIndx=uIndx,uIndxLU=uIndxLU))
}
#### 3434, 3473, 3785: no location have these as neighbors

### Description: update B and F.
updateBF=function(d, D,nIndx,nnIndxLU, CIndx, n, sigmaSq, tausq, phi, nu)#int covModel
{
  c=B=rep(NA,nIndx)
  C=rep(0,length(D))
  FF=rep(NA,n)
  FF[1]=sigmaSq+tausq
  ### C_(sij,sij): diag(sigma^2,n)
  ### c: c_(sij,N(sij)):
  for(i in 2:n)
  {
    ## nnIndxLU[n+i] is always > 0 except when i=1: NA
    for (k in 1:nnIndxLU[n+i])
    {
      dd=d[nnIndxLU[i]+k]
      c[nnIndxLU[i]+k] = sigmaSq*exp(-(dd*phi)) + (dd == 0) * tausq
      #c[nnIndxLU[i]+k] = sigmaSq*exp(-dd/phi)
      for(l in 1:nnIndxLU[n+i])
      {
        ddd=D[CIndx[i]+(l-1)*nnIndxLU[n+i]+k]
        C[CIndx[i]+(l-1)*nnIndxLU[n+i]+k] = sigmaSq*exp(-(ddd*phi)) + (dd == 0) * tausq
      }
    }
    temC=matrix(C[(CIndx[i]+1):(CIndx[i]+CIndx[n+i])],nnIndxLU[n+i],nnIndxLU[n+i])
    Cinv=chol2inv(chol(temC))
    B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv
    FF[i]=sigmaSq + tausq -B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%t(t(c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]))
    C[(CIndx[i]+1):(CIndx[i]+CIndx[n+i])]=Cinv
    #B[i]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv: a vector!!!
    #print(i)#:25 seconds
  }
  return(list(c=c,B=B,C=C,FF=FF))
}


######## Spherical:
updateBF_spherical=function(d, D,nIndx,nnIndxLU, CIndx, n, sigmaSq, tausq, phi, nu)#int covModel
{
  c=B=rep(NA,nIndx)
  C=rep(0,length(D))
  FF=rep(NA,n)
  FF[1]=sigmaSq + tausq
  ### C_(sij,sij): diag(sigma^2,n)
  ### c: c_(sij,N(sij)):
  for(i in 2:n)
  {
    ## nnIndxLU[n+i] is always > 0 except when i=1: NA
    for (k in 1:nnIndxLU[n+i])
    {
      dd=d[nnIndxLU[i]+k]
      c[nnIndxLU[i]+k] = sigmaSq*(1-(dd*phi)*2/3 + (dd*phi)^3/2 ) + (dd == 0) * tausq
      for(l in 1:nnIndxLU[n+i])
      {
        ddd=D[CIndx[i]+(l-1)*nnIndxLU[n+i]+k]
        C[CIndx[i]+(l-1)*nnIndxLU[n+i]+k] = sigmaSq*(1-(ddd*phi)*2/3 + (ddd*phi)^3/2 ) + (dd == 0) * tausq
      }
    }
    temC=matrix(C[(CIndx[i]+1):(CIndx[i]+CIndx[n+i])],nnIndxLU[n+i],nnIndxLU[n+i])
    Cinv=chol2inv(chol(temC))
    B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv
    FF[i]=sigmaSq + tausq -B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%t(t(c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]))
    C[(CIndx[i]+1):(CIndx[i]+CIndx[n+i])]=Cinv
    #B[i]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv: a vector!!!
    #print(i)#:25 seconds
  }
  return(list(c=c,B=B,C=C,FF=FF))
}


######## Matern:
updateBF_matern=function(d, D,nIndx,nnIndxLU, CIndx, n, sigmaSq, tausq, phi, nu)#int covModel
{
  c=B=rep(NA,nIndx)
  C=rep(0,length(D))
  FF=rep(NA,n)
  FF[1]=sigmaSq + tausq
  ### C_(sij,sij): diag(sigma^2,n)
  ### c: c_(sij,N(sij)):
  for(i in 2:n)
  {
    ## nnIndxLU[n+i] is always > 0 except when i=1: NA
    for (k in 1:nnIndxLU[n+i])
    {
      dd=d[nnIndxLU[i]+k]
      if (dd == 0)
      {
        c[nnIndxLU[i]+k] = sigmaSq + tausq
      }
      if (dd > 0)
      {
        c[nnIndxLU[i]+k] = sigmaSq*Matern(d=dd,alpha=phi,nu=nu)
      }
      for(l in 1:nnIndxLU[n+i])
      {
        ddd=D[CIndx[i]+(l-1)*nnIndxLU[n+i]+k]
        if (ddd == 0) { C[CIndx[i]+(l-1)*nnIndxLU[n+i]+k] = sigmaSq + tausq }
        if (ddd > 0) { C[CIndx[i]+(l-1)*nnIndxLU[n+i]+k] = sigmaSq*Matern(d=ddd,alpha=phi,nu=nu) }
      }
    }
    temC=matrix(C[(CIndx[i]+1):(CIndx[i]+CIndx[n+i])],nnIndxLU[n+i],nnIndxLU[n+i])
    Cinv=chol2inv(chol(temC))
    B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv
    FF[i]=sigmaSq + tausq -B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%t(t(c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]))
    C[(CIndx[i]+1):(CIndx[i]+CIndx[n+i])]=Cinv
    #B[i]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv: a vector!!!
    #print(i)#:25 seconds
  }
  return(list(c=c,B=B,C=C,FF=FF))
}


##### Gaussian:
updateBF_gaussian=function(d, D,nIndx,nnIndxLU, CIndx, n, sigmaSq, tausq, phi, nu)#int covModel
{
  c=B=rep(NA,nIndx)
  C=rep(0,length(D))
  FF=rep(NA,n)
  FF[1]=sigmaSq + tausq
  ### C_(sij,sij): diag(sigma^2,n)
  ### c: c_(sij,N(sij)):
  for(i in 2:n)
  {
    ## nnIndxLU[n+i] is always > 0 except when i=1: NA
    for (k in 1:nnIndxLU[n+i])
    {
      dd=d[nnIndxLU[i]+k]
      #c[nnIndxLU[i]+k] = sigmaSq*Matern(dd,range=1/phi,smoothness=nu)
      c[nnIndxLU[i]+k] = sigmaSq*exp(-((dd*phi)^2)) + (dd == 0) * tausq
      #c[nnIndxLU[i]+k] = sigmaSq*exp(-dd/phi)
      for(l in 1:nnIndxLU[n+i])
      {
        ddd=D[CIndx[i]+(l-1)*nnIndxLU[n+i]+k]
        #C[CIndx[i]+(l-1)*nnIndxLU[n+i]+k] = sigmaSq*Matern(ddd,range=1/phi,smoothness=nu)
        C[CIndx[i]+(l-1)*nnIndxLU[n+i]+k] = sigmaSq*exp(-((ddd*phi)^2)) + (dd == 0) * tausq
        #sigmaSq*spCor(D[CIndx[i]+l*nnIndxLU[n+i]+k], phi, nu, covModel, &bk[threadID*nb]); 
      }
    }
    temC=matrix(C[(CIndx[i]+1):(CIndx[i]+CIndx[n+i])],nnIndxLU[n+i],nnIndxLU[n+i])
    Cinv=chol2inv(chol(temC))
    B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv
    FF[i]=sigmaSq + tausq -B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%t(t(c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]))
    C[(CIndx[i]+1):(CIndx[i]+CIndx[n+i])]=Cinv
    #B[i]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv: a vector!!!
    #print(i)#:25 seconds
  }
  return(list(c=c,B=B,C=C,FF=FF))
}



######## Block M-H h:
hf=function(n,w,c,C,FF,nnIndx,nnIndxLU,CIndx,sigmaSq,tausq)
{
  #### 1st pixel:
  tem1=dnorm(w[1],mean=0,sd=sqrt(sigmaSq+tausq))
  tem2=sapply(2:n,FUN=function(xx){dnorm(w[xx],mean=c[(nnIndxLU[xx]+1):(nnIndxLU[xx]+nnIndxLU[n+xx])]%*%matrix(C[(CIndx[xx]+1):(CIndx[xx]+CIndx[n+xx])],nnIndxLU[n+xx],nnIndxLU[n+xx])%*%t(t(w[nnIndx[(nnIndxLU[xx]+1):(nnIndxLU[xx]+nnIndxLU[n+xx])]])),
                                            sd=sqrt(FF[xx]))})
  c(tem1,tem2)
}

hfnngp=function(n,w,c,C,FF,nnIndx,nnIndxLU,CIndx,sigmaSq)
{
  #### 1st pixel:
  tem1=dnorm(w[1],mean=0,sd=sqrt(sigmaSq))
  tem2=sapply(2:n,FUN=function(xx){dnorm(w[xx],mean=c[(nnIndxLU[xx]+1):(nnIndxLU[xx]+nnIndxLU[n+xx])]%*%matrix(C[(CIndx[xx]+1):(CIndx[xx]+CIndx[n+xx])],nnIndxLU[n+xx],nnIndxLU[n+xx])%*%t(t(w[nnIndx[(nnIndxLU[xx]+1):(nnIndxLU[xx]+nnIndxLU[n+xx])]])),
                                         sd=sqrt(FF[xx]))})
  c(tem1,tem2)
}

mixture.sampling <- function(M,G,O,V,w1,w2)
{
  qtemp=rnorm(n=1,mean=M,sd=sqrt(G))
  pa=(w1*dnorm(x=qtemp,mean=M,sd=sqrt(G)) + w2*dnorm(x=qtemp,mean=O,sd=sqrt(V)))/(w1*dnorm(x=qtemp,mean=M,sd=sqrt(G)))
  u=runif(n=1,min=0,max=1)
  while (u>pa) {
    qtemp=rnorm(n=1,mean=M,sd=sqrt(G))
    pa=(w1*dnorm(x=qtemp,mean=M,sd=sqrt(G)) + w2*dnorm(x=qtemp,mean=O,sd=sqrt(V)))/(w1*dnorm(x=qtemp,mean=M,sd=sqrt(G)))
    u=runif(n=1,min=0,max=1)
  }
  return(qtemp)
}

mixture.sampling.MG <- function(M,G)
{
  V=1/(1+1/G)
  O=V*M/G
  w1=1/( 1-(exp(-0.5*((M^2)/G - (O^2/V) ))/sqrt(1+G)) )
  w2=1-w1  
  qtemp=rnorm(n=1,mean=M,sd=sqrt(G))
  pa=(w1*dnorm(x=qtemp,mean=M,sd=sqrt(G)) + w2*dnorm(x=qtemp,mean=O,sd=sqrt(V)))/(w1*dnorm(x=qtemp,mean=M,sd=sqrt(G)))
  u=runif(n=1,min=0,max=1)
  while (u>pa) {
    qtemp=rnorm(n=1,mean=M,sd=sqrt(G))
    pa=(w1*dnorm(x=qtemp,mean=M,sd=sqrt(G)) + w2*dnorm(x=qtemp,mean=O,sd=sqrt(V)))/(w1*dnorm(x=qtemp,mean=M,sd=sqrt(G)))
    u=runif(n=1,min=0,max=1)
  }
  return(qtemp)
}


######### Block Gibbs Sampling for w_i:
qiupdate=function(sigmasq,tausq,lambda,FF,cancer,q,d,n,B,nIndx,nnIndx,nnIndxLU,uIndx,uIndxLU,uiIndx)
{
  m=M=G=V=O=pstar=numeric()
  #### 1st element:
  ii=1
  m[ii]=0
  G[ii]=1/(sigmasq+tausq)
  for (jj in 1:uIndxLU[n+ii])
  {
    t = uIndx[uIndxLU[ii]+jj]
    l=uiIndx[uIndxLU[ii]+jj]
    #a=q[t]-B[((nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t]))[-l]]%*%t(t(q[(nnIndx[(nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t])])[-l]]))
    a=q[t]-sum(B[((nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t]))[-l]]*q[(nnIndx[(nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t])])[-l]])
    m[ii]=m[ii]+B[nnIndxLU[t]+l]*a/FF[t]
    G[ii]=G[ii]+(B[nnIndxLU[t]+l]^2)/FF[t]
  }
  G[ii]=1/G[ii]
  M[ii]=m[ii]*G[ii]
  #V[ii]=1/(1+1/G[ii])
  #O[ii]=V[ii]*M[ii]/G[ii]
  if (cancer[ii]==1)
  {
    # Truncated Normal:
    q[ii]=rtruncnorm(n=1,a=-lambda[ii],b=Inf,mean=M[ii],sd=sqrt(G[ii]))
  }
  if (cancer[ii]==0)
  {
    q[ii]=rtruncnorm(n=1,a=-Inf,b=-lambda[ii],mean=M[ii],sd=sqrt(G[ii]))
  }
  pstar[ii]=1-pnorm(q=-lambda[ii],mean=M[ii],sd=sqrt(G[ii]))
  ##### 2~nth element: Cholesky decomposition:
  for (ii in 2:n)
  {
    if (uIndxLU[n+ii] == 0)
    {
      m[ii] = (1/FF[ii])* (B[(nnIndxLU[ii]+1):(nnIndxLU[ii]+nnIndxLU[n+ii])]%*%q[nnIndx[(nnIndxLU[ii]+1):(nnIndxLU[ii]+nnIndxLU[n+ii])]])
      G[ii] = FF[ii]
      M[ii]=m[ii]*G[ii]
      #V[ii]=1/(1+1/G[ii])
      #O[ii]=V[ii]*M[ii]/G[ii]
      if (cancer[ii]==1)
      {
        # Truncated Normal:
        q[ii]=rtruncnorm(n=1,a=-lambda[ii],b=Inf,mean=M[ii],sd=sqrt(G[ii]))
      }
      if (cancer[ii]==0)
      {
        q[ii]=rtruncnorm(n=1,a=-Inf,b=-lambda[ii],mean=M[ii],sd=sqrt(G[ii]))
      }
      pstar[ii]=1-pnorm(q=-lambda[ii],mean=M[ii],sd=sqrt(G[ii]))
    }
    ############
    if (uIndxLU[n+ii] >0)
    {
      m[ii] = (1/FF[ii])* B[(nnIndxLU[ii]+1):(nnIndxLU[ii]+nnIndxLU[n+ii])]%*%q[nnIndx[(nnIndxLU[ii]+1):(nnIndxLU[ii]+nnIndxLU[n+ii])]]
      G[ii] = 1/FF[ii]
      for (jj in 1:uIndxLU[n+ii])
      {
        t = uIndx[uIndxLU[ii]+jj]
        l=uiIndx[uIndxLU[ii]+jj]
        #a=q[t]-B[((nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t]))[-l]]%*%t(t(q[(nnIndx[(nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t])])[-l]]))
        a=q[t]-sum(B[((nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t]))[-l]]*q[(nnIndx[(nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t])])[-l]])
        m[ii]=m[ii]+B[nnIndxLU[t]+l]*a/FF[t]
        G[ii]=G[ii]+(B[nnIndxLU[t]+l]^2)/FF[t]
      }
      G[ii]=1/G[ii]
      M[ii]=m[ii]*G[ii]
      #V[ii]=1/(1+1/G[ii])
      #O[ii]=V[ii]*M[ii]/G[ii]
      if (cancer[ii]==1)
      {
        # Truncated Normal:
        q[ii]=rtruncnorm(n=1,a=-lambda[ii],b=Inf,mean=M[ii],sd=sqrt(G[ii]))
      }
      if (cancer[ii]==0)
      {
        q[ii]=rtruncnorm(n=1,a=-Inf,b=-lambda[ii],mean=M[ii],sd=sqrt(G[ii]))
      }
      pstar[ii]=1-pnorm(q=-lambda[ii],mean=M[ii],sd=sqrt(G[ii]))
    }
  }
  return(list(q=q,pstar=pstar,M=M,G=G))
}
######## update pstars:
updatepstar=function(M,G,q0)
{
  1-pnorm(q=-q0,mean=M,sd=sqrt(G))
}

calstar=function(n,nnIndx,nnIndxLU,B.update,w,sigmasq,tausq,FF.update)
{
  star=w[1]^2/(sigmasq+tausq)
  for (i in 2:n)
  {
    e=0
    for(j in 1:nnIndxLU[n+i])
    {
      e=e+B.update[nnIndxLU[i]+j]*w[nnIndx[nnIndxLU[i]+j]]
    }
    b=w[i]-e
    star=star+b^2/FF.update[i]
  }
  return(list(star=star))
}

########### Calculate B matrix for CAR model:
updateBi=function(x)
{
  coords=fillnew[[x]][,c("x","y")]
  Dinv <- 1/as.matrix(dist(coords))
  for (x in 1:nrow(Dinv))
  {
    Dinv[x,x]=0
  }
  Dinv <- Dinv/rowSums(Dinv)
  return(Dinv)
}


######### Updating q with CAR model: (old)
qicar=function(sigmasq,lambda,cancer,q,Dinv)
{ 
  pstar=numeric()
  for (x in 1:length(q))
  { 
    qmean=Dinv[x,-x]%*%q[-x]
    if (cancer[x]==1)
    { 
      # Truncated Normal:
      q[x]=rtruncnorm(n=1,a=-lambda[x],b=Inf,mean=qmean,sd=sqrt(sigmasq))
    }
    if (cancer[x]==0)
    { 
      q[x]=rtruncnorm(n=1,a=-Inf,b=-lambda[x],mean=qmean,sd=sqrt(sigmasq))
    }
    pstar[x]=1-pnorm(q=-lambda[x],mean=qmean,sd=sqrt(sigmasq))
  }
  return(list(q=q,pstar=pstar))
}


######### Updating w with CAR model: (new)
wicar=function(sigmasq,lambda,q,Dinv)
{ 
  pstar=numeric()
  for (x in 1:length(q))
  { 
    qmean=Dinv[x,-x]%*%q[-x]
    q[x]=rtruncnorm(n=1,mean=qmean,sd=sqrt(sigmasq))
    pstar[x]=1-pnorm(q=-lambda[x],mean=qmean,sd=sqrt(sigmasq))
  }
  return(list(q=q))
}



###### Cancer full conditional:
piupdate_nodelta=function(zon,y,mu,cov,pstar)
{
  indicat1=1*(zon==1)+2*(zon==0)
  indicat0=3*(zon==1)+4*(zon==0)
  p1=dmvnorm(x=y,mean=mu[[indicat1]],sigma=cov[[indicat1]])*pstar
  p0=dmvnorm(x=y,mean=mu[[indicat0]],sigma=cov[[indicat0]])*(1-pstar)
  ppcancer=p1/(p0+p1)
  return(list(ppcancer=ppcancer))
}


piupdate=function(zon,y,mu,Delta,cov,pstar)
{
  indicat1=1*(zon==1)+2*(zon==0)
  indicat0=3*(zon==1)+4*(zon==0)
  p1=dmvnorm(x=y,mean=mu[[indicat1]]+Delta,sigma=cov[[indicat1]])*pstar
  p0=dmvnorm(x=y,mean=mu[[indicat0]]+Delta,sigma=cov[[indicat0]])*(1-pstar)
  ppcancer=p1/(p0+p1)
  return(list(ppcancer=ppcancer))
}




qstar=function(n,nnIndx,nnIndxLU,B.update,q,sigmasq,FF.update)
{
  star=q[1]^2/sigmasq
  for (i in 2:n)
  {
    e=0
    for(j in 1:nnIndxLU[n+i])
    {
      e=e+B.update[nnIndxLU[i]+j]*q[nnIndx[nnIndxLU[i]+j]]
    }
    b=q[i]-e
    star=star+b^2/FF.update[i]
  }
  return(list(star=star))
}

q0update=function(Z,zon,pp1,pp0,cancer){
  ZONE=Z*zon + (1-Z)*(1-zon)
  lik1=((1-pp1)*(1-cancer)+pp1*cancer)[ZONE!=0]
  lik0=((1-pp0)*(1-cancer)+pp0*cancer)[ZONE!=0] 
  #lik=((1-pp)*(1-cancer)+pp*cancer)*ZONE
  prod(lik1[lik0!=0]/lik0[lik0!=0])
}

#q0update=function(Z,zon,pp1,pp0,cancer){
#  ZONE=Z*zon + (1-Z)*(1-zon)
#  lik1=((1-pp1)*(1-cancer)+pp1*cancer)[ZONE!=0]
#  lik0=((1-pp0)*(1-cancer)+pp0*cancer)[ZONE!=0]
#  #lik=((1-pp)*(1-cancer)+pp*cancer)*ZONE
#  prod(lik1/lik0)
#}
####################### Low Rank:
######## Block M-H for Theta:
hflr=function(n,w,ci,invCistar,tausq)
{
  dmvnorm(w,mean=rep(0,n),sigma=tausq*diag(1,n) + t(ci)%*%invCistar%*%ci)
}

################ update ci etc:
updateDi=function(x,knots)
{
  coordsi=fillnew[[x]][,c("x","y")]
  Di=as.matrix(dist(coordsi))
  diknots=Di[knots[x,],knots[x,]]
  #DIknots=as.matrix(dist(coordsi[knots[x,],]))
  direst=Di[knots[x,],]
  return(list(diknots=diknots,direst=direst))
}

updateci=function(sigmasq,tausq,phi,nu,Diknots,Direst)
{
  Cistar=materncovcpp(Diknots,sigmasq,phi,nu) 
  invCistar=solvecpp(Cistar)
  ci = materncovcpp(Direst,sigmasq,phi,nu) 
  #covmat=t(ci)%*%invCistar%*%ci
  ciinvCistar=t(ci)%*%invCistar
  return(list(ci=ci,Cistar=Cistar,ciinvCistar=ciinvCistar,invCistar=invCistar))
}

################# update mean, va and pstar of kij:
# updatelrkcpp(double sigmaSq, double tausq, vec lambda, vec cancer, vec wstar, int n, mat ci, mat ciinvCistar)




updatemecoef=function(x,n,ci,Cistar,covmat,invCistar,tausq)
{
  mecoef=matrix(NA,n,(n-1))
  va=numeric()
  diagmat=diag(1/tausq,(n-1))
  mecoefresult=mclapply(1:n,FUN=function(ii){update2(x,n,ci,Cistar,covmat,invCistar,tausq,diagmat)},mc.cores=n.cores2)
  for (ii in 1:n)
  {
    mecoef[ii,]=mecoefresult[[ii]]$me
  }
  va=sapply(1:n,FUN=function(ii){mecoefresult[[ii]]$va})
  return(list(mecoef=mecoef,va=va))
}

########## update mecoef and va layer 2:
################# update mecoef and va:
update2=function(ii,n,ci,Cistar,covmat,invCistar,tausq,diagmat)
{
    cijtilde=covmat[ii,-ii]
    ci_j=ci[,-ii]
    cijj=covmat[ii,ii]
    tempinv=solvecpp(Cistar+ci_j%*%t(ci_j)/tausq)
    tempbracket=cijtilde%*%(diagmat-(1/(tausq)^2)*t(ci_j)%*%tempinv%*%ci_j)
    me=tempbracket
    va=cijj-tempbracket%*%t(t(cijtilde))
    #print(ii)
  return(list(me=me,va=va))
}

######### Block Gibbs Sampling for w_i:
qiupdate_lr=function(sigmasq,tausq,q0,cancer,w,n,mecoef,va)
{
  pstar=numeric()
  M=numeric()
  for (ii in 1:n)
  {
    me=mecoef[ii,]%*%(w[-ii]-0)
    if (cancer[ii]==1)
    {
      w[ii]=rtruncnorm(n=1,a=-q0[ii],b=Inf,mean=me,sd=sqrt(va[ii]))
    }
    if (cancer[ii]==0)
    {
      w[ii]=rtruncnorm(n=1,a=-Inf,b=-q0[ii],mean=me,sd=sqrt(va[ii]))
    }
    pstar[ii]=1-pnorm(q=-q0[ii],mean=me,sd=sqrt(va[ii]))
    M[ii]=me
  }
  return(list(w=w,pstar=pstar,M=M,G=va))
}

######### Full model - Block Gibbs Sampling for w_i:
qiupdate_full=function(q0,cancer,w,n,meanfull,varfull)
{
  pstar=numeric()
  M=G=numeric()
  for (ii in 1:n)
  {
    M[ii]=meanfull[ii,]%*%(w[-ii]-0)
    if (cancer[ii]==1)
    {
      w[ii]=rtruncnorm(n=1,a=-q0[ii],b=Inf,mean=M[ii],sd=sqrt(varfull[ii]))
    }
    if (cancer[ii]==0)
    {
      w[ii]=rtruncnorm(n=1,a=-Inf,b=-q0[ii],mean=M[ii],sd=sqrt(varfull[ii]))
    }
    pstar[ii]=1-pnorm(q=-q0[ii],mean=M[ii],sd=sqrt(varfull[ii]))
  }
  return(list(w=w,pstar=pstar,M=M,G=varfull))
}

####### Calculate Dinv:
calcdinv=function(ii)
{
  temp=fillnew[[ii]]
  coords=temp[,c("x","y")]
  Din <- 1/as.matrix(dist(coords))
  for (x in 1:nrow(Din))
  {
    Din[x,x]=0
  }
  Dinv<- Din/rowSums(Din)
  return(Dinv)
}

##### select knots:
selectknot=function(tem)
{
 coor=tem[,c("x","y")]
 xcoor=sort(unique(coor[,"x"]))
 ycoor=sort(unique(coor[,"y"]))
 xwidth=floor((length(xcoor)-1)/6)
 knots.x=xcoor[seq(1,length(xcoor),by=xwidth)]
 ywidth=floor((length(ycoor)-1)/6)
 knots.y=ycoor[seq(1,length(ycoor),by=ywidth)]
 Knots=numeric()
 knot=1
 for (a in 1:nrow(tem))
 {
   if ((tem[a,"x"]%in%knots.x)&(tem[a,"y"]%in%knots.y))
   {
     Knots[knot]=a
     knot=knot+1
   }
 }
 Knots=sample(Knots,size=n.knots,replace=FALSE)
 return(Knots)
}

updatelrwstar=function(tausq, cancer, lambda, kij, n, m, ci, Cistar, invCistar)
{
  Va=solvecpp(invCistar%*%(ci%*%t(ci))%*%invCistar/tausq + invCistar)
  Me=Va%*%invCistar%*% (ci%*%(kij-lambda))/tausq
  wstar=rmvnorm(1,mean=Me,sigma=Va)
  return(wstar);
  #return(list(wstar=wstar,Me=Me,Va=Va));
}

proptheta_lr=function(n,m,Cistar,invCistar,ci,Cistar1,invCistar1,ci1,kij,wstar,wstar1,lambda)
{
  a=dmvnorm(wstar1,mean=rep(0,m),sigma=Cistar1)/dmvnorm(wstar,mean=rep(0,m),sigma=Cistar)
  tem0=as.numeric(invCistar%*%t(wstar))
  tem1=as.numeric(invCistar1%*%t(wstar1))
  b=prod(sapply(1:n,FUN=function(xx){dnorm(kij[xx],mean=lambda[xx]+sum(as.numeric(ci1[,xx])*tem1),sd=1)/dnorm(kij[xx],mean=lambda[xx]+sum(as.numeric(ci[,xx])*tem0),sd=1)}))
  rat=a*b
  return(rat)
}


####### SVC functions:
vupdate=function(n, nx, w, B, FF, nnIndx, nnIndxLU, CIndx, iden)
{
  mats=list()
  mats[[1]] = t(t(w[1,])) %*% t(w[1,])/ FF[1];
  for (i in 2:n)
  {
    idx = nnIndx[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])];
    Btem = t(B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i]),1]) %x% iden
    wtem = matrix(t(w[idx,]),nx*length(idx),1)
    tempmean = Btem %*% wtem
    mats[[i]] = (t(t(w[i,])) - tempmean) %*% (t(w[i,]) - t(tempmean)) / FF[i];
  }
  return(mats);
}


svcmgupdate=function(n, nx, V, FF, FFmat, w, d, B, nIndx, nnIndx, nnIndxLU, uIndx, uIndxLU, uiIndx, iden)
{
  m=M=G=list()
  Gvec=numeric()
  Vinv = solvecpp(V);
  #1st element:
  ii=1;
  m[[ii]]=matrix(c(0,0,0),3,1)
  Gvec[ii]=1/FF[1];
  for (jj in 1:uIndxLU[n+ii])
  {
    t = uIndx[uIndxLU[ii]+jj];
    l = uiIndx[uIndxLU[ii]+jj];
    indx = nnIndx[(nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t])];
    wtem = matrix(t(w[indx,]),nx*length(indx),1);
    a = t(t(w[t,]))-t(B[(nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t]),1]) %x% iden %*% wtem
    Bscale = B[nnIndxLU[t]+l,1];
    BFinv = B[nnIndxLU[t]+l,1] / FF[t]
    BFinvkron = BFinv %x% Vinv
    a = a + (Bscale %x% iden) %*% t(t(w[indx[l],]))
    m[[ii]] = m[[ii]] + BFinvkron %*% a
    Gvec[ii] = Gvec[ii] + (B[nnIndxLU[t]+l,1])^2/FF[t]
  }
  Gscale = 1/Gvec[ii]
  G[[ii]] = Gscale %x% V
  M[[ii]] = G[[ii]] %*% m[[ii]]
  
    for (ii in 2:n)
    {
      if (uIndxLU[n+ii] == 0)
      {
        idx2 = nnIndx[(nnIndxLU[ii]+1):(nnIndxLU[ii]+nnIndxLU[n+ii])]
        wtem = matrix(t(w[idx2,]),nx*length(idx2),1)
        aa = t(B[(nnIndxLU[ii]+1):(nnIndxLU[ii]+nnIndxLU[n+ii]),1]) %x% iden %*% wtem
        Finv = solvecpp(FFmat[[ii]])
        m[[ii]] =  Finv %*% aa
        G[[ii]] = FFmat[[ii]]
        M[[ii]] = G[[ii]] %*% m[[ii]]
      }
      if (uIndxLU[n+ii] > 0)
      {
        idx2 = nnIndx[(nnIndxLU[ii]+1):(nnIndxLU[ii]+nnIndxLU[n+ii])]
        wtem = matrix(t(w[idx2,]),nx*length(idx2),1)
        aa = t(B[(nnIndxLU[ii]+1):(nnIndxLU[ii]+nnIndxLU[n+ii]),1]) %x% iden %*% wtem
        Finv = solvecpp(FFmat[[ii]])
        m[[ii]] =  Finv %*% aa
        Gvec[ii] = 1/FF[ii]
        for (jj in 1:uIndxLU[n+ii])
        {
          t = uIndx[uIndxLU[ii]+jj];
          l = uiIndx[uIndxLU[ii]+jj];
          indx = nnIndx[(nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t])]
          wtem = matrix(t(w[indx,]),nx*length(indx),1)
          a = t(t(w[t,]))-t(B[(nnIndxLU[t]+1):(nnIndxLU[t]+nnIndxLU[n+t]),1]) %x% iden %*% wtem
          Bscale = B[nnIndxLU[t]+l,1]
          BFinv = B[nnIndxLU[t]+l,1] / FF[t]
          BFinvkron = BFinv %x% Vinv
          a = a + (Bscale %x% iden) %*% t(t(w[indx[l],]))
          m[[ii]] = m[[ii]] + BFinvkron %*% a
          Gvec[ii] = Gvec[ii] + (B[nnIndxLU[t]+l,1])^2/FF[t]
        }
        Gscale = 1/Gvec[ii]
        G[[ii]] = Gscale %x% V
        M[[ii]] = G[[ii]] %*% m[[ii]]
      }
    }
  return (list(M=M, G=G))
}



svcwsupdate=function(n, nx, V, w, kij, M, G, X, loc_indx)
{
  postmean=postvar=list()
  for (ii in 1:n)
  {
    Xtemp = X[[ii]];
    Ginv = solvecpp(G[[ii]])
    Covtemp = solvecpp( t(Xtemp) %*% Xtemp + Ginv)
    locindex = loc_indx[[ii]]
    loc_index = locindex
    Mtem = M[[ii]];
    mtemp = t(Xtemp) %*% t(t(kij[loc_index])) + Ginv %*% Mtem
    postvar[[ii]] = Covtemp
    Mutemp = Covtemp %*% mtemp
    postmean[[ii]] = Mutemp
    w[ii,] = rmvnorm(1, Mutemp, Covtemp)
  }
  return(list(ws = w, postmean = postmean, postvar = postvar))
}


svckijupdate=function(nwhole, n, ws, cancer, Xij, indexij)
{
  kij=numeric()
  for (ii in 1:nwhole)
  {
    Xtemp = Xij[ii,];
    wij = ws[indexij[ii],]
    meanij = sum(Xtemp*ws)
    if (cancer[ii]==1)
    {
      kij[ii] = rtruncnorm(n=1,a=0,b=Inf,mean=meanij, sd=1);
    }
    if (cancer[ii]==0)
    {
      kij[ii] = rtruncnorm(n=1,a=-Inf,b=0,mean=meanij, sd=1);
    }
  }
  return(kij);
}
