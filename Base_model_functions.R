loghconst=function(delta, Omega){
  delta/2*as.numeric(determinant(Omega, logarithm = TRUE)$modulus)-0.5*nrow(Omega)*delta*log(2)- logGamma_p(0.5*delta, p=nrow(Omega))
}

logGamma_p=function(a, p){
  p*(p-1)/4*log(pi)+sum(lgamma(a+(1-1:p)/2))
}

####
logdet=function(Matrix)
{
  o=determinant(Matrix,logarithm=T)$modulus[1]
  return(o)
}

###################################################################################
# DataFit: Estimate the adjusted mean y.tilde and covariance S.tile in (7) and (8)#
###################################################################################
DataFit=function(dat){
  #################################################################################
  ##  return y.tilde and S.tilde
  ## correlation function implied by variable struc  specified at the beginning
  #################################################################################
  
  Npat=nrow(dat)
  mu.hat=vector(mode="list", 4)
  sigma.hat=vector(mode="list", 4)
  phidet.hat=vector(mode="list", 4)
  r=vector(mode="list", 4)
  Y.Y=list()
  for(z in 1:4){
    r.all=phidet=NULL
    y.tilde=rep(0, p)
    S.tilde=array(0,dim=c(p,p))
    
    for (i in 1:Npat)
    {
      keep=dat[[i,"region"]]==(z)
      n.i=sum(keep)
      ###################
      if(n.i>0){Y.i=as.matrix((dat[[i,"par"]])[keep,])}
      if(n.i==0){
        y.tilde.i=rep(0, p) 
        S.tilde.i=array(0, dim=c(p,p))
        phidet.i=NULL
        r.i=NULL
      }
      if (n.i>=1)
      {
        y.tilde.i=colSums(Y.i)
        S.tilde.i= t(Y.i)%*%Y.i
        phidet.i=1 
      }
      y.tilde=y.tilde+y.tilde.i 
      S.tilde=S.tilde+S.tilde.i
      r.all=c(r.all, rep(1,n.i))	
      phidet=c(phidet, phidet.i)
    }
    y.tilde=y.tilde/sum(r.all)
    Y.Y[[z]]=S.tilde
    S.tilde=S.tilde-y.tilde%*%t(y.tilde)*sum(r.all)
    mu.hat[z]=list(y.tilde)
    sigma.hat[z]=list(S.tilde)
    phidet.hat[z]=list(phidet)
    r[z]=list(r.all)
  }
  output=list(mu.hat=mu.hat, sigma.hat=sigma.hat, phidet=phidet.hat,r=r,Y.Y=Y.Y)
  return(output)
}

##################################################################
# LogPredDen: posterior predictive density
##################################################################
LogPredDen=function(data.train, data.pred, d, delta, omega, Parmfit){
  ############ (17) in Section 4.1
  stopifnot(nrow(data.pred[[2]])==length(d))
  data.pred[[1,4]]=as.numeric(d)
  data <- rbind(data.train, data.pred)
  e1=(data.pred[[1,4]]==1); e2=(data.pred[[1,"predzone"]]==2)
  z1=1*e1*e2+2*e1*(1-e2)+3*(1-e1)*e2+4*(1-e1)*(1-e2); z2=3*e1*e2+4*e1*(1-e2)+1*(1-e1)*e2+2*(1-e1)*(1-e2)
  
  r.1=ParmFit$r[[z1]] 
  r.2=ParmFit$r[[z2]] 
  ytilde.z=(ParmFit$mu.hat[[z1]]*length(r.1)+data.pred[[1,2]])/(length(r.1)+1)
  Stilde.z=ParmFit$sigma.hat[[z1]]+(as.numeric(data.pred[[1,2]]))%*%t(as.numeric(data.pred[[1,2]]))+(length(r.1))*ParmFit$mu.hat[[z1]]%*%t(ParmFit$mu.hat[[z1]])-(length(r.1)+1)*as.numeric(ytilde.z)%*%t(as.numeric(ytilde.z))
  output=-p/2*(log(sum(r.2)))-loghconst(delta+length(r.2)-1,  ParmFit$sigma.hat[[z2]]+omega[[z2]])-p/2*(log(sum(r.1)+1))-loghconst(delta+length(r.1)+1-1, Stilde.z+omega[[z1]])
  return(output)
}
