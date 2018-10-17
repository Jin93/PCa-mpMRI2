#include <stdlib.h>
#include <iostream>
#include <string>
//#include <nmath.h>
//#include <bessel.h>
//#include <Rcpp.h>
#include <RcppArmadillo.h> 
// if you include RcppArmadillo.h you shouldn't include Rcpp.h anymore.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppTN.h> 
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>
//#include <Rmath.h>
//#include <R.h>

// [[Rcpp::depends(RcppTN)]]
using namespace arma; //comments
using namespace Rcpp; //comments

// [[Rcpp::export]]
double sum_vec(vec x) {
  double result = 0; /* have to specify the type of the variable first, and types have to match. */
  for(unsigned int i = 0; i < x.n_elem; i++){  /* n_elem: armadillo; '.': means $ */
    result = result + x(i); /* parenthesis not a racket */
  }
  return(result);
}

// [[Rcpp::export]]
double sum_2vec(vec x, vec y) {
  double result = 0; /* have to specify the type of the variable first, and types have to match. */
  for(unsigned int i = 0; i < x.n_elem; i++){  /* n_elem: armadillo; '.': means $ */
    result = result + x(i)*y(i); /* parenthesis not a racket */
  }
  return(result);
}

// [[Rcpp::export]]
double dist2cpp(double a1, double a2, double b1, double b2){
  return(sqrt(pow(a1-b1,2)+pow(a2-b2,2)));
}

// [[Rcpp::export]]
vec ordercpp(vec x, int m){
  int l = x.n_elem;
  vec y(m);
  y.fill(0.5);
  vec a = sort(x);
  vec index(l);
  for(int i = 0; i < l; i++){
    index(i) = i;
  }
  for(int i = 0; i < m; i++)
  {
    for(int j = 0; j < l; j++)
    {
      if((x(j) == a(index(i)))&(y(i) == 0.5)){
        y(i) = j; // order in a C++ vector! not j+1.
        x(j) = 10;
      }
    }
  }
  return(y);
}
  

// iNNIndx: the number of indexes before the starting point of the (i+1)th location
// iNN: # of nearest neighbors for the (i+1)th location

// [[Rcpp::export]]
List getNNIndxcpp(int i, int m){
  int iNNIndx = 0;
  int iNN = 0;
  if(i == 0){
    iNNIndx = 0;//this should never be accessed
    iNN = 0;
  }else if(i < m){
    iNNIndx = i*(i-1)/2;
    iNN = i; // i-1 in R
  }else{
    iNNIndx = m*(m-1)/2+(i-m)*m;
    iNN = m;
  } 
  return List::create(
    _["iNNIndx"] = iNNIndx, // more readable, can put them in one line
    _["iNN"] = iNN
  );
}

// Allocated for the nearest neighbor index vector (note, first location has no neighbors), nnDist=d
// return: nnIndx, nnDist, nnIndxLU

// [[Rcpp::export]]
List mkNNIndxTree0cpp(int n, int m, mat coords){
  int nIndx = (1+m)/2*m+(n-m-1)*m; // calculate nIndx
  vec nnIndx(nIndx);
  vec nnDist(nIndx);
  vec nnIndxLU(2*n);
  nnIndx.fill(0);
  nnIndxLU.fill(0);
  nnDist.fill(0);
  // For (i+1)th location:
  for (int i = 1; i < m; i++)
  {
    int iNNIndx = getNNIndxcpp(i,m)["iNNIndx"]; 
    // the "real" location BEFORE the starting location!! 
    // (not the location in Rcpp vector =-=)
    int iNN = getNNIndxcpp(i,m)["iNN"];
    for (int j = 0; j < i; j++)
    {
      nnDist(iNNIndx+j) = dist2cpp(coords(i,0),coords(i,1),coords(j,0),coords(j,1));
      nnIndx(iNNIndx+j) = j+1;
    }
    nnIndxLU(i) = iNNIndx; 
    nnIndxLU(n+i) = iNN;
  }
  // i>m: nearest neighbor:
  for (int i = m; i < n; i++)
  {
    int iNNIndx = getNNIndxcpp(i,m)["iNNIndx"]; 
    // the "real" location BEFORE the starting location!! 
    // (not the location in Rcpp vector =-=)
    int iNN = getNNIndxcpp(i,m)["iNN"];
    // get the row indicator in coords of the n.neighbors:
    vec dis(i); 
    for (int j = 0; j < i; j++)
    {
      dis(j) = dist2cpp(coords(i,0),coords(i,1),coords(j,0),coords(j,1));
    }
    vec indx = ordercpp(dis,m);
    for (int j = 0; j < m; j++)
    {
      nnDist(m*(m-1)/2+j) = dis(indx(j));
      nnIndx(m*(m-1)/2+j) = indx(j)+1;
    }
    nnIndxLU(i) = iNNIndx; // the location BEFORE the starting location!!
    nnIndxLU(n+i) = iNN;
      // print(i) #40 seconds, can do parallel
  }
  return List::create(
    _["nnIndx"] = nnIndx, // more readable, can put them in one line
    _["nnIndxLU"] = nnIndxLU,
    _["nnDist"] = nnDist
  ); 
}

// Update B and F

// [[Rcpp::export]]
List updateBFcpp(int m, vec d, vec D, int nIndx, vec nnIndxLU, vec CIndx, int n, double sigmaSq, double tausq, double phi, double nu)
{
  vec c(nIndx);
  vec B(nIndx);
  vec FF(n);
  c.fill(0);
  B.fill(0);
  vec C(D.n_elem);
  FF(0) = sigmaSq;
  vec tem;
  double dd = 0;
  double ddd = 0;
  mat Cinv;
  // C_(sij,sij): diag(sigma^2,n)
  // c: c_(sij,N(sij)):
  for(int i = 1; i < n; i++)
  {
    int iNN = getNNIndxcpp(i, m)["iNN"];
    // nnIndxLU[n+i] is always > 0 except when i=1: NA
    for (int k = 0; k < nnIndxLU(n+i); k++)
    {
      dd = d(nnIndxLU(i)+k);
      if (dd == 0) 
      {c(nnIndxLU(i)+k) = sigmaSq;}
      if (dd > 0)
      {c(nnIndxLU(i)+k) = sigmaSq*exp(-(dd*phi));}
      for(int l = 0; l < nnIndxLU(n+i); l++)
      {
        ddd = D(CIndx(i)+l*nnIndxLU(n+i)+k);
        if (ddd == 0)
        {C(CIndx(i)+l*nnIndxLU(n+i)+k) = sigmaSq;}
        if (ddd > 0)
        {C(CIndx(i)+l*nnIndxLU(n+i)+k) = sigmaSq*exp(-(ddd*phi));}
      }
    }
    mat temC(iNN,iNN);
    tem = C.rows(CIndx(i),(CIndx(i)+CIndx(n+i)-1));
    for (int j = 0; j < iNN; j++)
    {
      temC.col(j) = tem.rows(iNN*j,(iNN*(j+1)-1)); //(3*i),(3*(i+1)-1)
    }
      Cinv=inv(temC);
      B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)) = (c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t() * Cinv).t();
      mat tempe = B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t() * c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1));
      FF(i) = sigmaSq - tempe(0,0);
      //B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv
      //FF[i]=sigmaSq + tausq -B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%t(t(c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]))
      for (int j = 0; j < iNN; j++)
      {
        for (int k = 0; k < iNN; k++)
        {
          C(CIndx(i)+j+iNN*k) = Cinv(j,k); 
        }
      }
      //print(i)#:25 seconds
  }
  return List::create(
    _["c"] = c, // more readable, can put them in one line
    _["B"] = B,
    _["C"] = C,
    _["FF"] = FF
  ); 
}

// [[Rcpp::depends(BH)]]


// [[Rcpp::export]]
List updateBFcpp_tau(int m, vec d, vec D, int nIndx, vec nnIndxLU, vec CIndx, int n, double sigmaSq, double tausq, double phi, double nu)
{
  vec c(nIndx);
  vec B(nIndx);
  vec FF(n);
  c.fill(0);
  B.fill(0);
  vec C(D.n_elem);
  FF(0) = sigmaSq + tausq;
  vec tem;
  double dd = 0;
  double ddd = 0;
  mat Cinv;
  // C_(sij,sij): diag(sigma^2,n)
  // c: c_(sij,N(sij)):
  for(int i = 1; i < n; i++)
  {
    int iNN = getNNIndxcpp(i, m)["iNN"];
    // nnIndxLU[n+i] is always > 0 except when i=1: NA
    for (int k = 0; k < nnIndxLU(n+i); k++)
    {
      dd = d(nnIndxLU(i)+k);
      if (dd == 0) 
      {c(nnIndxLU(i)+k) = sigmaSq + tausq;}
      if (dd > 0)
      {c(nnIndxLU(i)+k) = sigmaSq*exp(-(dd*phi));}
      for(int l = 0; l < nnIndxLU(n+i); l++)
      {
        ddd = D(CIndx(i)+l*nnIndxLU(n+i)+k);
        if (ddd == 0)
        {C(CIndx(i)+l*nnIndxLU(n+i)+k) = sigmaSq + tausq;}
        if (ddd > 0)
        {C(CIndx(i)+l*nnIndxLU(n+i)+k) = sigmaSq*exp(-(ddd*phi));}
      }
    }
    mat temC(iNN,iNN);
    tem = C.rows(CIndx(i),(CIndx(i)+CIndx(n+i)-1));
    for (int j = 0; j < iNN; j++)
    {
      temC.col(j) = tem.rows(iNN*j,(iNN*(j+1)-1)); //(3*i),(3*(i+1)-1)
    }
    Cinv=inv(temC);
    B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)) = (c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t() * Cinv).t();
    mat tempe = B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t() * c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1));
    FF(i) = sigmaSq + tausq - tempe(0,0);
    //B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]=c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%Cinv
    //FF[i]=sigmaSq + tausq -B[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]%*%t(t(c[(nnIndxLU[i]+1):(nnIndxLU[i]+nnIndxLU[n+i])]))
    for (int j = 0; j < iNN; j++)
    {
      for (int k = 0; k < iNN; k++)
      {
        C(CIndx(i)+j+iNN*k) = Cinv(j,k); 
      }
    }
    //print(i)#:25 seconds
  }
  return List::create(
    _["c"] = c, // more readable, can put them in one line
    _["B"] = B,
    _["C"] = C,
    _["FF"] = FF
  ); 
}


// [[Rcpp::export]]
List updateBF_materncpp(int m, vec d, vec D, int nIndx, vec nnIndxLU, vec CIndx, int n, double sigmaSq, double tausq, double phi, double nu)
{
  vec c(nIndx);
  vec B(nIndx);
  vec FF(n);
  c.fill(0);
  B.fill(0);
  vec C(D.n_elem);
  FF(0) = sigmaSq;
  vec tem;
  double dd = 0;
  double ddd = 0;
  mat Cinv;
  // C_(sij,sij): diag(sigma^2,n)
  // c: c_(sij,N(sij)):
  for(int i = 1; i < n; i++)
  {
    int iNN = getNNIndxcpp(i, m)["iNN"];
    // nnIndxLU[n+i] is always > 0 except when i=1: NA
    for (int k = 0; k < nnIndxLU(n+i); k++)
    {
      dd = d(nnIndxLU(i)+k);
      if (dd == 0) 
      {c(nnIndxLU(i)+k) = sigmaSq;}
      if (dd > 0)
      {c(nnIndxLU(i)+k) = boost::math::cyl_bessel_k(nu,pow(2*nu,0.5)*dd*phi) * sigmaSq*pow(2,1-nu) *pow(pow(2*nu,0.5)*dd*phi,nu)/R::gammafn(nu);}
      for(int l = 0; l < nnIndxLU(n+i); l++)
      {
        ddd = D(CIndx(i)+l*nnIndxLU(n+i)+k);
        if (ddd == 0)
        {C(CIndx(i)+l*nnIndxLU(n+i)+k) = sigmaSq;}
        if (ddd > 0)
        {C(CIndx(i)+l*nnIndxLU(n+i)+k) = sigmaSq * pow(2,1-nu) *pow(pow(2*nu,0.5)*ddd*phi,nu) * boost::math::cyl_bessel_k(nu,pow(2*nu,0.5)*ddd*phi)/R::gammafn(nu);}
      }
    }
    mat temC(iNN,iNN);
    tem = C.rows(CIndx(i),(CIndx(i)+CIndx(n+i)-1));
    for (int j = 0; j < iNN; j++)
    {
      temC.col(j) = tem.rows(iNN*j,(iNN*(j+1)-1)); //(3*i),(3*(i+1)-1)
    }
    Cinv=inv(temC);
    B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)) = Cinv * c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1));
    mat tempe = B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t() * c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1));
    FF(i) = sigmaSq - tempe(0,0);
    for (int j = 0; j < iNN; j++)
    {
      for (int k = 0; k < iNN; k++)
      {
        C(CIndx(i)+j+iNN*k) = Cinv(j,k); 
      }
    }
    //print(i)#:25 seconds
  }
  return List::create(
    _["c"] = c, // more readable, can put them in one line
    _["B"] = B,
    _["C"] = C,
    _["FF"] = FF
  ); 
}

// [[Rcpp::export]]
List updateBF_materncpp_tau(int m, vec d, vec D, int nIndx, vec nnIndxLU, vec CIndx, int n, double sigmaSq, double tausq, double phi, double nu)
{
  vec c(nIndx);
  vec B(nIndx);
  vec FF(n);
  c.fill(0);
  B.fill(0);
  vec C(D.n_elem);
  FF(0) = sigmaSq + tausq;
  vec tem;
  double dd = 0;
  double ddd = 0;
  mat Cinv;
  // C_(sij,sij): diag(sigma^2,n)
  // c: c_(sij,N(sij)):
  for(int i = 1; i < n; i++)
  {
    int iNN = getNNIndxcpp(i, m)["iNN"];
    // nnIndxLU[n+i] is always > 0 except when i=1: NA
    for (int k = 0; k < nnIndxLU(n+i); k++)
    {
      dd = d(nnIndxLU(i)+k);
      if (dd == 0) 
      {c(nnIndxLU(i)+k) = sigmaSq + tausq;}
      if (dd > 0)
      {c(nnIndxLU(i)+k) = boost::math::cyl_bessel_k(nu,pow(2*nu,0.5)*dd*phi) * sigmaSq*pow(2,1-nu) *pow(pow(2*nu,0.5)*dd*phi,nu)/R::gammafn(nu);}
      for(int l = 0; l < nnIndxLU(n+i); l++)
      {
        ddd = D(CIndx(i)+l*nnIndxLU(n+i)+k);
        if (ddd == 0)
        {C(CIndx(i)+l*nnIndxLU(n+i)+k) = sigmaSq + tausq;}
        if (ddd > 0)
        {C(CIndx(i)+l*nnIndxLU(n+i)+k) = sigmaSq * pow(2,1-nu) *pow(pow(2*nu,0.5)*ddd*phi,nu) * boost::math::cyl_bessel_k(nu,pow(2*nu,0.5)*ddd*phi)/R::gammafn(nu);}
      }
    }
    mat temC(iNN,iNN);
    tem = C.rows(CIndx(i),(CIndx(i)+CIndx(n+i)-1));
    for (int j = 0; j < iNN; j++)
    {
      temC.col(j) = tem.rows(iNN*j,(iNN*(j+1)-1)); //(3*i),(3*(i+1)-1)
    }
    Cinv=inv(temC);
    B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)) = Cinv * c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1));
    mat tempe = B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t() * c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1));
    FF(i) = sigmaSq + tausq - tempe(0,0);
    for (int j = 0; j < iNN; j++)
    {
      for (int k = 0; k < iNN; k++)
      {
        C(CIndx(i)+j+iNN*k) = Cinv(j,k); 
      }
    }
    //print(i)#:25 seconds
  }
  return List::create(
    _["c"] = c, // more readable, can put them in one line
    _["B"] = B,
    _["C"] = C,
    _["FF"] = FF
  ); 
}


//  M-H for h:
// [[Rcpp::export]]
vec hfcpp(int n,vec w,vec c,vec C,vec FF,uvec nnIndx,vec nnIndxLU,vec CIndx,double sigmaSq)
{
  vec output(n);
  // 1st pixel:
  output(0) = R::dnorm(w(0),0,sqrt(sigmaSq),0);
  for (int i = 1; i < n; i++){
    int iNN = sqrt(CIndx(n+i));
    mat tem1(nnIndxLU(n+i),nnIndxLU(n+i));
    for (int k = 0; k < nnIndxLU(n+i); k++)\
    {
      tem1.col(k) = C.rows(CIndx(i)+k*iNN,(CIndx(i)+(k+1)*iNN-1));
    }
    uvec idx = nnIndx.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)) - 1;
    mat tempmean = c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t() * tem1 * w.rows(idx);
    output(i) = R::dnorm(w(i),tempmean(0,0), sqrt(FF(i)), 0);
  }
  return(output);
}

// [[Rcpp::export]]
vec hfnngpcpp(int n,vec w,vec c,vec C,vec FF,uvec nnIndx,vec nnIndxLU,vec CIndx,double sigmaSq)
{
  vec output(n);
  // 1st pixel:
  output(0) = R::dnorm(w(0),0,sqrt(sigmaSq),0);
  for (int i = 1; i < n; i++){
    int iNN = sqrt(CIndx(n+i));
    mat tem1(nnIndxLU(n+i),nnIndxLU(n+i));
    for (int k = 0; k < nnIndxLU(n+i); k++)
    {
      tem1.col(k) = C.rows(CIndx(i)+k*iNN,(CIndx(i)+(k+1)*iNN-1));
    }
    uvec idx = nnIndx.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)) - 1;
    mat tempmean = c.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t() * tem1 * w.rows(idx);
    output(i) = R::dnorm(w(i),tempmean(0,0), sqrt(FF(i)), 0);
  }
  return(output);
}

// [[Rcpp::export]]
vec hfnngp2cpp(int n,vec w,vec B,vec FF,uvec nnIndx,vec nnIndxLU,vec CIndx,double sigmaSq)
{
  vec output(n);
  // 1st pixel:
  output(0) = R::dnorm(w(0),0,sqrt(sigmaSq),0);
  for (int i = 1; i < n; i++){
    uvec idx = nnIndx.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)) - 1;
    mat tempmean = B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t() * w.rows(idx);
    output(i) = R::dnorm(w(i),tempmean(0,0), sqrt(FF(i)), 0);
  }
  return(output);
}


// [[Rcpp::export()]]
double get1TN(double mean, double sd, double low, double high) {
  double draw = 0.0 ;
  bool valid = false ;
  while (!valid) {
    NumericVector cand1 = rnorm(1, mean, sd) ;
    if ((cand1(0) <= high) & (cand1(0) >= low)) {
      draw = cand1(0) ;
      valid = true ;
    }
  }
  return(draw) ;
}

// [[Rcpp::export]]
List qiupdatecpp(double sigmaSq, double tausq, vec lambda, vec FF, vec cancer, vec q, vec d, int n, vec B,
              vec nIndx, uvec nnIndx, vec nnIndxLU, vec uIndx, vec uIndxLU, vec uiIndx)
{
  vec m(n),M(n),G(n),V(n),O(n),pstar(n);
  // 1st element:
  int ii=0;
  m(ii)=0;
  G(ii)=1/(sigmaSq+tausq);
  for (int jj = 0; jj < uIndxLU(n+ii); jj++)
  {
    int t = uIndx(uIndxLU(ii)+jj) - 1;
    int l = uiIndx(uIndxLU(ii)+jj) - 1;
    uvec indx = nnIndx.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)) - 1;
    double a = q(t)-sum_2vec(B.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)), q.rows(indx));
    //double indxx = indx(l);
    a = a + B(nnIndxLU(t)+l) * q(indx(l));
    m(ii) = m(ii) + B(nnIndxLU(t)+l)*a/FF(t);
    G(ii) = G(ii) + pow(B(nnIndxLU(t)+l),2)/FF(t);
  }
  G(ii) = 1/G(ii);
  M(ii) = m(ii) * G(ii);
  if (cancer(ii) == 1)
  {
    // Truncated Normal:
    // q(ii) = get1TN(M(ii),sqrt(G(ii)),-lambda(ii),INFINITY);
    q(ii) = RcppTN::rtn1(M(ii),sqrt(G(ii)),-lambda(ii),INFINITY);
  }
  if (cancer[ii]==0)
  {
    // q(ii) = get1TN(M(ii),sqrt(G(ii)),-INFINITY,-lambda(ii));
    q(ii) = RcppTN::rtn1(M(ii),sqrt(G(ii)),-INFINITY,-lambda(ii)); 
  }
  pstar(ii) = 1-R::pnorm(-lambda(ii),M(ii),sqrt(G(ii)),1,0);
  // 2~nth element: Cholesky decomposition:
    for (int ii = 1; ii < n; ii++)
    {
      if (uIndxLU(n+ii) == 0)
      {
        uvec idx2 = nnIndx.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)) - 1;
        // mat aa = B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)).t() * q.rows(idx2);
        // m(ii) =  aa(0,0)/FF(ii);
        double aa =  sum_2vec(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)), q.rows(idx2));
        m(ii) =  aa/FF(ii);
        G(ii) = FF(ii);
        M(ii) = m(ii) * G(ii);
        if (cancer(ii) == 1)
        { 
          // q(ii) = get1TN(M(ii), sqrt(G(ii)), -lambda(ii), INFINITY);
          q(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), -lambda(ii), INFINITY);
        }
        if (cancer[ii]==0)
        {
          // q(ii) = get1TN(M(ii), sqrt(G(ii)), -INFINITY, -lambda(ii));
          q(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), -INFINITY, -lambda(ii));
        }
        pstar(ii)=1-R::pnorm(-lambda(ii), M(ii), sqrt(G(ii)), 1, 0);
      }
      if (uIndxLU(n+ii) >0)
      {
        uvec idx2 = nnIndx.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)) - 1;
        // mat aa = B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)).t() * q.rows(idx2);
        // m(ii) =  aa(0,0)/FF(ii);
        double aa =  sum_2vec(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)), q.rows(idx2));
        m(ii) =  aa/FF(ii);
        G(ii) = 1/FF(ii);
        for (int jj = 0; jj < uIndxLU(n+ii); jj++)
        {
          int t = uIndx(uIndxLU(ii)+jj) - 1;
          int l = uiIndx(uIndxLU(ii)+jj) - 1;
          uvec indx = nnIndx.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)) - 1;
          double temmat = sum_2vec(B.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)), q.rows(indx));
          double a = q(t)-temmat;
          // double indxx = indx(l);
          a = a + B(nnIndxLU(t)+l) * q(indx(l));
          m(ii) = m(ii) + B(nnIndxLU(t)+l)*a/FF(t);
          G(ii) = G(ii) + pow(B(nnIndxLU(t)+l),2)/FF(t);
        }
        G(ii) = 1/G(ii);
        M(ii) = m(ii) * G(ii);
        if (cancer(ii)==1)
        {
          // q(ii) = get1TN(M(ii), sqrt(G(ii)), -lambda(ii), INFINITY);
          q(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), -lambda(ii), INFINITY);
        }
        if (cancer(ii)==0)
        {
          // q(ii) = get1TN(M(ii), sqrt(G(ii)), -INFINITY, -lambda(ii));
          q(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), -INFINITY, -lambda(ii));
        }
        pstar(ii) = 1-R::pnorm(-lambda[ii], M(ii), sqrt(G(ii)), 1, 0);
      }
    }
    return List::create(_["q"] = q,_["pstar"] = pstar,_["M"] = M,_["G"] = G);
}


// [[Rcpp::export]]
List qiupdatecpp_new(double sigmaSq, double tausq, vec lambda, vec FF, vec cancer, vec q, vec d, int n, vec B,
                 vec nIndx, uvec nnIndx, vec nnIndxLU, vec uIndx, vec uIndxLU, vec uiIndx)
{
  vec m(n),M(n),G(n),V(n),O(n),pstar(n);
  // 1st element:
  int ii=0;
  m(ii)=0;
  G(ii)=1/sigmaSq;
  for (int jj = 0; jj < uIndxLU(n+ii); jj++)
  {
    int t = uIndx(uIndxLU(ii)+jj) - 1;
    int l = uiIndx(uIndxLU(ii)+jj) - 1;
    uvec indx = nnIndx.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)) - 1;
    double a = q(t)-sum_2vec(B.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)), q.rows(indx));
    //double indxx = indx(l);
    a = a + B(nnIndxLU(t)+l) * q(indx(l));
    m(ii) = m(ii) + B(nnIndxLU(t)+l)*a/FF(t);
    G(ii) = G(ii) + pow(B(nnIndxLU(t)+l),2)/FF(t);
  }
  G(ii) = 1/G(ii);
  M(ii) = m(ii) * G(ii);
  G(ii) = G(ii) + tausq;
  if (cancer(ii) == 1)
  {
    // Truncated Normal:
    // q(ii) = get1TN(M(ii),sqrt(G(ii)),-lambda(ii),INFINITY);
    q(ii) = RcppTN::rtn1(M(ii),sqrt(G(ii)),-lambda(ii),INFINITY);
  }
  if (cancer[ii]==0)
  {
    // q(ii) = get1TN(M(ii),sqrt(G(ii)),-INFINITY,-lambda(ii));
    q(ii) = RcppTN::rtn1(M(ii),sqrt(G(ii)),-INFINITY,-lambda(ii)); 
  }
  pstar(ii) = 1-R::pnorm(-lambda(ii),M(ii),sqrt(G(ii)),1,0);
  // 2~nth element: Cholesky decomposition:
  for (int ii = 1; ii < n; ii++)
  {
    if (uIndxLU(n+ii) == 0)
    {
      uvec idx2 = nnIndx.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)) - 1;
      // mat aa = B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)).t() * q.rows(idx2);
      // m(ii) =  aa(0,0)/FF(ii);
      double aa =  sum_2vec(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)), q.rows(idx2));
      m(ii) =  aa/FF(ii);
      G(ii) = FF(ii);
      M(ii) = m(ii) * G(ii);
      G(ii) = G(ii) + tausq;
      if (cancer(ii) == 1)
      { 
        // q(ii) = get1TN(M(ii), sqrt(G(ii)), -lambda(ii), INFINITY);
        q(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), -lambda(ii), INFINITY);
      }
      if (cancer[ii]==0)
      {
        // q(ii) = get1TN(M(ii), sqrt(G(ii)), -INFINITY, -lambda(ii));
        q(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), -INFINITY, -lambda(ii));
      }
      pstar(ii)=1-R::pnorm(-lambda(ii), M(ii), sqrt(G(ii)), 1, 0);
    }
    if (uIndxLU(n+ii) >0)
    {
      uvec idx2 = nnIndx.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)) - 1;
      // mat aa = B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)).t() * q.rows(idx2);
      // m(ii) =  aa(0,0)/FF(ii);
      double aa =  sum_2vec(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)), q.rows(idx2));
      m(ii) =  aa/FF(ii);
      G(ii) = 1/FF(ii);
      for (int jj = 0; jj < uIndxLU(n+ii); jj++)
      {
        int t = uIndx(uIndxLU(ii)+jj) - 1;
        int l = uiIndx(uIndxLU(ii)+jj) - 1;
        uvec indx = nnIndx.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)) - 1;
        double temmat = sum_2vec(B.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)), q.rows(indx));
        double a = q(t)-temmat;
        // double indxx = indx(l);
        a = a + B(nnIndxLU(t)+l) * q(indx(l));
        m(ii) = m(ii) + B(nnIndxLU(t)+l)*a/FF(t);
        G(ii) = G(ii) + pow(B(nnIndxLU(t)+l),2)/FF(t);
      }
      G(ii) = 1/G(ii);
      M(ii) = m(ii) * G(ii);
      G(ii) = G(ii) + tausq;
      if (cancer(ii)==1)
      {
        // q(ii) = get1TN(M(ii), sqrt(G(ii)), -lambda(ii), INFINITY);
        q(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), -lambda(ii), INFINITY);
      }
      if (cancer(ii)==0)
      {
        // q(ii) = get1TN(M(ii), sqrt(G(ii)), -INFINITY, -lambda(ii));
        q(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), -INFINITY, -lambda(ii));
      }
      pstar(ii) = 1-R::pnorm(-lambda[ii], M(ii), sqrt(G(ii)), 1, 0);
    }
  }
  return List::create(_["q"] = q,_["pstar"] = pstar,_["M"] = M,_["G"] = G);
}

// [[Rcpp::export]]
List nngpwupdatecpp(double sigmaSq, double tausq, vec lambda, vec FF, vec cancer, vec q, vec d, int n, vec B,
                     vec nIndx, uvec nnIndx, vec nnIndxLU, vec uIndx, vec uIndxLU, vec uiIndx, vec kij)
{
  vec m(n),M(n),G(n),V(n),O(n),pstar(n),postmean(n),postvar(n);
  // 1st element:
  int ii=0;
  m(ii)=0;
  G(ii)=1/sigmaSq;
  for (int jj = 0; jj < uIndxLU(n+ii); jj++)
  {
    int t = uIndx(uIndxLU(ii)+jj) - 1;
    int l = uiIndx(uIndxLU(ii)+jj) - 1;
    uvec indx = nnIndx.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)) - 1;
    double a = q(t)-sum_2vec(B.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)), q.rows(indx));
    //double indxx = indx(l);
    a = a + B(nnIndxLU(t)+l) * q(indx(l));
    m(ii) = m(ii) + B(nnIndxLU(t)+l)*a/FF(t);
    G(ii) = G(ii) + pow(B(nnIndxLU(t)+l),2)/FF(t);
  }
  G(ii) = 1/G(ii);
  M(ii) = m(ii) * G(ii);
  //G(ii) = G(ii) + tausq;
  postvar(ii) = 1/(1/tausq+1/G(ii));
  postmean(ii) = postvar(ii) * ((kij(ii)-lambda(ii))/tausq + M(ii)/G(ii));
  mat rr = rnorm(1, postmean(ii), sqrt(postvar(ii)));
  q(ii) = rr(0);
  // 2~nth element: Cholesky decomposition:
  for (int ii = 1; ii < n; ii++)
  {
    if (uIndxLU(n+ii) == 0)
    {
      uvec idx2 = nnIndx.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)) - 1;
      double aa =  sum_2vec(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)), q.rows(idx2));
      m(ii) =  aa/FF(ii);
      G(ii) = FF(ii);
      M(ii) = m(ii) * G(ii);
      //G(ii) = G(ii) + tausq;
      postvar(ii) = 1/(1/tausq+1/G(ii));
      postmean(ii) = postvar(ii) * ((kij(ii)-lambda(ii))/tausq + M(ii)/G(ii));
      mat rr = rnorm(1, postmean(ii), sqrt(postvar(ii)));
      q(ii) = rr(0);
    }
    if (uIndxLU(n+ii) >0)
    {
      uvec idx2 = nnIndx.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)) - 1;
      double aa =  sum_2vec(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)), q.rows(idx2));
      m(ii) =  aa/FF(ii);
      G(ii) = 1/FF(ii);
      for (int jj = 0; jj < uIndxLU(n+ii); jj++)
      {
        int t = uIndx(uIndxLU(ii)+jj) - 1;
        int l = uiIndx(uIndxLU(ii)+jj) - 1;
        uvec indx = nnIndx.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)) - 1;
        double temmat = sum_2vec(B.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)), q.rows(indx));
        double a = q(t)-temmat;
        a = a + B(nnIndxLU(t)+l) * q(indx(l));
        m(ii) = m(ii) + B(nnIndxLU(t)+l)*a/FF(t);
        G(ii) = G(ii) + pow(B(nnIndxLU(t)+l),2)/FF(t);
      }
      G(ii) = 1/G(ii);
      M(ii) = m(ii) * G(ii);
      //G(ii) = G(ii) + tausq;
      postvar(ii) = 1/(1/tausq+1/G(ii));
      postmean(ii) = postvar(ii) * ((kij(ii)-lambda(ii))/tausq + M(ii)/G(ii));
      mat rr = rnorm(1, postmean(ii), sqrt(postvar(ii)));
      q(ii) = rr(0);
      }
  }
  return List::create(_["wij"] = q,_["postmean"] = postmean,_["postvar"] = postvar,_["M"] = M,_["G"] = G);
}


// [[Rcpp::export]]
List nngpkijupdatecpp(int n, vec lambda, vec wij, vec cancer, double tausq)
{
  vec kij(n);
  vec pstar(n);
  for (int ii = 0; ii < n; ii++)
  {
    if (cancer(ii)==1)
    {
      kij(ii) = RcppTN::rtn1(lambda(ii) + wij(ii), sqrt(tausq), 0, INFINITY);
    }
    if (cancer(ii)==0)
    {
      kij(ii) = RcppTN::rtn1(lambda(ii) + wij(ii), sqrt(tausq), -INFINITY, 0);
    }
    pstar(ii) = 1-R::pnorm(0, lambda(ii) + wij(ii), sqrt(tausq), 1, 0);
  }
  return List::create(_["kij"] = kij,_["pstar"] = pstar);
}


const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
vec Mahalanobis(mat x, rowvec center, mat cov) {
  int n = x.n_rows;
  mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

// [[Rcpp::export]]
vec dmvnormcpp(mat x, rowvec mean, mat sigma) { 
  vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(log(eig_sym(sigma)));
  vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  return(exp(logretval));
}

// [[Rcpp::export]]
vec logdmvnormcpp(mat x, rowvec mean, mat sigma) { 
  vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(log(eig_sym(sigma)));
  vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  return(logretval);
}

// ###### Cancer full conditional:

// [[Rcpp::export]]
double piupdatecpp(int zon, vec y, List mu, vec Delta, List cov, double pstar)
{
  mat yy(1,y.n_elem);
  yy.row(0) = y.t();
  int indicat1 = 1 * (zon == 1) + 2 * (zon == 0);
  int indicat0 = 3 * (zon == 1) + 4 * (zon == 0);
  vec mu1 = mu[indicat1-1];
  mat cov1 = cov[indicat1-1];
  vec mu0 = mu[indicat0-1];
  mat cov0 = cov[indicat0-1];
  rowvec m1(Delta.n_elem);
  rowvec m0(Delta.n_elem);
  m1.row(0) = (mu1+Delta).t();
  m0.row(0) = (mu0+Delta).t();
  vec tem1 = dmvnormcpp(yy,m1,cov1);
  vec tem0 = dmvnormcpp(yy,m0,cov0);
  double p1 = tem1(0) * pstar;
  double p0 = tem0(0) * (1-pstar);
  double ppcancer=p1/(p0+p1);
  return (ppcancer);
}

// [[Rcpp::export]]
double piupdatecpp_nodelta(int zon, vec y, List mu, List cov, double pstar)
{
  mat yy(1,y.n_elem);
  yy.row(0) = y.t();
  int indicat1 = 1 * (zon == 1) + 2 * (zon == 0);
  int indicat0 = 3 * (zon == 1) + 4 * (zon == 0);
  vec mu1 = mu[indicat1-1];
  mat cov1 = cov[indicat1-1];
  vec mu0 = mu[indicat0-1];
  mat cov0 = cov[indicat0-1];
  rowvec m1(mu1.n_elem);
  rowvec m0(mu1.n_elem);
  m1.row(0) = (mu1).t();
  m0.row(0) = (mu0).t();
  vec tem1 = dmvnormcpp(yy,m1,cov1);
  vec tem0 = dmvnormcpp(yy,m0,cov0);
  double p1 = tem1(0) * pstar;
  double p0 = tem0(0) * (1-pstar);
  double ppcancer=p1/(p0+p1);
  return (ppcancer);
}

// ###### Simulation Data Generating:

// A function that can return distance matrix?

// // [[Rcpp::export]]
//mat simulatemri(mat coords)
//{
//  mat D = R::dist(coords);
//  return(D);
//}

// rmvn function:

// [[Rcpp::export]]
mat rmvncpp(int n, vec mu, mat V){
  int p = mu.n_elem;
  //if(any(is.na(match(dim(V),p))))
  //stop("Dimension problem!")
    mat D = chol(V);
    mat T1(n,p);
    mat T2(n,p);
    for(int co = 0; co < p; co ++)
    {
      T1.col(co) = (vec)rnorm(n,0,1);
      vec tem(n);
      tem.fill(mu(co));
      T2.col(co) = tem;
    }
    mat T = (T1 * D + T2).t();
    return(T);
}

// [[Rcpp::export]]
mat cholcpp(mat V){
  mat D = chol(V);
  return(D);
}

// [[Rcpp::export]]
mat solvecpp(mat V){
  mat D = inv(V);
  return(D);
}

// [[Rcpp::export]]
mat multiplycpp(mat a, mat b){
  mat D = a * b;
  return(D);
}

// [[Rcpp::export]]
double rbinomcpp(int zon, vec y, List mu, vec Delta, List cov, double pstar){
  double prob = piupdatecpp(zon,y,mu,Delta,cov,pstar);
  double a = R::rbinom(1,prob);
  return(a);
}

// [[Rcpp::export]]
double rbinomcpp_nodelta(int zon, vec y, List mu, List cov, double pstar){
  double prob = piupdatecpp_nodelta(zon,y,mu,cov,pstar);
  double a = R::rbinom(1,prob);
  return(a);
}


// For loops:

// [[Rcpp::export]]
double hflrcpp(int n,vec w,mat ci,mat invCistar,double tausq, mat diagmat)
{
  mat ww(1,n);
  ww.row(0) = w.t();
  vec mu(n);
  mu.fill(0);
  rowvec Mu(n);
  Mu.row(0) = mu.t();
  mat Cov = tausq * diagmat + ci.t() * invCistar * ci;
  vec lik = dmvnormcpp(ww,Mu,Cov);
  return(lik(0));
}


// update mecoef and va:

// [[Rcpp::export]]
List update2cpp(int ii, int n, mat ci, mat Cistar, mat covmat, mat invCistar, double tausq, mat diagmat)
{
  mat cc = covmat;
  cc.shed_col(ii);
  mat cijtilde=cc.row(ii);
  mat ci_j = ci;
  ci_j.shed_col(ii);
  double cijj = covmat(ii,ii);
  
  mat tempinv = solvecpp(Cistar+ci_j * ci_j.t()/tausq);
  mat tempbracket = (diagmat-ci_j.t() * tempinv * ci_j/pow(tausq,2) );
  mat me = cijtilde * tempbracket;
  mat temp = cijtilde * tempbracket * cijtilde.t();
  double va = cijj - temp(0);
  return List::create(_["me"] = me, _["va"] = va);
}

// [[Rcpp::export]]
List updatemecoefcpp(int x,int n,mat ci,mat Cistar,mat covmat,mat invCistar,double tausq, mat diagn_1)
{
  //diagn_1=diag(1/tausq,(n-1))
  mat mecoef(n,(n-1));
  vec va(n);
  for (int ii = 0; ii < n; ii ++)
  {
    mat cc = covmat;
    cc.shed_col(ii);
    mat cijtilde=cc.row(ii);
    mat ci_j = ci;
    ci_j.shed_col(ii);
    double cijj = covmat(ii,ii);
    mat tempinv = inv(Cistar+ci_j*ci_j.t()/tausq);
    mat tempbracket = (diagn_1-pow(tausq,-2)*ci_j.t() * tempinv * ci_j);
    mecoef.row(ii) = cijtilde * tempbracket;
    mat temp = cijtilde * tempbracket * cijtilde.t();
    va(ii) = cijj - temp(0);
  }
  return List::create(_["mecoef"] = mecoef, _["va"] = va);
}


// [[Rcpp::export]]
List updatelrmscpp(int x,int n,mat ci, mat Cistar, mat covmat, mat invCistar, double tausq)
{
  mat mecoef(n,(n-1));
  vec va(n);
  vec vinvtausq(n-1);
  vinvtausq.fill(1/tausq);
  mat diagm = diagmat(vinvtausq);
  for (int ii = 0; ii < n; ii ++)
    {
      mat cijtilde = covmat.row(ii);
      cijtilde.shed_col(ii);
      mat ci_j = ci;
      ci_j.shed_col(ii);
      double cijj = covmat(ii,ii);
      mat tempinv = inv(Cistar+ci_j * (ci_j).t()/tausq);
      mat tempbracket =  diagm- (ci_j).t() * tempinv * ci_j/pow(tausq,2);
      mat temp2 = cijtilde * tempbracket;
      mecoef.row(ii) = temp2;
      mat tt = temp2 * cijtilde.t();
      va(ii)=cijj - tt(0);
    }
  return List::create(_["mecoef"] = mecoef, _["va"] = va);
}


// [[Rcpp::export]]
mat creatediag(vec v)
{
  mat m1 = diagmat(v); // to test if the diagmat can turn a vector into the diagonal entries of a diagonal matrix.
  return(m1);
}


// [[Rcpp::export]]
List qicarcpp(double sigmasq, vec lambda,vec cancer,vec q, mat Dinv)
{ 
  vec pstar(cancer.n_elem);
  vec M(cancer.n_elem);
  vec G(cancer.n_elem);
  for (int x = 0; x < cancer.n_elem; x ++)
  { 
    mat Qmean = Dinv.row(x) * q - Dinv(x,x) * q(x);
    double qmean = Qmean(0);
    if (cancer(x)==1)
    { 
      q(x)=RcppTN::rtn1(qmean,sqrt(sigmasq),-lambda(x),INFINITY);
    }
    if (cancer[x]==0)
    { 
      q(x)=RcppTN::rtn1(qmean,sqrt(sigmasq),-INFINITY,-lambda(x));
    }
    pstar(x)=1-R::pnorm(-lambda(x), qmean, sqrt(sigmasq), 1, 0);
    M(x) = qmean;
    G(x) = sigmasq;
  }
  return List::create(_["q"] = q, _["pstar"] = pstar, _["M"] = M, _["G"] = G);
}

// [[Rcpp::export]]
List qicarnewcpp(double sigmasq, vec lambda, vec q, mat Dinv)
{ 
  vec pstar(q.n_elem);
  vec M(q.n_elem);
  vec G(q.n_elem);
  for (int x = 0; x < q.n_elem; x ++)
  { 
    mat Qmean = Dinv.row(x) * q - Dinv(x,x) * q(x);
    double qmean = Qmean(0);
    mat rr = rnorm(1, qmean, sqrt(sigmasq));
    q(x) = rr(0);
    M(x) = qmean;
    G(x) = sigmasq;
  }
  return List::create(_["q"] = q, _["M"] = M, _["G"] = G);
}

// [[Rcpp::export]]
List qiupdate_lrcpp(double sigmasq, double tausq, vec q0, vec cancer, vec w, int n, mat mecoef, vec va)
{
  vec pstar(n);
  vec M(n);
  for (int ii = 0; ii < n; ii ++)
  {
    vec wtemp = w;
    wtemp.shed_row(ii);
    mat Me = mecoef.row(ii) * wtemp;
    double me = Me(0);
    if (cancer(ii)==1)
    {
      w(ii) = RcppTN::rtn1(me,sqrt(va(ii)),-q0(ii),INFINITY);
    }
    if (cancer(ii)==0)
    {
      w(ii) = RcppTN::rtn1(me,sqrt(va(ii)),-INFINITY,-q0(ii));
    }
    pstar(ii) = 1-R::pnorm(-q0(ii), me, sqrt(va(ii)), 1, 0);
    M(ii) = me;
  }
  return List::create(_["w"] = w, _["pstar"] = pstar, _["M"] = M, _["G"] = va);
}

// [[Rcpp::export]]
List qiupdate_fullcpp(vec q0,vec cancer, vec w, int n, mat meanfull, vec varfull)
{
  vec pstar(n);
  vec M(n);
  for (int ii = 0; ii < n; ii ++)
  {
    vec wtemp = w;
    wtemp.shed_row(ii);
    mat Me = meanfull.row(ii) * wtemp;
    M(ii) = Me(0);
    if (cancer(ii)==1)
    {
      w(ii)=RcppTN::rtn1(M(ii),sqrt(varfull(ii)),-q0(ii),INFINITY);
    }
    if (cancer(ii)==0)
    {
      w(ii)=RcppTN::rtn1(M(ii),sqrt(varfull(ii)),-INFINITY,-q0(ii));
    }
    pstar(ii)=1-R::pnorm(-q0(ii),M(ii),sqrt(varfull(ii)), 1, 0);
  }
  return List::create(_["w"] = w, _["pstar"] = pstar, _["M"] = M, _["G"] = varfull);
}

//notice: cannot use "." in the name of the variables of a function.

// [[Rcpp::export]]
double calstarcpp(int n, vec nnIndx, vec nnIndxLU, mat Bupdate, vec w, double sigmasq, double tausq, vec FFupdate)
{
  double star = pow(w(0),2)/(sigmasq+tausq);
  for (int i = 1; i < n; i ++)
  {
    double e = 0;
    for(int j = 0; j < nnIndxLU(n+i); j ++)
    {
      e = e + Bupdate(nnIndxLU(i)+j,0) * w(nnIndx(nnIndxLU(i)+j)-1);
    }
    double b = w(i) - e;
    star = star + pow(b,2)/FFupdate(i);
  }
  return(star);
}


// [[Rcpp::export]]
double calstarcpp_new(int n, vec nnIndx, vec nnIndxLU, mat Bupdate, vec w, double sigmasq, vec FFupdate)
{
  double star = pow(w(0),2)/(sigmasq);
  for (int i = 1; i < n; i ++)
  {
    double e = 0;
    for(int j = 0; j < nnIndxLU(n+i); j ++)
    {
      e = e + Bupdate(nnIndxLU(i)+j,0) * w(nnIndx(nnIndxLU(i)+j)-1);
    }
    double b = w(i) - e;
    star = star + pow(b,2)/FFupdate(i);
  }
  return(star);
}


// [[Rcpp::export]]
mat materncovcpp(mat D, double sigmasq, double phi, double nu)
{ 
  mat COV = mat(D.n_rows,D.n_cols);
  for (int i = 0; i < D.n_rows; i ++)
  { 
    for (int j = 0; j < D.n_cols; j ++)
    { 
      if (D(i,j) == 0) {COV(i,j) = sigmasq;}
      if (D(i,j) > 0) {
        COV(i,j)=sigmasq * pow(2,1-nu) *pow(pow(2*nu,0.5)*D(i,j)*phi,nu) * boost::math::cyl_bessel_k(nu,pow(2*nu,0.5)*D(i,j)*phi)/R::gammafn(nu);
        //COV(i,j) = boost::math::cyl_bessel_k(nu,pow(2*nu,0.5)*D(i,j)*phi) * sigmasq*pow(2,1-nu) *pow(pow(2*nu,0.5)*D(i,j)*phi,nu)/R::gammafn(nu);
      }
    }
  }
  return(COV);
}

// [[Rcpp::export]]
double maternvarcpp(double d, double sigmasq, double phi, double nu)
{ 
  double varout;
      if (d == 0) {varout = sigmasq;}
      if (d > 0) {
        varout=sigmasq * pow(2,1-nu) *pow(pow(2*nu,0.5)*d*phi,nu) * boost::math::cyl_bessel_k(nu,pow(2*nu,0.5)*d*phi)/R::gammafn(nu);
        //COV(i,j) = boost::math::cyl_bessel_k(nu,pow(2*nu,0.5)*D(i,j)*phi) * sigmasq*pow(2,1-nu) *pow(pow(2*nu,0.5)*D(i,j)*phi,nu)/R::gammafn(nu);
      }
  return(varout);
}

// [[Rcpp::export]]
List besselkcpp(double d, double sigmasq, double phi, double nu)
{ 
  double bvalue1 = boost::math::cyl_bessel_k(nu,pow(2*nu,0.5)*d*phi);
  double bvalue2 = boost::math::cyl_bessel_k(pow(2*nu,0.5)*d*phi,nu);
  double gammavalue = R::gammafn(nu);
  return List::create(_["bvalue1"] = bvalue1, _["bvalue2"] = bvalue2, _["gammavalue"] = gammavalue);
}

// [[Rcpp::export]]
List updatecicpp(double sigmasq, double tausq, double phi, double nu, mat Diknots, mat Direst, mat nuggetcov)
{
  mat Cistar = materncovcpp(Diknots, sigmasq, phi, nu);
  mat invCistar = solvecpp(Cistar);
  mat ci = materncovcpp(Direst, sigmasq, phi, nu);
  // nuggetcov=diag(tausq,nrow(fillnew[[x]]))
  mat covmat = ci.t() * invCistar * ci + nuggetcov;
  return List::create(_["ci"] = ci, _["Cistar"] = Cistar, _["covmat"] = covmat, _["invCistar"] = invCistar);
}

// [[Rcpp::export]]
List updatecietc(int N, vec n, double sigmasq,double tausq, double phi, double nu, List Diknots, List Direst)
{
  Rcpp::List out(N);
  for (int i = 0; i < N; i ++) 
  {
    mat diagmatrix = mat(n(i),n(i));
    diagmatrix.eye();
    out[i]=updatecicpp(sigmasq,tausq,phi,nu,Diknots[i],Direst[i],diagmatrix*tausq);
  }
  return(out);
}



// New LR update:

// [[Rcpp::export]]
List updatecicpp_new(double sigmasq, double tausq, double phi, double nu, mat Diknots, mat Direst)
{
  mat Cistar = materncovcpp(Diknots, sigmasq, phi, nu);
  mat invCistar = solvecpp(Cistar);
  mat ci = materncovcpp(Direst, sigmasq, phi, nu);
  // nuggetcov=diag(tausq,nrow(fillnew[[x]]))
  mat covmat = ci.t() * invCistar * ci;
  return List::create(_["ci"] = ci, _["Cistar"] = Cistar, _["covmat"] = covmat, _["invCistar"] = invCistar);
}

// [[Rcpp::export]]
List updatelrkcpp(double sigmaSq, double tausq, vec lambda, vec cancer, vec wstar, int n, mat ci, mat ciinvCistar)
{
  vec M(n),G(n),k(n),pstar(n);
  // for ii-th element: 
  for (int ii = 0; ii < n; ii++)
  {
    mat mtem = ciinvCistar.row(ii) * wstar;
    M(ii) =  lambda(ii) + mtem(0);
    G(ii) = tausq;
    if (cancer(ii)==1)
    {
      k(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), 0, INFINITY);
    }
    if (cancer(ii)==0)
    {
      k(ii) = RcppTN::rtn1(M(ii), sqrt(G(ii)), -INFINITY, 0);
    }
    pstar(ii) = 1-R::pnorm(0, M(ii), sqrt(G(ii)), 1, 0); 
  }
  return List::create(_["k"] = k,_["pstar"] = pstar,_["M"] = M,_["G"] = G);
}

// [[Rcpp::export]]
mat updatelrwstarcpp(double tausq, vec cancer, vec lambda, vec kij, int n, int m, mat ci, mat Cistar, mat invCistar)
{
  vec Me(m);
  mat Va(m,m);
  vec menumerator(m);
  mat medenominator(m,m);
  vec t(m);
  mat ttinv(m,m); 
  int ii = 0;
  vec cij = ci.col(ii);
  mat tem = invCistar * cij;
  medenominator = tem * tem.t();
  menumerator = tem * (kij(ii)-lambda(ii));
  for (int ii = 1; ii < n; ii++)
  {
    cij = ci.col(ii);
    tem = invCistar * cij;
    menumerator = menumerator +  tem * (kij(ii)-lambda(ii));
    medenominator = tem * tem.t();
  }
  menumerator = menumerator/tausq;
  medenominator = medenominator/tausq + invCistar;
  Va = inv(medenominator);
  Me = Va * menumerator;
  mat wstar = rmvncpp(1, Me, Va);
  return(wstar);
  //return List::create(_["wstar"] = wstar,_["Me"] = Me,_["Va"] = Va);
  
}

// SVC algorithms:

// [[Rcpp::export]]
vec thetaupdatecpp(int n, int nx, mat w, vec B, List FF, uvec nnIndx, vec nnIndxLU, vec CIndx, mat iden)
{
  vec multivariate(n);
  // 1st pixel:
  rowvec mu0(nx);
  mu0.fill(0);
  mat ww(1,nx);
  ww.row(0) = w.row(0);
  vec tem1=logdmvnormcpp(ww,mu0,FF[0]);
  multivariate(0) = tem1(0);
  for (int i = 1; i < n; i++){
    uvec idx = nnIndx.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)) - 1;
    mat Btem = kron(B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t(),iden);
    mat wtem = vectorise(w.rows(idx),1);
    mat tempmean = Btem * wtem.t();
    rowvec meantem(nx);
    meantem.row(0) = tempmean.t();
    mat ww(1,nx);
    ww.row(0) = w.row(i);
    vec tem1=logdmvnormcpp(ww,meantem,FF[i]);
    multivariate(i) = tem1(0);
  }
  return(multivariate);
}

// [[Rcpp::export]]
List vupdatecpp(int n, int nx, mat w, vec B, vec FF, uvec nnIndx, vec nnIndxLU, vec CIndx, mat iden)
{
  // 1st pixel:
  List mats(n);
  mats[0] = w.row(0).t() * w.row(0) / FF(0);
  for (int i = 1; i < n; i++){
    uvec idx = nnIndx.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)) - 1;
    mat Btem = kron(B.rows(nnIndxLU(i),(nnIndxLU(i)+nnIndxLU(n+i)-1)).t(),iden);
    mat wtem=vectorise(w.rows(idx),1);
    mat tempmean = Btem * wtem.t();
    mats[i] = (w.row(i).t() - tempmean) * (w.row(i) - tempmean.t()) / FF(i);
  }
  return(mats);
}

// [[Rcpp::export]]
List svcmgupdatecpp(int n, int nx, mat V, vec FF, List FFmat, mat w, vec d, vec B,
                    vec nIndx, uvec nnIndx, vec nnIndxLU, vec uIndx, vec uIndxLU, vec uiIndx, mat iden)
{
  List m(n),M(n),G(n);
  vec Gvec(n);
  mat Vinv = solvecpp(V);
  // 1st element:
  int ii=0;
  mat zerovec(nx,1);
  zerovec.fill(0);
  m[ii]=zerovec;
  Gvec(ii)=1/FF(0);
  for (int jj = 0; jj < uIndxLU(n+ii); jj++)
  {
    int t = uIndx(uIndxLU(ii)+jj) - 1;
    int l = uiIndx(uIndxLU(ii)+jj) - 1;
    uvec indx = nnIndx.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)) - 1;
    //double a = w.row(t)-sum_2vec(B.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)), q.rows(indx));
    mat wtem = vectorise(w.rows(indx),1);
    mat a = w.row(t).t()-kron(B.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)).t(), iden) * wtem.t();
    mat Bscale(1,1); 
    Bscale(0,0)  = B(nnIndxLU(t)+l);
    mat BFinv(1,1);
    BFinv(0,0) = B(nnIndxLU(t)+l) / FF(t);
    mat BFinvkron = kron(BFinv, Vinv);
    a = a + kron(Bscale, iden) * w.row(indx(l)).t();
    mat mprevious = m[ii];
    m[ii] = mprevious + BFinvkron * a;
    Gvec(ii) = Gvec(ii) + pow(B(nnIndxLU(t)+l),2)/FF(t); // for G
  }
  mat Gscale(1,1);
  Gscale(0,0) = 1/Gvec(ii);
  G[ii] = kron(Gscale, V);
  mat gmattemp = G[ii];
  mat mtemp = m[ii];
  M[ii] = gmattemp * mtemp;
  
  // 2~nth element: Cholesky decomposition:
  for (int ii = 1; ii < n; ii++)
  {
    if (uIndxLU(n+ii) == 0)
    {
      uvec idx2 = nnIndx.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)) - 1;
      //double aa =  sum_2vec(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)), q.rows(idx2));
      mat wtem = vectorise(w.rows(idx2),1);
      mat aa = kron(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)).t(), iden) * wtem.t();
      mat Finv = solvecpp(FFmat[ii]);
      m[ii] =  Finv * aa;
      G[ii] = FFmat[ii];
      mat gmattemp = G[ii];
      mat mtemp = m[ii];
      M[ii] = gmattemp * mtemp;
    }
    if (uIndxLU(n+ii) > 0)
    {
      uvec idx2 = nnIndx.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)) - 1;
      //double aa =  sum_2vec(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)), q.rows(idx2));
      mat wtem = vectorise(w.rows(idx2),1);
      mat aa = kron(B.rows(nnIndxLU(ii),(nnIndxLU(ii)+nnIndxLU(n+ii)-1)).t(), iden) * wtem.t();
      mat Finv = solvecpp(FFmat[ii]);
      m[ii] =  Finv * aa;
      Gvec(ii) = 1/FF(ii);
      for (int jj = 0; jj < uIndxLU(n+ii); jj++)
      {
        int t = uIndx(uIndxLU(ii)+jj) - 1;
        int l = uiIndx(uIndxLU(ii)+jj) - 1;
        uvec indx = nnIndx.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)) - 1;
        mat wtem = vectorise(w.rows(indx),1);
        mat a = w.row(t).t()-kron(B.rows(nnIndxLU(t),(nnIndxLU(t)+nnIndxLU(n+t)-1)).t(), iden) * wtem.t();
        mat Bscale(1,1); 
        Bscale(0,0)  = B(nnIndxLU(t)+l);
        mat BFinv(1,1);
        BFinv(0,0) = B(nnIndxLU(t)+l) /FF(t);
        mat BFinvkron = kron(BFinv, Vinv);
        a = a + kron(Bscale, iden) * w.row(indx(l)).t();
        mat mprevious = m[ii];
        m[ii] = mprevious + BFinvkron * a;
        Gvec(ii) = Gvec(ii) + pow(B(nnIndxLU(t)+l),2)/FF(t); // for G
      }
      mat Gscale(1,1);
      Gscale(0,0) = 1/Gvec(ii);
      G[ii] = kron(Gscale, V);
      mat gmattemp = G[ii];
      mat mtemp = m[ii];
      M[ii] = gmattemp * mtemp;
    }
  }
  
  return List::create(_["M"] = M,_["G"] = G);
}

// [[Rcpp::export]]
List svcwsupdatecpp(int n, int nx, mat V, mat w, vec kij, List M, List G, List X, List loc_indx)
{
  List postmean(n);
  List postvar(n);
  for (int ii = 0; ii < n; ii++)
  {
    mat Xtemp = X[ii];
    mat Ginv = solvecpp(G[ii]);
    mat Covtemp = solvecpp( Xtemp.t() * Xtemp + Ginv);
    uvec locindex = loc_indx[ii];
    uvec loc_index = locindex -1;
    mat Mtem = M[ii];
    mat mtemp = Xtemp.t() * kij.rows(loc_index) + Ginv * Mtem;
    postvar[ii] = Covtemp;
    mat Mutemp = Covtemp * mtemp;
    postmean[ii] = Mutemp;
    vec wii = rmvncpp(1, Mutemp, Covtemp);
    w.row(ii) = wii.t();
  }
  return List::create(_["ws"] = w, _["postmean"] = postmean, _["postvar"] = postvar);
}

// [[Rcpp::export]]
List svckijupdatecpp(int nwhole, int n, mat ws, vec cancer, mat Xij, uvec indexij)
{
  vec kij(nwhole),meanij(nwhole),pstar(nwhole);
  for (int ii = 0; ii < nwhole; ii++)
  {
    mat Xtemp = Xij.row(ii);
    mat wij = ws.row(indexij(ii)-1);
    mat meanmat = Xtemp * wij.t();
    meanij(ii) = meanmat(0,0);
    //double meanij = sum_2vec(Xtemp, wij);
    //double meanij = sum_2vec(Xij.row(ii), ws.row(indexij(ii)-1));
    if (cancer(ii)==1)
    {
      kij(ii) = RcppTN::rtn1(meanij(ii), 1, 0, INFINITY);
    }
    if (cancer(ii)==0)
    {
      kij(ii) = RcppTN::rtn1(meanij(ii), 1, -INFINITY, 0);
    }
    pstar(ii) = 1-R::pnorm(0, meanij(ii), 1, 1, 0);
  }
  return List::create(_["kij"] = kij, _["meanij"] = meanij,_["pstar"] = pstar);
}


// [[Rcpp::export]]
vec svccancerupdatecpp(vec kij, vec pstar){
  vec a(kij.n_elem);
  for (int ii = 0; ii < kij.n_elem; ii++)
  {
    a(ii) = R::rbinom(1,pstar(ii));
  }
  return(a);
}
