#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
using namespace Rcpp;

//' Compute mu in C for AR(1) model
//' 
//' Void function that gives the forced response for the AR(1) model with 
//' time-dependent memory in C.
//' 
//' @param mu Vector that will be overwritten with the forced response.
//' @param n Integer denoting the dimension of \code{mu}.
//' @param tforcing Vector with the shifted and scaled known forcing.
//' @param phis Vector describing the time-dependent lag-one correlation.
//' 
//' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
//' @seealso \code{\link{mu.computer}}
//' @keywords C ar1 forcing
//' @export
// [[Rcpp::export]]
void compute_mu_ar1(NumericVector mu, int n, NumericVector tforcing, 
                    NumericVector phis) {
  
  for(int i=0;i<n;i=i+1){
    mu[i] = 0;
    double lambda = phis[i]-1;
    for(int j=0; j<=i; j=j+1){
      mu[i] += pow(M_E,lambda*((double)i-(double)j +0.5 ))*tforcing[j];
    }
  }
  
}

//' Compute mu in C for fGn model
//' 
//' Void function that gives the forced response for the fractional Gaussian noise 
//' model with time-dependent memory in C.
//' 
//' @param mu Vector that will be overwritten with the forced response.
//' @param n Integer denoting the dimension of \code{mu}.
//' @param tforcing Vector with the shifted and scaled known forcing.
//' @param Hs Vector describing the time-dependent Hurst exponent.
//' 
//' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
//' @seealso \code{\link{mu.computer}}
//' @keywords C fgn forcing
//' @export
// [[Rcpp::export]]
void compute_mu_fgn(NumericVector mu, int n, NumericVector tforcing, 
                    NumericVector Hs) {
  
  for(int i=0;i<n;i=i+1){
    mu[i] = 0;
    for(int j=0; j<=i;j=j+1){
      mu[i] += pow((double)i-(double)j+0.5,Hs[i]-1.5)*tforcing[j];
    }
  }
  
}

//' Multiply a number by two
//' 
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}