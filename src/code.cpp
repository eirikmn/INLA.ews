#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
#include <R.h>

using namespace Rcpp;

//' Compute mu in C for AR(1) model
//' 
//' Void function that gives the forced response for the AR(1) model with 
//' time-dependent memory in C.
//' 
//' @param means Vector that will be overwritten with the forced response.
//' @param z Vector with the shifted and scaled known forcing.
//' @param n Integer denoting the dimension of \code{mu}.
//' @param lambda Vector describing the lambdas at each sampled time step.
//' @param time Vector describing the time steps.
//' 
//' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
//' @seealso \code{\link{mucomputer}}
//' @keywords C ar1 forcing
//' @export
// [[Rcpp::export]]
void c_mu_ar1(NumericVector means, NumericVector z,int n, NumericVector lambda, NumericVector time)
{
    for(int i=0;i<n;i=i+1){
        for(int j=0;j<=i;j=j+1){
			    //means[i] +=    z[j]*structure[i-j];
			    means[i] += z[j]*pow(M_E,-lambda[i]*(time[i]-time[j]) );
		  }
    }
    
	
}
