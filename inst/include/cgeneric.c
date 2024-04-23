#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>

#include "cgeneric.h"

#define Calloc(n_, type_)  (type_ *)calloc((n_), sizeof(type_))
#define SQR(x) ((x)*(x))


double *inla_cgeneric_ar1_model(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
	// this reimplement `inla.rgeneric.ar1.model` using cgeneric

	double *ret = NULL, prec, lprec, rho, rho_intern;

	if (theta) {
		lprec = theta[0];
		prec = exp(lprec);
		rho_intern = theta[1];
		rho = 2.0 * exp(rho_intern) / (1.0 + exp(rho_intern)) - 1.0;
	} else {
		prec = lprec = rho = rho_intern = NAN;
	}

	assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	int N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
	}
		break;

	case INLA_CGENERIC_GRAPH:
	{
		// return a vector of indices with format
		// c(N, M, ii, jj)
		// where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
		// so like the loop
		// for i=0, ...
		// for j=i, ...
		// G_ij = 
		// and M is the length of ii

		int M = N + N - 1, offset, i, k;
		ret = Calloc(2 + 2 * M, double);
		assert(ret);
		offset = 2;
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
		for (k = i = 0; i < N; i++) {
			ret[offset + k] = i;		       /* i */
			ret[offset + M + k++] = i;	       /* j */
			if (i < N - 1) {
				ret[offset + k] = i;	       /* i */
				ret[offset + M + k++] = i + 1; /* j */
			}
		}
	}
		break;

	case INLA_CGENERIC_Q:
	{
		// optimized format
		// return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
		// where M is the length of Qij

		double param = prec / (1.0 - SQR(rho));
		int M = N + N - 1;
		int offset, i, k;
		ret = Calloc(2 + M, double);
		assert(ret);
		offset = 2;
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;
		for (i = k = 0; i < N; i++) {
			ret[offset + k++] = param * (i == 0 || i == N - 1 ? 1.0 : (1.0 + SQR(rho)));
			if (i < N - 1) {
				ret[offset + k++] = -param * rho;
			}
		}
	}
		break;

	case INLA_CGENERIC_MU:
	{
		// return (N, mu)
		// if N==0 then mu is not needed as its taken to be mu[]==0

		ret = Calloc(1, double);
		assert(ret);
		ret[0] = 0;
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters

		ret = Calloc(3, double);
		assert(ret);
		ret[0] = 2;
		ret[1] = 1.0;
		ret[2] = 1.0;
	}
		break;

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		// return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself

		double prec_innovation = prec / (1.0 - SQR(rho));
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = N * (-0.5 * log(2.0 * M_PI) + 0.5 * log(prec_innovation)) + 0.5 * log(1.0 - SQR(rho));
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR)

		ret = Calloc(1, double);
		assert(ret);
		ret[0] = -prec + lprec - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(rho_intern);
	}
		break;

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}




double *inla_cgeneric_timedep(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
	// this reimplement `inla.rgeneric.ar1.model` using cgeneric

	//double *ret = NULL, prec, lprec, rho, rho_intern;
    double *ret = NULL, kappa_eps, lkappa_eps, b, b_intern, a, a_intern; //, kappa_f, lkappa_f, F0;

    assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	int N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

    assert(!strcasecmp(data->doubles[0]->name, "time"));	       // time
	inla_cgeneric_vec_tp *time = data->doubles[0];
	assert(time->len == N);

    

	if (theta) {
		lkappa_eps = theta[0];
		kappa_eps = exp(lkappa_eps);
		b_intern = theta[1];
        b = -1.0 + 2.0/(1.0+exp(-b_intern));
		//rho = 2.0 * exp(rho_intern) / (1.0 + exp(rho_intern)) - 1.0;
        a_intern = theta[2];
        if(b<0){
            a = -b + (1+b)/(1+exp(-a_intern));
        }else{
            a = 0 + (1-b)/(1+exp(-a_intern));
        }

        double *time = Calloc(N, double);
        double *phi0 = Calloc(N, double);
        double *lambdas = Calloc(N, double);
        double *phi  = Calloc(N, double);
        double *kappas  = Calloc(N, double);
        
        double cc;
        cc = 1/(N-1);
        for(int i=0;i<N;i++){
            time[i] = data->doubles[0]->doubles[i];
            lambdas[i] = -log(a+b*time[i]);
            kappas[i] = kappa_eps*2*lambdas[i];
            phi0[i] = a+b*time[i];
            if(i==0){
                phi0[i] = 0;
                phi[i] = 0;
            } else {
                phi[i] = exp(-lambda[i] * (time[i]-time[i-1])/cc);
            }
        }

	} else {
        kappa_eps = lkappa_eps = b = b_intern = a = a_intern = NAN;
		//prec = lprec = rho = rho_intern = NAN;
	}

	

	switch (cmd) {
	case INLA_CGENERIC_VOID:
	{
		assert(!(cmd == INLA_CGENERIC_VOID));
	}
		break;

	case INLA_CGENERIC_GRAPH:
	{
		// return a vector of indices with format
		// c(N, M, ii, jj)
		// where ii<=jj, ii is non-decreasing and jj is non-decreasing for the same ii
		// so like the loop
		// for i=0, ...
		// for j=i, ...
		// G_ij = 
		// and M is the length of ii

		int M = N + N - 1, offset, i, k;
		ret = Calloc(2 + 2 * M, double);
		assert(ret);
		offset = 2;
		ret[0] = N;				       /* dimension */
		ret[1] = M;				       /* number of (i <= j) */
		for (k = i = 0; i < N; i++) {
			ret[offset + k] = i;		       /* i */
			ret[offset + M + k++] = i;	       /* j */
			if (i < N - 1) {
				ret[offset + k] = i;	       /* i */
				ret[offset + M + k++] = i + 1; /* j */
			}
		}
	}
		break;

	case INLA_CGENERIC_Q:
	{
		// optimized format
		// return c(-1, M, Qij) in the same order as defined in INLA_CGENERIC_GRAPH
		// where M is the length of Qij

		//double param = prec / (1.0 - SQR(rho));
		int M = N + N - 1;
		int offset, i, k;
		ret = Calloc(2 + M, double);
		assert(ret);
		offset = 2;
		ret[0] = -1;				       /* REQUIRED */
		ret[1] = M;
		for (i = k = 0; i < N; i++) {
            if(i < N-1){
                ret[offset + k++] = kappas[i] + kappas[i+1]*phi[i]*phi[i];
                ret[offset + k++] = -kappas[i] * phi[i];
            }else{
                ret[offset + k++] = kappas[i];
            }
            
			
            //ret[offset + k++] = param * (i == 0 || i == N - 1 ? 1.0 : (1.0 + SQR(rho)));
			//if (i < N - 1) {
		//		ret[offset + k++] = -param * rho;
		//	}
		}
	}
		break;

	case INLA_CGENERIC_MU:
	{
		// return (N, mu)
		// if N==0 then mu is not needed as its taken to be mu[]==0

		ret = Calloc(1, double);
		assert(ret);
		ret[0] = 0;
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters

		ret = Calloc(3, double);
		assert(ret);
		ret[0] = 3;
		ret[1] = 1.0;
		ret[2] = 1.0;
        ret[3] = 1.0;
	}
		break;

	case INLA_CGENERIC_LOG_NORM_CONST:
	{
		// return c(NORM_CONST) or a NULL-pointer if INLA should compute it by itself
        /*
        double prec_innovation = prec / (1.0 - SQR(rho));
		ret = Calloc(1, double);
		assert(ret);
		ret[0] = N * (-0.5 * log(2.0 * M_PI) + 0.5 * log(prec_innovation)) + 0.5 * log(1.0 - SQR(rho));
        */
		//ret = Calloc(1, double);
		//assert(ret);
		//ret[0] = 0;
        ret = nullptr;
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR)

		ret = Calloc(1, double);
		assert(ret);
		//ret[0] = -prec + lprec - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(rho_intern);
        ret[0] = -kappa_eps + lkappa_eps - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(a_intern) - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(b_intern);
	}
		break;

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}