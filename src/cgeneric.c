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
	// gcc -Wall -fpic -g -O -c -o cgeneric-demo.o cgeneric-demo.c
	// gcc -shared -o cgeneric-demo.so cgeneric-demo.o
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

		

        double *phi0 = Calloc(N, double);
        double *lambdas = Calloc(N, double);
        double *phi  = Calloc(N, double);
        double *kappas  = Calloc(N, double);
        
        double cc = 1/(N-1);
		//inla_cgeneric_vec_tp *psigma = data->doubles[0];
		//psigma->doubles[1])

		double time0; 
		time0 = (double)(time->doubles[0]);
		double timen; 
		timen = (double)(time->doubles[N-1]);
		assert(time0>0);
		assert(time0<2);
		assert(timen>0);
		assert(timen<2*N);
		double *time_norm  = Calloc(N, double);

		assert(time->doubles[0] > 0);
		assert(time->doubles[0] < 2*N);
		assert(time->doubles[0]>0);

		//cast time->doubles[0] to double
		assert( time->doubles[0] == time0 ); 
		
		//double tmp = time->doubles[0];
		double tmp;

		for(int i=0;i<N;i++){
			tmp = time->doubles[i];
			assert(tmp == time->doubles[i]);
			time_norm[i] = (tmp-time0)/(timen-time0);
			//assert(tmp-time_norm[i]);
			//assert(tmp-time0);
		}
		
		
        for(int i=0;i<N;i++){
            //time[i] = data->doubles[0]->doubles[i];
			phi0[i] = a+b*(time_norm[i]);
            lambdas[i] = -log(a+b * (time_norm[i]) );
            kappas[i] = kappa_eps*2*lambdas[i];

			/*
			double tmp = time_norm[i];
			
			assert(tmp<0);
			assert(a+b*tmp>0);
			assert(a+b*tmp<1);
			assert(lambdas[i]>0);

			assert(lambdas[i]>0);
			
			assert(-log(a+b*tmp)>0);
			*/
			
            
            if(i==0){
                phi0[i] = 0;
                phi[i] = 0;
            } else {
                phi[i] = exp(-lambdas[i] * ((time_norm[i])-(time_norm[i-1]))/cc);
            }
        }
		
		
		/*
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
		*/
		for (i = k = 0; i < N; i++) {
            if(i < N-1){
                ret[offset + k++] = kappas[i] + kappas[i+1]*phi[i]*phi[i];
                ret[offset + k++] = -kappas[i] * phi[i];
            }else{
                ret[offset + k++] = kappas[i];
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
        //ret = nullptr;
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



double *inla_cgeneric_timedep_forcing(inla_cgeneric_cmd_tp cmd, double *theta, inla_cgeneric_data_tp *data)
{
	// this reimplement `inla.rgeneric.ar1.model` using cgeneric

	//double *ret = NULL, prec, lprec, rho, rho_intern;
    double *ret = NULL, kappa_eps, lkappa_eps, b, b_intern, a, a_intern, kappa_f, lkappa_f, F0; //, kappa_f, lkappa_f, F0;

    assert(!strcasecmp(data->ints[0]->name, "n"));	       // this will always be the case
	int N = data->ints[0]->ints[0];			       // this will always be the case
	assert(N > 0);

    assert(!strcasecmp(data->doubles[0]->name, "time"));	       // time
	inla_cgeneric_vec_tp *time = data->doubles[0];
	assert(time->len == N);

	assert(!strcasecmp(data->doubles[1]->name, "forcing"));	       // time
	inla_cgeneric_vec_tp *forcing = data->doubles[1];
	assert(forcing->len == N);
	assert(forcing->doubles[N-1]);

	
	
    

	if (theta) {
		/*
		theta[0] = 1.3862944;
		theta[1] = 0.4054651;
		theta[2] = -0.5108256;
		theta[3] = 9.2103404;
		theta[4] = 0.0;

		*/
		
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
		lkappa_f = theta[3];
		kappa_f = exp(lkappa_f);
		F0 = theta[4];

		

        

	} else {
        kappa_eps = lkappa_eps = b = b_intern = a = a_intern = lkappa_f = kappa_f = F0 = NAN;
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

		

        double *phi0 = Calloc(N, double);
        double *lambdas = Calloc(N, double);
        double *phi  = Calloc(N, double);
        double *kappas  = Calloc(N, double);
        
        double cc = 1/(N-1);
		//inla_cgeneric_vec_tp *psigma = data->doubles[0];
		//psigma->doubles[1])

		double time0; 
		time0 = (double)(time->doubles[0]);
		double timen; 
		timen = (double)(time->doubles[N-1]);
		assert(time0>0);
		assert(time0<2);
		assert(timen>0);
		assert(timen<2*N);
		double *time_norm  = Calloc(N, double);

		assert(time->doubles[0] > 0);
		assert(time->doubles[0] < 2*N);
		assert(time->doubles[0]>0);

		//cast time->doubles[0] to double
		assert( time->doubles[0] == time0 ); 
		
		//double tmp = time->doubles[0];
		double tmp;

		for(int i=0;i<N;i++){
			tmp = time->doubles[i];
			assert(tmp == time->doubles[i]);
			time_norm[i] = (tmp-time0)/(timen-time0);
			assert(tmp-time_norm[i]);
			//assert(tmp-time0);
		}
		
		
        for(int i=0;i<N;i++){
            //time[i] = data->doubles[0]->doubles[i];
			phi0[i] = a+b*(time_norm[i]);
            lambdas[i] = -log(a+b * (time_norm[i]) );
            kappas[i] = kappa_eps*2*lambdas[i];

            if(i==0){
                phi0[i] = a;
                phi[i] = 1;
				phi[i] = exp(-lambdas[i]);
            } else {
                phi[i] = exp(-lambdas[i] * ((time_norm[i])-(time_norm[i-1]))/cc);
            }
        }
		
		
		for (i = k = 0; i < N; i++) {
            if(i < N-1){
				if(i == 0){
					//ret[offset + k++] = kappas[i]*(1-phi[i]*phi[i]) + kappas[i+1]*phi[i+1]*phi[i+1];
					ret[offset + k++] = kappas[i] + kappas[i+1]*phi[i+1]*phi[i+1];
                	ret[offset + k++] = -kappas[i+1] * phi[i+1];
				}else{
					ret[offset + k++] = kappas[i] + kappas[i+1]*phi[i+1]*phi[i+1];
                	ret[offset + k++] = -kappas[i+1] * phi[i+1];
				}
                
            }else{
                ret[offset + k++] = kappas[i];
            }
			
		}
		assert(ret[2+M-1]);
	}
		break;

	case INLA_CGENERIC_MU:
	{
		int offset = 1;
		ret = Calloc(N+offset, double);
		assert(ret);
		ret[0] = N;

		
		double *z  = Calloc(N, double);
		//double *phi0 = Calloc(N, double);
        double *lambdas = Calloc(N, double);
        //double *phi  = Calloc(N, double);
        //double *kappas  = Calloc(N, double);
		double *kappafs  = Calloc(N, double);
		double *time_norm  = Calloc(N, double);
		double tmp;

		double time0;
		time0 = (double)(time->doubles[0]);
		double timen; 
		timen = (double)(time->doubles[N-1]);
		
		/*
		a=0.3;
		b=0.2;
		F0=0.0;
		kappa_eps=4.0;
		kappa_f=10000.0;
		*/
		
		double tmpf;
		for(int i=0;i<N;i++){
			tmp = time->doubles[i];
			assert(tmp == time->doubles[i]);
			time_norm[i] = (tmp-time0)/(timen-time0);
			assert(time_norm[i]>=0);
			assert(time_norm[i]<=1);
			//phi0[i] = a+b*(time_norm[i]);
            lambdas[i] = -log(a+b * (time_norm[i]) );
            kappafs[i] = kappa_f*2*lambdas[i];
			
			tmpf = forcing->doubles[i];
			z[i] = (tmpf+F0)/sqrt(kappafs[i]);
			
		}
		/*
		assert(z[0] < 0.018);
		assert(z[0] > 0.017);
		assert(z[N-1] < -0.028);
		assert(z[N-1] > -0.029);

		*/
		
		/*
		assert(tmpf == forcing->doubles[N-1]);
		*/
	double tmpp;
		for(int i=0;i<N;i++){
			ret[i + offset] = 0;
			for(int s=0;s<=i;s++){
			//	tmp = z[s]*pow(M_E,-lambdas[i]*(time_norm[i]-time_norm[s]) );
				//tmp = z[s]*pow(M_E,-lambdas[i]*(time_norm[i]-time_norm[s]) );
				//assert((time_norm[i]-time_norm[s]) <=1);
				//assert((time_norm[i]-time_norm[s]) >=0);
				//assert(-lambdas[i] * (time_norm[i]-time_norm[s]) >=0);
				//assert(-lambdas[i] * (time_norm[i]-time_norm[s]) <= -log(a+b) );
				tmpp = (double)(z[s])* (double)(pow(M_E,- (double)(lambdas[i])*( (double)(time_norm[i])-(double)(time_norm[s]) ) )) ;
				//assert(tmpp<10);
				ret[i + offset] += tmpp;

				
			}
			//assert(tmp-time0);
			//assert(ret[i] < 10000);
			//assert(ret[i] < 1000);
		}
		/*
		assert(tmp == z[N-1]);
		assert(ret[0] == z[0]);
		assert(ret[0] < 0.0179);
		assert(ret[0] > 0.0178);
		assert(ret[N-1] < 2.29);
		assert(ret[N-1] > 2.28);
		*/
		
		
		/*
		
		//ret = Calloc(N, double);

		for(int i=0;i<N;i++){
			ret[i] = 0;
			for(int s=0;s<=i;s++){
				tmp = z[s]*pow(M_E,-lambdas[i]*(time_norm[i]-time_norm[s]) );
				ret[i] += tmp;
			}
			//assert(tmp-time0);
		}
		assert(tmp == z[N-1]);

		*/
		
		//assert(ret[0] == 0);
		// return (N, mu)
		// if N==0 then mu is not needed as its taken to be mu[]==0

		
		//assert(ret[N-1]);
	}
		break;

	case INLA_CGENERIC_INITIAL:
	{
		// return c(M, initials)
		// where M is the number of hyperparameters

		ret = Calloc(6, double);
		assert(ret);
		ret[0] = 5;
		ret[1] = 0.1;
		ret[2] = 0.1;
        ret[3] = 0.1;
		ret[4] = 0.1;
		ret[5] = 0.1;
		assert(ret[5]);
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
        //ret = nullptr;
	}
		break;

	case INLA_CGENERIC_LOG_PRIOR:
	{
		// return c(LOG_PRIOR)

		ret = Calloc(1, double);
		assert(ret);
		//ret[0] = -prec + lprec - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(rho_intern);

		double k1 = 1/SQR(0.5);
     	double k2 = 1/SQR(0.5);
     	double k3 = 1/SQR(0.5);
     	double k4 = 1/SQR(0.5);
     	double k5 = 1/SQR(0.5);
		double lprior;

		lprior = -2.5*log(2.0*M_PI) + 0.5*log(k1*k2*k3*k4*k5) -0.5*(k1 * SQR(lkappa_eps) +k2*SQR(b_intern) +k3 * SQR(a_intern) +k4 * SQR(lkappa_f) +k5 * SQR(F0));
		 // -b_intern +log(1-pow(M_E,-b_intern)) -2*log(1-pow(M_E,-a_intern));
		ret[0] = lprior;
		/*
		lprior = -2.5*log(2.0*M_PI) - log(s1*s2*s3*s4*s5) - 0.5 * (  )
	  0.5*(s1*theta[1]^2 + s2*theta[2]^2 + s3*theta[3]^2 +
			 s4*theta[4]^2 + s5*theta[5]^2);

    lprior = -2.5*log(2.0*pi) - log(s1*s2*s3*s4*s5) -
      0.5*(s1*theta[1]^2 + s2*theta[2]^2 + s3*theta[3]^2 +
             s4*theta[4]^2 + s5*theta[5]^2)


        ret[0] = -kappa_eps + lkappa_eps - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(a_intern) - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(b_intern) -kappa_f + lkappa_f - 0.5 * log(2.0 * M_PI) - 0.5 * SQR(F0);
	
		*/
	}
		
		break;

	case INLA_CGENERIC_QUIT:
	default:
		break;
	}

	return (ret);
}