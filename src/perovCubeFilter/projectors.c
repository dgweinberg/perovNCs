/*****************************************************************************/
//This file contains the function for generating the projectors used for the
//non-local part of the hamiltonian

#include "fd.h"

/*****************************************************************************/
#define calcBessel(x,x1)   ((x) < EPS ? 0 : (sin((x)) * ((x1)*(x1)) - cos((x)) * (x1)))

#define EPS 1e-10
/*****************************************************************************/
void gen_SO_projectors(double dx, double rcut, long nproj, double*  projectors, double* vr){
	long long N  = 1024;
	double kmax = 1.5 * PIE/dx;
	double dk = kmax/((double) N);
	double k, ki, kj;
	double dr = rcut / ((double) N);
	double r;
	double sum, preFactor;
	int i,j,rpoint,kpoint,projector;

	// Width of the gaussian fucntion used for the SO potential. 
	double width = 0.7;



	double *A = calloc(N*N, sizeof(double));	

	for ( rpoint = 0; rpoint < N; rpoint++){
		vr[rpoint] = (double) rpoint * dr;
	}
	
	//calculate the matrix of integrals for each ki and kj (symmetry assisted for i<-->j)
	#pragma omp parallel for private(i,j,rpoint,ki,kj,sum,r)
	for ( i = 0; i < N; i++) {
		ki = (double) i * dk;
		for ( j = 0; j <= i; j++) {
		  	kj = (double) j * dk;
		  	
		  	// perform numerical integration 
		  	sum = 0.0;
		  	for ( rpoint = 1; rpoint < N; rpoint++) {
		    	r = vr[rpoint];
		    	sum += (sqr(r) * calcBessel(ki*r, (1.00/(ki*r + EPS)) )
		      		* exp(-1.0 * (sqr(r/width))) * calcBessel(kj*r, (1.00/(kj*r + EPS))) * dr);                       
	  		}
		  A[i*N+j] = sum;
		  A[j*N+i] = sum;
		}
	}

	//setup call to lapack
	char JOBZ = 'V';
	char UPLO = 'U';
	long long LDA = N;
	double* W = calloc( N, sizeof(double));
	double* WORK = calloc( 3 * N, sizeof(double));
	long long LWORK = 3 * N;
	long long INFO = 0;
	//diagonalize the matrix
	dsyev_(&JOBZ, &UPLO, &N, &A[0], &LDA, &W[0], &WORK[0], &LWORK, &INFO);


	for ( projector = 0; projector <nproj; projector++){
		preFactor =  (2.00 / PIE) * sqrt(W[N-projector-1]);
		for ( rpoint = 0; rpoint < N; rpoint++){
			r = vr[rpoint];
			
			//perform numerical integration
			sum = 0;
			#pragma omp parallel for private(k,kpoint) reduction(+:sum)
			for ( kpoint = 0; kpoint < N; kpoint++){
				k = (double) kpoint * dk;
				sum += k * k * calcBessel(k*r, 1/(k*r + EPS)) * A[N*(N-projector-1)+kpoint] * dk;
			}
			projectors[N * projector + rpoint] = sum * preFactor;
			
		}
	}

	free(W); free(WORK); free(A);
	return;
}
/*****************************************************************************/
void gen_nlc_projectors(double dx, double rcut, long nproj, double*  projectors,int* sgnProj, double* vr, atm_st *atm,long jatom,  long_st ist){
	long long N  = 1024;
	double kmax = 1.5 * PIE/dx;
	double dk = kmax/((double) N);
	double k, ki, kj;
	double dr = rcut / ((double) N);
	double r;
	double sum, preFactor;
	int i,j,rpoint,kpoint,projector;
	double lam1, lam2;

	// Width of the gaussian fucntion used for the SO potential. 
	double width = 1.0;
	double shift = 1.5;
	double *A = calloc(N*N, sizeof(double));
	double* W = calloc( N, sizeof(double));
	double* WORK = calloc( 3 * N, sizeof(double));
	
	lam1= atm[jatom].nlcPar[0];
	lam2= atm[jatom].nlcPar[1];
			
	printf("For atom %ld with Zval %d we have l1=%g l2=%f\n", jatom, atm[jatom].Zval, lam1, lam2);
	if (lam1==0 && lam2==0){
		for ( projector = 0; projector <nproj; projector++){
			for(rpoint=0;rpoint<N;rpoint++){
				projectors[ N * projector + rpoint] = 0.0;
				sgnProj[projector] = 0.0;
			}

		}
		return;
	}

	for ( rpoint = 0; rpoint < N; rpoint++){
		vr[rpoint] = (double) rpoint * dr;
	}

	//calculate the matrix of integrals for each ki and kj (symmetry assisted for i<-->j)
	#pragma omp parallel for private(i,j,rpoint,ki,kj,sum,r)
	for ( i = 0; i < N; i++) {
		ki = (double) i * dk;
		for ( j = 0; j <= i; j++) {
		  	kj = (double) j * dk;
		  	
		  	// perform numerical integration 
		  	sum = 0.0;
		  	for ( rpoint = 1; rpoint < N; rpoint++) {
		    	r = vr[rpoint];
		    	sum += (sqr(r) * calcBessel(ki*r, (1.00/(ki*r + EPS)) )
		      		* (lam1*exp(-1.0 * (sqr(r/width))) +lam2*exp(-(sqr((r-shift)/width))))* calcBessel(kj*r, (1.00/(kj*r + EPS))) * dr);                       
	  		}
		  A[i*N+j] = sum;
		  A[j*N+i] = sum;
		}
	}

	//setup call to lapack
	char JOBZ = 'V';
	char UPLO = 'U';
	long long LDA = N;
	long long LWORK = 3 * N;
	long long INFO = 0;
	//diagonalize the matrix
	dsyev_(&JOBZ, &UPLO, &N, &A[0], &LDA, &W[0], &WORK[0], &LWORK, &INFO);
	


	double eigs[nproj];
	int i_eigs[nproj];
	for (i=0;i<nproj;i++){
		i_eigs[i]=0;
		eigs[i]=0.0;
	}


	for (i=0;i<N;i++){
		if (fabs(W[i])>eigs[nproj-1]){
			for(j=nproj-1;j>=0;j--){
				if(j==0 || eigs[j-1]>fabs(W[i])){
					eigs[j] = fabs(W[i]);
					i_eigs[j]= i;
					break;
				}
				else{
					eigs[j]=eigs[j-1];
					i_eigs[j]=i_eigs[j-1];
				}
			} 
		}
	}



	for ( projector = 0; projector <nproj; projector++){
		if(W[i_eigs[projector]]<0.0){
			W[i_eigs[projector]]=fabs(W[i_eigs[projector]]);
			sgnProj[projector]=-1;
		}
		else{
			sgnProj[projector]=1;
		}

		preFactor =  (2.00 / PIE) * sqrt(W[i_eigs[projector]]);
		for ( rpoint = 0; rpoint < N; rpoint++){
			r = vr[rpoint];
			
			//perform numerical integration
			sum = 0;
			#pragma omp parallel for private(k,kpoint) reduction(+:sum)
			for ( kpoint = 0; kpoint < N; kpoint++){
				k = (double) kpoint * dk;
				sum += k * k * calcBessel(k*r, 1/(k*r + EPS)) * A[N*(i_eigs[projector])+kpoint] * dk;
			}
			projectors[N * projector + rpoint] = sum * preFactor;
		}
	}
	

	free(W); free(WORK); free(A);
	return;
}
/*****************************************************************************/
