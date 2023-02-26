#include "fd.h"

/*****************************************************************************/

long portho(MKL_Complex16 *psi,double dv,long_st ist)
{
  long long lwork;  long long info, one=1, i, cutoff;
  long long ngrid = (long long)(ist.nspinngrid), mstot = (long long)(2*ist.mstot);
  double *S, *rwork; MKL_Complex16 *work;
  
  lwork = 5*(long long)(mstot*mstot+ngrid);
  S = (double*) malloc(mstot * sizeof(double));
  work = (MKL_Complex16*) malloc(lwork * sizeof(MKL_Complex16));
  rwork = (double*) malloc(5*mstot * sizeof(double));


  //Call lapack function
  zgesvd_("O","N",&(ngrid),&(mstot),&(psi[0]),&(ngrid),&(S[0]),
	  NULL,&one,NULL,&one,&(work[0]),&(lwork),&(rwork[0]),&info);
  if (info != 0) {printf("error in zgesvd(1) %lld, exiting",info); exit(0);}

  for (cutoff = mstot, i=0; i<mstot; i++) {
    if ((S[i] / S[0]) < 1.0e-10) {
      cutoff = i;
      break;
    }
  }
  printf("cutoff is %lld\n",cutoff);
  
  free(work); 
  free(S);
  free(rwork);

  return (cutoff);
}

/*****************************************************************************/

