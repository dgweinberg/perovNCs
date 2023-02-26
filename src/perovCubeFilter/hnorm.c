#include "fd.h"

/*****************************************************************************/

void hnorm(zomplex *psim1,zomplex *psin,double *potl,double *ksqr,par_st par,nlc_st *nlc,long *nl,double zm1,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long i;
  hamiltonian(psin,psim1,potl,ksqr,ist,par,nlc,nl,planfw,planbw,fftwpsi);

  for (i = 0; i < ist.nspinngrid; i++){
    /*** par.dE_1 = 4.0 / par.dE and therefore I don't multiply by 4 ***/
    psin[i].re = par.dE_1 * psin[i].re - (2.0 + zm1 + par.Vmin * par.dE_1) * psim1[i].re;
    psin[i].im = par.dE_1 * psin[i].im - (2.0 + zm1 + par.Vmin * par.dE_1) * psim1[i].im;
  }
  
  return;
}

/*****************************************************************************/


/*****************************************************************************/
