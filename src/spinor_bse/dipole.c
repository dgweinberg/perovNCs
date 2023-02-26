/****************************************************************************/

#include "fd.h"

/****************************************************************************/

void dipole(double *vx,double *vy,double *vz,zomplex *psi,zomplex *mux,zomplex *muy,zomplex *muz,double *eval,long_st ist,par_st par)
{
  FILE *pf1, *pf2, *pf3, *pf; 
  long a, i, jx, jy, jz, jgridup,jgriddn, jyz; 
  zomplex sumX, sumY, sumZ, tmp;
  double dz, dx, dy, ev, os;
  
  // Output will be written to these files
  pf = fopen("OS0.dat" , "w"); 
  pf1 = fopen("ux.dat", "w"); pf2 = fopen("uy.dat", "w"); pf3 = fopen("uz.dat", "w");

  // Make sure mux, muy and muz are zero to begin with
  for (i = 0; i < ist.totallumo*ist.totalhomo; i++) mux[i].re = muy[i].re = muz[i].re = mux[i].im = muy[i].im = muz[i].im = 0.0;

  // Main computional work of function performed here - must loop over all electron-hole (i-a) pairs
  // compute $$\vec{\mu} = <i|\vec{r}|a> = \int dr \sum_\simga \psi_{i}^{*}(r,\sigma) \vec{r} \psi_a*(r,\sigma)$$
  for (i = 0; i < ist.totalhomo; i++, fprintf(pf1,"\n"),fprintf(pf2,"\n"),fprintf(pf3,"\n")){
    for (a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
      sumX.re = sumY.re = sumZ.re = sumX.im = sumY.im = sumZ.im = 0.0;
      for (jz = 0; jz < ist.nz; jz++) {
        dz = vz[jz];
        for (jy = 0; jy < ist.ny; jy++) {
          dy = vy[jy];
          jyz = ist.nx * (ist.ny * jz + jy);
          for (jx = 0; jx < ist.nx; jx++) {
            dx = vx[jx];
            jgridup = jyz + jx;
            jgriddn = jgridup+ist.ngrid;
      	    tmp.re =  par.dv*(psi[i*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgridup].re + psi[i*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgridup].im
                             +psi[i*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgriddn].re + psi[i*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgriddn].im );
            tmp.im = par.dv * (-psi[i*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgridup].re + psi[i*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgridup].im
                              -psi[i*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgriddn].re + psi[i*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgriddn].im) ;
      	    

            sumX.re +=  dx*tmp.re; sumX.im +=  dx*tmp.im;
      	    sumY.re +=  dy*tmp.re; sumY.im +=  dy*tmp.im;
      	    sumZ.re +=  dz*tmp.re; sumZ.im +=  dz*tmp.im;
          }
        }
      }
      mux[i*ist.totallumo+(a-ist.nlumo)].re = sumX.re;  mux[i*ist.totallumo+(a-ist.nlumo)].im = sumX.im; 
      muy[i*ist.totallumo+(a-ist.nlumo)].re = sumY.re;  muy[i*ist.totallumo+(a-ist.nlumo)].im = sumY.im;
      muz[i*ist.totallumo+(a-ist.nlumo)].re = sumZ.re;  muz[i*ist.totallumo+(a-ist.nlumo)].im = sumZ.im;
      fprintf (pf1,"%ld %ld %g %g\n",i,a,mux[i*ist.totallumo+(a-ist.nlumo)].re, mux[i*ist.totallumo+(a-ist.nlumo)].im);
      fprintf (pf2,"%ld %ld %g %g\n",i,a,muy[i*ist.totallumo+(a-ist.nlumo)].re, muy[i*ist.totallumo+(a-ist.nlumo)].im);
      fprintf (pf3,"%ld %ld %g %g\n",i,a,muz[i*ist.totallumo+(a-ist.nlumo)].re, muz[i*ist.totallumo+(a-ist.nlumo)].im);

      os=((sqr(sumX.re)+sqr(sumX.im)) + (sqr(sumY.re)+sqr(sumY.im)) + (sqr(sumZ.re)+sqr(sumZ.im)));
      ev = eval[a] - eval[i];
      fprintf(pf,"%ld %ld %.8f %.12f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n", i, a, sqrt(os), ev, (2.0/3.0)*ev*os,
	       mux[i*ist.totallumo+(a-ist.nlumo)].re, mux[i*ist.totallumo+(a-ist.nlumo)].im,
	       muy[i*ist.totallumo+(a-ist.nlumo)].re, muy[i*ist.totallumo+(a-ist.nlumo)].im,
	       muz[i*ist.totallumo+(a-ist.nlumo)].re, muz[i*ist.totallumo+(a-ist.nlumo)].im);
    }
  }
  fclose(pf); fclose(pf1); fclose(pf2); fclose(pf3);

  return;
}


/****************************************************************************/
