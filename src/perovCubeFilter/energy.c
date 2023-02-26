#include "fd.h"

/*****************************************************************************/

double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,nlc_st *nlc,long *nl,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long i; 
  double ene = 0.0;

  memcpy(&phi[0], &psi[0], ist.nspinngrid*sizeof(phi[0]));
  hamiltonian(phi, psi, potl, ksqr, ist, par, nlc, nl, planfw, planbw, fftwpsi);

  for (i = 0; i < ist.nspinngrid; i++) {
    ene += (psi[i].re * phi[i].re + psi[i].im * phi[i].im);
  }
  ene *= par.dv;

  return (ene);
}

/***************************************************************************/

void energy_all(zomplex *psi,zomplex *phi,zomplex *psims,double *potl,double *ksqr,double *ene,long_st ist,par_st par,nlc_st *nlc,long *nl,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long jgrid, jms, jmsg;

  for (jms = 0; jms < ist.ms; jms++) {
    jmsg = jms * ist.nspinngrid;

    for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
      psi[jgrid].re = psims[jmsg+jgrid].re;
      psi[jgrid].im = psims[jmsg+jgrid].im;
    }
    
    memcpy(&phi[0], &psi[0], ist.nspinngrid*sizeof(phi[0]));
    hamiltonian(phi, psi, potl, ksqr, ist, par, nlc, nl, planfw, planbw, fftwpsi);
    
    for (ene[jms] = 0.0, jgrid = 0; jgrid < ist.nspinngrid; jgrid++)
      ene[jms] += (psi[jgrid].re * phi[jgrid].re + psi[jgrid].im * phi[jgrid].im);
    
    ene[jms] *= par.dv;
  }

  return;
}

/***************************************************************************/

void get_energy_range(zomplex *psi,zomplex *phi,double *potl,double *vx,double *vy,double *vz,double *ksqr,long_st ist,par_st *par,nlc_st *nlc,long *nl,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  FILE *pf;
  long i, ispn, jgrid; 
  long idum = -874917403;
  double ene_old; 
  double norma, Emin, Emax, tau = 0.025;

  /*** calculate the energy range ***/
  /*** Emin ***/
  pf = fopen("Emin-init.dat" , "w");
  for (ispn = 0; ispn < ist.nspin; ispn++)  {
    init_psi(&phi[ispn*ist.ngrid], ist, *par, &idum);
  }
  clock_t timers[4];

  Emin = (ene_old = 0.0) + 0.1;
  for (i = 0; (fabs((Emin-ene_old)/Emin)>1.0e-6) && (i < 500); i++){
    memcpy(&psi[0], &phi[0], ist.nspinngrid*sizeof(phi[0]));
    hamiltonian_t(phi, psi, potl, ksqr, ist, *par, nlc, nl, planfw, planbw, fftwpsi, &timers[0]);
    for (ispn = 0; ispn < ist.nspin; ispn++) {
      for (jgrid = 0; jgrid < ist.ngrid; jgrid++) {
        phi[ispn*ist.ngrid+jgrid].re = psi[ispn*ist.ngrid+jgrid].re - tau*phi[ispn*ist.ngrid+jgrid].re;
        phi[ispn*ist.ngrid+jgrid].im = psi[ispn*ist.ngrid+jgrid].im - tau*phi[ispn*ist.ngrid+jgrid].im;
      }
    }
    norma = normalize(&phi[0], par->dv, ist.nspinngrid, ist.nthreads);

    ene_old = Emin;
    Emin = energy(phi, psi, potl, ksqr, ist, *par, nlc, nl, planfw, planbw, fftwpsi);
    fprintf(pf, "%ld %.16g %.16g %.16g\n", i, ene_old, Emin, norma);
    fflush(pf);
  }
  fclose(pf);

  /*** Emax ***/
  for (ispn = 0; ispn < ist.nspin; ispn++) init_psi(&psi[ispn*ist.ngrid],ist,*par,&idum);
  Emax = (ene_old = 0.0) + 0.1;
  pf = fopen("Emax-init.dat" , "w");
  for (i = 0; (fabs((Emax-ene_old)/Emax)>1.0e-6) & (i < 500); i++){
    memcpy(&phi[0], &psi[0], ist.nspinngrid*sizeof(phi[0]));
    hamiltonian_t(psi, phi, potl, ksqr, ist, *par, nlc, nl, planfw, planbw, fftwpsi, &timers[0]);
    norma = normalize(&psi[0],par->dv,ist.nspinngrid,ist.nthreads);

    ene_old = Emax;
    Emax = energy(psi,phi,potl,ksqr,ist,*par,nlc,nl,planfw,planbw,fftwpsi);
    fprintf (pf,"%ld %.16g %.16g %.16g\n",i,ene_old,Emax,norma);
    fflush(pf);
  }
  fclose(pf);
  printf("Timing hamiltonian:\n");
  int KEtime = 1000*timers[0]/CLOCKS_PER_SEC;
  int SOPEtime = 1000*timers[1]/CLOCKS_PER_SEC;
  int NLPEtime = 1000*timers[2]/CLOCKS_PER_SEC;
  int LOCPEtime = 1000*timers[3]/CLOCKS_PER_SEC;
  int TOTtime = KEtime+SOPEtime+NLPEtime+LOCPEtime;
  printf("Kinetic time: %d msec / %d msec\n",KEtime,TOTtime);
  printf("SO Pot  time: %d msec / %d msec\n",SOPEtime,TOTtime);
  printf("NonLoc  time: %d msec / %d msec\n",NLPEtime,TOTtime);
  printf("Loc pot time: %d msec / %d msec\n",LOCPEtime,TOTtime);
  Emax *= 1.2;
  Emin -= 0.2 *fabs(Emin);
  
  par->Vmin = Emin;
  par->dE = (Emax - Emin);
  par->dE_1 = 4.0 / par->dE;

  printf("Emax, Emin, dE = %g %g %g\n", Emax, Emin, par->dE);
  fflush(stdout);

  return;
}


/****************************************************************************************/

void calc_sigma_E(zomplex *psi,zomplex *phi,zomplex *psitot,double *potl,double *ksqr,double *eval2,long_st ist,par_st par,nlc_st *nlc,long *nl ,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long jgrid, ims;
  double eval;
  
  for (ims = 0; ims < ist.mstot; ims++) {
    
    for (jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
      psi[jgrid].re = psitot[ims*ist.nspinngrid+jgrid].re;
      psi[jgrid].im = psitot[ims*ist.nspinngrid+jgrid].im;
    }
    memcpy(&phi[0],&psi[0],ist.nspinngrid*sizeof(phi[0]));
    hamiltonian(phi,psi,potl,ksqr,ist,par,nlc,nl,planfw,planbw,fftwpsi);

    for (eval = 0.0, jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
      eval += (psitot[ims*ist.nspinngrid+jgrid].re * phi[jgrid].re + psitot[ims*ist.nspinngrid+jgrid].im * phi[jgrid].im);
    }
    eval *= par.dv;
    
    memcpy(&psi[0], &phi[0], ist.nspinngrid*sizeof(psi[0]));
    hamiltonian(phi, psi, potl, ksqr, ist, par, nlc, nl, planfw, planbw, fftwpsi);

    for (eval2[ims] = 0.0, jgrid = 0; jgrid < ist.nspinngrid; jgrid++) {
      eval2[ims] += (psitot[ims*ist.nspinngrid+jgrid].re * phi[jgrid].re + psitot[ims*ist.nspinngrid+jgrid].im * phi[jgrid].im);
    }
    eval2[ims] *= par.dv;
    
    eval2[ims] -= sqr(eval);
    eval2[ims] = sqrt(fabs(eval2[ims]));
  }

  return;
}

/*****************************************************************************/
