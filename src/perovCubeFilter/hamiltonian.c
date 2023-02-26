/*****************************************************************************/
// File contains the main functions for calcuating the action of the hamiltonian on 
// a wavefunction including non-local spin-orbit contributions 

#include "fd.h"

/*****************************************************************************/

void hamiltonian_t(zomplex *phi, zomplex *psi, double *potl, double *ksqr, long_st ist, par_st par, nlc_st *nlc,
                long *nl, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,clock_t *timers ) {
  long i, itmp;
  clock_t start,diff;
  
  // Copy psi into phi
  memcpy(&psi[0], &phi[0], ist.nspinngrid*sizeof(psi[0]));
  

  start = clock();
  // Calculate the action of the kinetic energy part of the Hamiltonian on psi: |phi> = T|psi>
  kinetic(&phi[0*ist.ngrid], ksqr, planfw, planbw, fftwpsi, ist); //spin up
  kinetic(&phi[1*ist.ngrid], ksqr, planfw, planbw, fftwpsi, ist); //spin dn
  diff = clock()-start;
  timers[0]+=diff;

  //add the nonlocal potential into phi: |phi> += Vso|psi>

  if(ist.flagSO==1){
    start = clock();
    spin_orbit_proj_pot(phi, psi, ist, par, nlc, nl);
    diff = clock()-start;
    timers[1]+=diff;

    start = clock();
    nonlocal_proj_pot(phi, psi, ist, par, nlc, nl);
    diff = clock()-start;
    timers[2]+=diff;

  }

  
  start = clock();
  // Calculate the action of the local potential energy part of the Hamiltonian on psi
  // and add them all together: |phi> = T|psi> + Vlocal|psi> + Vso|psi>
  for (i = 0; i < ist.ngrid; i++) {
    phi[i].re += (potl[i] * psi[i].re);
    phi[i].im += (potl[i] * psi[i].im);

    itmp = i + ist.ngrid;
    phi[itmp].re += (potl[i] * psi[itmp].re);
    phi[itmp].im += (potl[i] * psi[itmp].im);
  }
  diff = clock()-start;
  timers[3]+=diff;
  


  return;
}

/*****************************************************************************/
void hamiltonian(zomplex *phi, zomplex *psi, double *potl, double *ksqr, long_st ist, par_st par, nlc_st *nlc,
                long *nl, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi) {
  long i, itmp;
  
  // Copy psi into phi
  memcpy(&psi[0], &phi[0], ist.nspinngrid*sizeof(psi[0]));
  

  
  // Calculate the action of the kinetic energy part of the Hamiltonian on psi: |phi> = T|psi>
  kinetic(&phi[0*ist.ngrid], ksqr, planfw, planbw, fftwpsi, ist); //spin up
  kinetic(&phi[1*ist.ngrid], ksqr, planfw, planbw, fftwpsi, ist); //spin dn
  

  //add the nonlocal potential into phi: |phi> += Vso|psi>

  if(ist.flagSO==1){
    
    spin_orbit_proj_pot(phi, psi, ist, par, nlc, nl);


    
    nonlocal_proj_pot(phi, psi, ist, par, nlc, nl);


  }

  
  
  // Calculate the action of the local potential energy part of the Hamiltonian on psi
  // and add them all together: |phi> = T|psi> + Vlocal|psi> + Vso|psi>
  for (i = 0; i < ist.ngrid; i++) {
    phi[i].re += (potl[i] * psi[i].re);
    phi[i].im += (potl[i] * psi[i].im);

    itmp = i + ist.ngrid;
    phi[itmp].re += (potl[i] * psi[itmp].re);
    phi[itmp].im += (potl[i] * psi[itmp].im);
  }



  return;
}

/*****************************************************************************/
// Calculates T|psi> via FFT and stores result in psi 

void kinetic(zomplex *psi, double *ksqr, fftw_plan_loc planfw, fftw_plan_loc planbw,
              fftw_complex *fftwpsi,long_st ist)
{
  long i;

  // Copy psi to fftwpsi
  memcpy(&fftwpsi[0], &psi[0], ist.ngrid*sizeof(fftwpsi[0]));
  
  // FT from r-space to k-space
  fftw_execute(planfw);
  
  // Kinetic energy is diagonal in k-space, just multiply fftwpsi by k^2
  for (i = 0; i < ist.ngrid; i++) {
    fftwpsi[i][0] *= ksqr[i];
    fftwpsi[i][1] *= ksqr[i];
  }

  // Inverse FT back to r-space
  fftw_execute(planbw);
  
  // Copy fftwpsi to psi to store T|psi> into |psi>
  memcpy(&psi[0], &fftwpsi[0], ist.ngrid*sizeof(psi[0]));

  return;
}

/*****************************************************************************/
// Calculates the action of the spin-orbit nonlocal potential using separable iprojs 
// and adds the result to phi: |phi> = |phi> + Vso|psi>
void spin_orbit_proj_pot(zomplex *phi, zomplex *psi, long_st ist, par_st par, nlc_st *nlc, long* nl ){
  zomplex Lx[3][3], Ly[3][3], Lz[3][3];
  zomplex Sx[2][2], Sy[2][2], Sz[2][2];
double sq2 = sqrt(0.5);

Lx[0][0].re=0.00; Lx[0][0].im=0.00;     Lx[0][1].re=sq2;  Lx[0][1].im=0.00;     Lx[0][2].re=0.00; Lx[0][2].im=0.00;
Lx[1][0].re=sq2;  Lx[1][0].im=0.00;     Lx[1][1].re=0.00; Lx[1][1].im=0.00;     Lx[1][2].re=sq2;  Lx[1][2].im=0.00;
Lx[2][0].re=0.00; Lx[2][0].im=0.00;     Lx[2][1].re=sq2;  Lx[2][1].im=0.00;     Lx[2][2].re=0.00; Lx[2][2].im=0.00;

Ly[0][0].re=0.00; Ly[0][0].im=0.00;     Ly[0][1].re=0.00; Ly[0][1].im=1.0*sq2;       Ly[0][2].re=0.00; Ly[0][2].im=0.00;
Ly[1][0].re=0.00; Ly[1][0].im=-1.0*sq2;  Ly[1][1].re=0.00; Ly[1][1].im=0.00;           Ly[1][2].re=0.00; Ly[1][2].im=1.0*sq2;
Ly[2][0].re=0.00; Ly[2][0].im=0.00;     Ly[2][1].re=0.00; Ly[2][1].im=-1.0*sq2;        Ly[2][2].re=0.00; Ly[2][2].im=0.00;
//Define Lz 
Lz[0][0].re=-1.00; Lz[0][0].im=0.00;     Lz[0][1].re=0.00; Lz[0][1].im=0.00;     Lz[0][2].re=0.00;  Lz[0][2].im=0.00;
Lz[1][0].re=0.00; Lz[1][0].im=0.00;     Lz[1][1].re=0.00; Lz[1][1].im=0.00;     Lz[1][2].re=0.00;  Lz[1][2].im=0.00;
Lz[2][0].re=0.00; Lz[2][0].im=0.00;     Lz[2][1].re=0.00; Lz[2][1].im=0.00;     Lz[2][2].re=1.00; Lz[2][2].im=0.00;

//Define Sx
Sx[0][0].re=0.00; Sx[0][0].im=0.00;     Sx[0][1].re=0.50; Sx[0][1].im=0.00;
Sx[1][0].re=0.50; Sx[1][0].im=0.00;     Sx[1][1].re=0.00; Sx[1][1].im=0.00;
//Define Sy
Sy[0][0].re=0.00; Sy[0][0].im=0.00;     Sy[0][1].re=0.00; Sy[0][1].im=-0.50;
Sy[1][0].re=0.00; Sy[1][0].im=0.50;     Sy[1][1].re=0.00; Sy[1][1].im=0.00;
//Define Sz
Sz[0][0].re=0.50; Sz[0][0].im=0.00;     Sz[0][1].re=0.00;  Sz[0][1].im=0.00;
Sz[1][0].re=0.00; Sz[1][0].im=0.00;     Sz[1][1].re=-0.50; Sz[1][1].im=0.00;


  


  zomplex proj;
  long jatom, j1, index1, index1_psi;
  int iproj, spin_p, m_p, spin, m;
  for ( jatom =0; jatom<ist.nnonlocal; jatom++){
    for ( iproj = 0; iproj<ist.nproj; iproj++){
      for ( spin_p = 0; spin_p<2; spin_p++){
        for ( m_p = 0; m_p < 3; m_p++){
          proj.re = proj.im = 0.00;
          for ( j1=0; j1<nl[jatom]; j1++){
            index1 = jatom * ist.nnlc + j1;
            index1_psi = nlc[index1].jxyz + (ist.ngrid)*spin_p;
            
            //weird signs b/c of Y_{lm}^*            
            proj.re += psi[index1_psi].re * nlc[index1].y1[m_p].re * nlc[index1].proj[iproj];
            proj.re += psi[index1_psi].im * nlc[index1].y1[m_p].im * nlc[index1].proj[iproj];
           
            proj.im += psi[index1_psi].im * nlc[index1].y1[m_p].re * nlc[index1].proj[iproj];
            proj.im -= psi[index1_psi].re * nlc[index1].y1[m_p].im * nlc[index1].proj[iproj];
            
          }
          proj.re *= par.dv;
          proj.im *= par.dv;   
          for (spin = 0; spin<2; spin++){
            for (m = 0; m<3; m++){
              //get L_{m,m'}\cdot S_{s,s'}*P_{n,m',s'} = PLS_{n,m,m',s,s'}
              zomplex LS, PLS;

              LS.re = LS.im = 0.00;
              PLS.re = PLS.im = 0.00;
              //checked ordering of m and m_p, printed LdotS and checked against ref. 
              LS.re += (  Lx[m][m_p].re * Sx[spin][spin_p].re - Lx[m][m_p].im * Sx[spin][spin_p].im);
              LS.im += (  Lx[m][m_p].re * Sx[spin][spin_p].im + Lx[m][m_p].im * Sx[spin][spin_p].re);

              LS.re += (  Ly[m][m_p].re * Sy[spin][spin_p].re - Ly[m][m_p].im * Sy[spin][spin_p].im);
              LS.im += (  Ly[m][m_p].re * Sy[spin][spin_p].im + Ly[m][m_p].im * Sy[spin][spin_p].re);

              LS.re += (  Lz[m][m_p].re * Sz[spin][spin_p].re - Lz[m][m_p].im * Sz[spin][spin_p].im);
              LS.im += (  Lz[m][m_p].re * Sz[spin][spin_p].im + Lz[m][m_p].im * Sz[spin][spin_p].re);              
              //printf("m:%i m':%i s:%i s':%i  LS: %f+i*%f\n", m, m_p, spin, spin_p, LS.re, LS.im);

              PLS.re = LS.re * proj.re - LS.im * proj.im;
              PLS.im = LS.re * proj.im + LS.im * proj.re;

              for (j1=0; j1<nl[jatom]; j1++){
                index1 = jatom * ist.nnlc + j1;
                index1_psi = nlc[index1].jxyz + (ist.ngrid)*spin;

                phi[index1_psi].re += nlc[index1].proj[iproj] * nlc[index1].y1[m].re * PLS.re;
                phi[index1_psi].re -= nlc[index1].proj[iproj] * nlc[index1].y1[m].im * PLS.im;

                phi[index1_psi].im += nlc[index1].proj[iproj] * nlc[index1].y1[m].re * PLS.im;
                phi[index1_psi].im += nlc[index1].proj[iproj] * nlc[index1].y1[m].im * PLS.re;
                
              }
            } 
          }
        }
      }
    }

  }
}
/*****************************************************************************/

// Calculates the action of the spin-orbit nonlocal potential using separable iprojs 
// and adds the result to phi: |phi> = |phi> + Vso|psi>
void nonlocal_proj_pot(zomplex *phi, zomplex *psi, long_st ist, par_st par, nlc_st *nlc, long* nl ){
  zomplex proj;
  long jatom, j1, index1, index1_psi;
  int iproj, spin, m;
  for ( jatom =0; jatom<ist.nnonlocal; jatom++){
    for ( iproj = 0; iproj<ist.nproj; iproj++){
      for ( spin = 0; spin<2; spin++){
        for ( m = 0; m < 3; m++){
          proj.re = proj.im = 0.00;
          for ( j1=0; j1<nl[jatom]; j1++){
            index1 = jatom * ist.nnlc + j1;
            index1_psi = nlc[index1].jxyz + (ist.ngrid)*spin;
            
            //weird signs b/c of Y_{lm}^*            
            proj.re += psi[index1_psi].re * nlc[index1].y1[m].re * nlc[index1].nlProj[iproj];
            proj.re += psi[index1_psi].im * nlc[index1].y1[m].im * nlc[index1].nlProj[iproj];
           
            proj.im += psi[index1_psi].im * nlc[index1].y1[m].re * nlc[index1].nlProj[iproj];
            proj.im -= psi[index1_psi].re * nlc[index1].y1[m].im * nlc[index1].nlProj[iproj];
            
          }
          proj.re *= par.dv*nlc[index1].nlProjSgn[iproj];
          proj.im *= par.dv*nlc[index1].nlProjSgn[iproj];   
          
          for (j1=0; j1<nl[jatom]; j1++){
            index1 = jatom * ist.nnlc + j1;
            index1_psi = nlc[index1].jxyz + (ist.ngrid)*spin;

            phi[index1_psi].re += nlc[index1].nlProj[iproj] * nlc[index1].y1[m].re * proj.re;
            phi[index1_psi].re -= nlc[index1].nlProj[iproj] * nlc[index1].y1[m].im * proj.im;

            phi[index1_psi].im += nlc[index1].nlProj[iproj] * nlc[index1].y1[m].re * proj.im;
            phi[index1_psi].im += nlc[index1].nlProj[iproj] * nlc[index1].y1[m].im * proj.re;
            

          }
        }
      }
    }

  }
}
/*****************************************************************************/


void timeReverseAll(zomplex *psitot, zomplex *dest,long_st ist){
  long ms, ims,i;
  long ngrid = ist.ngrid;
  omp_set_dynamic(0);
  omp_set_num_threads(ist.nthreads);
#pragma omp parallel for private(ms, ims,i)
  for (ms = 0; ms < ist.mstot; ms++) {
    for (ims = ngrid*ms,i=0; i<ngrid;i++){
      
      //dest dn = (psi up)^*
      dest[ims+i+ngrid].re = psitot[ims+i].re;
      dest[ims+i+ngrid].im = -1.0*psitot[ims+i].im;

      //dest up = -(psi dn)^*
      dest[ims+i].re = -1.0*psitot[ims+i+ngrid].re;
      dest[ims+i].im = psitot[ims+i+ngrid].im;

    }
       
  } 

}

