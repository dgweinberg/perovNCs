#include "fd.h"
#include <float.h>

/***************************************************************************************/

void single_coulomb_openmp(zomplex       *psi, 
                           zomplex       *potq,
                           zomplex       *potqx,
                           zomplex       *poth,
                           double        *eval,
                           long_st        ist,
                           par_st        par,
                           fftw_plan_loc *planfw,
                           fftw_plan_loc *planbw,
                           fftw_complex  *fftwpsi,
                           zomplex       *bsmat,
                           zomplex       *direct,
                           zomplex       *exchange,
                           double       *h0mat)
{
    FILE   *pf;  
    long   flag, i, j, a, b, ibs, jbs, igrid, ispingrid; 
    int    tid, ispin; 
    long   *listibs;
    double ene, ene1, ene2;
    zomplex *rho, sum1, sum2, tmp;

    rho = (zomplex *) calloc(ist.ngrid * ist.nthreads, sizeof(zomplex));
    listibs = (long *) calloc(ist.ms2, sizeof(long));

    for (ibs = 0, a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
        for (i = 0; i < ist.totalhomo; i++, ibs++) {
            listibs[(a - ist.nlumo) * ist.totalhomo + i] = ibs;
        }
    }

    omp_set_dynamic(0);
    omp_set_num_threads(ist.nthreads);

    pf = fopen("direct.dat" , "w");
    /*** vabji direct ***/
    //loop over electron states a
    for (a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
        
        //loop over electron states b
#pragma omp parallel for private(sum1,ibs,jbs,ene1,ene2,ene,tid,igrid,ispin,ispingrid,b,i,j)
        for (b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {
            tid = omp_get_thread_num();
            
            //get joint density \rho_{ab}(r) = \sum_{\sigma} psi_{a}^{*}(r,\sigma) psi_{b}(r,\sigma)
            for (igrid = 0; igrid < ist.ngrid; igrid++) {
                rho[tid*ist.ngrid+igrid].re = rho[tid*ist.ngrid+igrid].im = 0.00;
                for (ispin = 0; ispin<2; ispin++){
                    ispingrid=igrid+ist.ngrid*ispin;
    	            rho[tid*ist.ngrid+igrid].re += psi[a*ist.nspinngrid+ispingrid].re * psi[b*ist.nspinngrid+ispingrid].re
                                                         + psi[a*ist.nspinngrid+ispingrid].im * psi[b*ist.nspinngrid+ispingrid].im;
    	            rho[tid*ist.ngrid+igrid].im += psi[a*ist.nspinngrid+ispingrid].re * psi[b*ist.nspinngrid+ispingrid].im
                                                         - psi[a*ist.nspinngrid+ispingrid].im * psi[b*ist.nspinngrid+ispingrid].re;
                }          
            }
            
            //this should populate the poth array with h_d(r) = \int W(r,r') \rho_{ab}(r') d^3r' via fourier transform
            hartree(&rho[tid*ist.ngrid], potqx, &poth[tid*ist.ngrid], ist, planfw[tid], planbw[tid], &fftwpsi[tid*ist.ngrid]);            
            
            //loop over hole states i
            for (i = 0; i < ist.totalhomo; i++) {
	            
                //loop over hole states j
                for (j = 0; j < ist.totalhomo; j++) {
	                //get pair state excitation energy
                    ene1 = eval[a] - eval[i];
	                ene2 = eval[b] - eval[j];
	                ene = ene1 - ene2;
                    
                    //integrate the effective potential to get K^d_{ai,bj}=\int h_d(r) \sum_\sigma psi_{i}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
                    sum1.re = sum1.im = 0.0;
					for (igrid = 0; igrid < ist.ngrid; igrid++) 
                        for (ispin = 0; ispin<2; ispin++){
                            ispingrid=igrid+ist.ngrid*ispin;
				            tmp.re = (psi[j*ist.nspinngrid + ispingrid].re * psi[i*ist.nspinngrid + ispingrid].re
                                     +psi[j*ist.nspinngrid + ispingrid].im * psi[i*ist.nspinngrid + ispingrid].im);
                            tmp.im = (psi[j*ist.nspinngrid + ispingrid].re * psi[i*ist.nspinngrid + ispingrid].im
                                     -psi[j*ist.nspinngrid + ispingrid].im * psi[i*ist.nspinngrid + ispingrid].re);
                            
                            sum1.re += poth[tid*ist.ngrid + igrid].re * tmp.re 
                                    -  poth[tid*ist.ngrid + igrid].im * tmp.im;
                            
                            sum1.im += poth[tid*ist.ngrid + igrid].re * tmp.im
                                    +  poth[tid*ist.ngrid + igrid].im * tmp.re;                                     
                        }
	                sum1.re *= par.dv;
                    sum1.im *= par.dv;

                    //get the matrix indicies for {ai,bj} and set bsmat
	                ibs = listibs[(a - ist.nlumo)*ist.totalhomo + i];
	                jbs = listibs[(b - ist.nlumo)*ist.totalhomo + j];
	                bsmat[ibs * ist.ms2 + jbs].re = sum1.re;
                    bsmat[ibs * ist.ms2 + jbs].im = sum1.im;
                    
                    direct[ibs * ist.ms2 + jbs].re = sum1.re;
                    direct[ibs * ist.ms2 + jbs].im = sum1.im;
	                
                    //if diagonal put the energy difference in the h0mat
                    if (ibs == jbs) 
                        h0mat[ibs*ist.ms2+jbs] = eval[a] - eval[i];
	                else 
                        h0mat[ibs*ist.ms2+jbs] = 0.0;
	                fprintf(pf,"%ld %ld %ld %ld %ld %ld %.*g %.*g %.*g %.*g\n",a,i,b,j,
		                     listibs[(a-ist.nlumo)*ist.totalhomo+i],
		                     listibs[(b-ist.nlumo)*ist.totalhomo+j],
		                     DBL_DIG, ene1, DBL_DIG, ene2, DBL_DIG, sum1.re, DBL_DIG, sum1.im);	            }
            }
        }
    }
    
    fclose(pf);
    /*** vjbai exchange ***/
    
    pf = fopen("exchange.dat" , "w");
    //loop over electron states a
    for (a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) {
#pragma omp parallel for private(sum2,ibs,jbs,ene1,ene2,ene,tid,igrid,b,i,j)
        
        //loop over hole states i
        for (i = 0; i < ist.totalhomo; i++) {
            tid = omp_get_thread_num();	
            ene1 = eval[a] - eval[i];
            
            //get joint density \rho_{ai}(r) = \sum_{\sigma} psi_{a}^*(r,\sigma) psi_{i}(r,\sigma)
            for (igrid = 0; igrid < ist.ngrid; igrid++) {
                rho[tid*ist.ngrid+igrid].re = rho[tid*ist.ngrid+igrid].im = 0.00;
                for (ispin = 0; ispin<2; ispin++){
                    ispingrid=igrid+ist.ngrid*ispin;
                    rho[tid*ist.ngrid+igrid].re += psi[a*ist.nspinngrid+ispingrid].re * psi[i*ist.nspinngrid+ispingrid].re
                                                 + psi[a*ist.nspinngrid+ispingrid].im * psi[i*ist.nspinngrid+ispingrid].im;
                    rho[tid*ist.ngrid+igrid].im += psi[a*ist.nspinngrid+ispingrid].re * psi[i*ist.nspinngrid+ispingrid].im
                                                 - psi[a*ist.nspinngrid+ispingrid].im * psi[i*ist.nspinngrid+ispingrid].re;
                }          
            }
            //this should populate the poth array with h_x(r) = \int v(r,r') \rho_{ai}(r') d^3r' but i don't understand how
            hartree(&rho[tid*ist.ngrid], potq, &poth[tid*ist.ngrid], ist, planfw[tid], planbw[tid], &fftwpsi[tid*ist.ngrid]);
            
            //loop over electron states b
            for (b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) {

                //loop over hole states j
	            for (j = 0; j < ist.totalhomo; j++) {
	                ene2 = eval[b] - eval[j];
	                ene = ene1 - ene2;
	                //TODO seems to be some bug with spin symmetry!
                    //integrate the effective potential to get K^x_{ai,bj}=\int h_x(r) \sum_\sigma psi_{b}(r,\sigma) psi_{j}^{*}(r,\sigma) d^3r
                    sum2.re = sum2.im = 0.0;
                    for ( igrid = 0; igrid < ist.ngrid; igrid++){
                        for (ispin = 0; ispin<2; ispin++){
                            ispingrid=igrid+ist.ngrid*ispin;
                            tmp.re = (psi[j*ist.nspinngrid + ispingrid].re * psi[b*ist.nspinngrid + ispingrid].re
                                     +psi[j*ist.nspinngrid + ispingrid].im * psi[b*ist.nspinngrid + ispingrid].im);
                            tmp.im = (psi[j*ist.nspinngrid + ispingrid].re * psi[b*ist.nspinngrid + ispingrid].im
                                     -psi[j*ist.nspinngrid + ispingrid].im * psi[b*ist.nspinngrid + ispingrid].re);

                            sum2.re += poth[tid*ist.ngrid + igrid].re * tmp.re
                                    -  poth[tid*ist.ngrid + igrid].im * tmp.im;
                            
                            sum2.im += poth[tid*ist.ngrid + igrid].re * tmp.im
                                    +  poth[tid*ist.ngrid + igrid].im * tmp.re;                                        
                        }
                    }
                    sum2.re *= par.dv;
                    sum2.im *= par.dv;

                    ibs = listibs[(a-ist.nlumo)*ist.totalhomo+i];
	                jbs = listibs[(b-ist.nlumo)*ist.totalhomo+j];
                    //NOTE: scalar version has a 2 as to calc for bright only. Don't want for full matrix. Took out 11/10 --DW
					bsmat[ibs*ist.ms2+jbs].re -=  sum2.re;
                    bsmat[ibs*ist.ms2+jbs].im -=  sum2.im;

                    exchange[ibs*ist.ms2+jbs].re = -1.0* sum2.re;
                    exchange[ibs*ist.ms2+jbs].im = -1.0* sum2.im;

                    if(ibs==jbs){bsmat[ibs*ist.ms2+jbs].im = 0.0;}
	                fprintf(pf,"%ld %ld %ld %ld %ld %ld %.*g %.*g %.*g %.*g\n",
                            a,i,b,j,ibs,jbs,
                            DBL_DIG, ene1, DBL_DIG, ene2, DBL_DIG, sum2.re, DBL_DIG, sum2.im);
	                fflush(0);
	            }
            }
        }
    }
    
	fclose(pf);
    
	free(rho); free(listibs);
	
    return;
}

/***************************************************************************************/
