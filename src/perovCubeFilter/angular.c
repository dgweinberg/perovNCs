#include "fd.h"

/************************************************************/
void calcAngularExp(zomplex* psitot, double* vx, double* vy, double* vz,
 	fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,long_st ist, par_st par, int start, int stop){
	int i, ipsi; long jx,jy,jz,jyz,jxyz;
	double Jxexp, Jxsqr, Jyexp, Jysqr,Jzexp, Jzsqr, Jsqrexp, norm;
	zomplex *psi, *Jxpsi, *Jypsi, *Jzpsi;
	FILE* pf = fopen("angular.dat", "w");
	fprintf(pf, "###\tJx\t\tJx^2\t\tJy\t\tJy^2\t\tJz\t\tJz^2\t\tJ^2\t\tnorm\n");
	
	//allocate memory for the operated wavefucntions
	if ((Jxpsi  = (zomplex*)calloc(ist.nspinngrid,sizeof(zomplex)))==NULL)nerror("Jxpsi");
  	if ((Jypsi  = (zomplex*)calloc(ist.nspinngrid,sizeof(zomplex)))==NULL)nerror("Jypsi");
  	if ((Jzpsi  = (zomplex*)calloc(ist.nspinngrid,sizeof(zomplex)))==NULL)nerror("Jzpsi");
  	if ((psi = (zomplex*)calloc(ist.nspinngrid,sizeof(zomplex)))==NULL)nerror("psi");

	for(ipsi = start; ipsi < stop; ipsi++){	  	
	  	//psi = &psitot[ipsi*ist.nspinngrid];
	  	for (jz = 0; jz < ist.nz; jz++) { 
	  		for (jy = 0;jy<ist.ny;jy++){
	  			jyz = ist.nx * (ist.ny * jz + jy);
	  			for(jx = 0; jx<ist.nx;jx++){
					jxyz = jyz + jx;
					psi[jxyz].re = psitot[ipsi*ist.nspinngrid+jxyz].re;
					psi[jxyz].im = psitot[ipsi*ist.nspinngrid+jxyz].im;
					psi[jxyz+ist.ngrid].re = psitot[ipsi*ist.nspinngrid+jxyz+ist.ngrid].re;
					psi[jxyz+ist.ngrid].im = psitot[ipsi*ist.nspinngrid+jxyz+ist.ngrid].im;
				}
			}		
		}

		
		lowPassFilter(&psi[0],planfw,planbw,fftwpsi, ist, par);
		lowPassFilter(&psi[ist.ngrid],planfw,planbw,fftwpsi, ist, par);
		
		normalize(psi,par.dv,ist.nspinngrid,ist.nthreads);

		char filename[20];
		double* rho = calloc(ist.ngrid,sizeof(double));
		
		for (long igrid = 0;igrid<ist.ngrid; igrid++){
    		rho[igrid] = (psi[igrid].re);
		 }
		sprintf(filename, "smpsi%iUpRe.cube", ipsi);
	  	writeCubeFile(rho, par,ist, filename);

	  	for (long igrid = 0;igrid<ist.ngrid; igrid++){
    		rho[igrid] = (psi[igrid].im);
	  	}
	  	sprintf(filename, "smpsi%iUpIm.cube", ipsi);
	  	writeCubeFile(rho, par,ist, filename);

	  	for (long igrid = 0;igrid<ist.ngrid; igrid++){
   	 		rho[igrid] = (psi[ist.ngrid+igrid].re);
	  	}
	  	sprintf(filename, "smpsi%iDnRe.cube", ipsi);
	  	writeCubeFile(rho, par,ist, filename);

	  	for (long igrid = 0;igrid<ist.ngrid; igrid++){
	    	rho[igrid] = (psi[ist.ngrid+igrid].im);
	  	}
 	 	sprintf(filename, "smpsi%iDnIm.cube", ipsi);
	  	writeCubeFile(rho, par,ist, filename);
		  

		

		//compute the action of the J operator
	  	jOpp( Jxpsi,  Jypsi,  Jzpsi, psi, vx,  vy,  vz,	 planfw,  planbw, fftwpsi, ist,  par);

	  	//lOpp( &Jxpsi[0], &Jypsi[0], &Jzpsi[0], &psi[0], vx,  vy,  vz, planfw, planbw, fftwpsi, ist, par);
	  	//lOpp( &Jxpsi[ist.ngrid], &Jypsi[ist.ngrid], &Jzpsi[ist.ngrid], &psi[ist.ngrid], vx,  vy,  vz, planfw, planbw, fftwpsi, ist, par);

	  	//compute Ji expectation as <psi|Jipsi> and J^2 as \sum_i <Jipsi|Jipsi>
	  	Jxexp = Jxsqr = Jyexp = Jysqr = Jzexp = Jzsqr = Jsqrexp = norm = 0.00;
	  	omp_set_dynamic(0);
	  	omp_set_num_threads(ist.nthreads);
		#pragma omp parallel for reduction (+:Jxexp, Jxsqr, Jyexp, Jysqr,Jzexp, Jzsqr, Jsqrexp, norm)
		for (i = 0; i < ist.nspinngrid; i++){
			//<psi|Jxpsi>
			Jxexp += (psi[i].re * Jxpsi[i].re) + (psi[i].im * Jxpsi[i].im);
			Jxsqr += (Jxpsi[i].re * Jxpsi[i].re) + (Jxpsi[i].im * Jxpsi[i].im);

			//<psi|Jypsi>
			Jyexp += (psi[i].re * Jypsi[i].re) + (psi[i].im * Jypsi[i].im);
			Jysqr += (Jypsi[i].re * Jypsi[i].re) + (Jypsi[i].im * Jypsi[i].im);

			//<psi|Jzpsi>
			Jzexp += (psi[i].re * Jzpsi[i].re) + (psi[i].im * Jzpsi[i].im);	
			Jzsqr += (Jzpsi[i].re * Jzpsi[i].re) + (Jzpsi[i].im * Jzpsi[i].im);	

			Jsqrexp += (Jxpsi[i].re*Jxpsi[i].re) + (Jxpsi[i].im*Jxpsi[i].im) 
					+ (Jypsi[i].re*Jypsi[i].re) + (Jypsi[i].im*Jypsi[i].im)
					+ (Jzpsi[i].re*Jzpsi[i].re) + (Jzpsi[i].im*Jzpsi[i].im);
		
			norm += (psi[i].re*psi[i].re)+(psi[i].im*psi[i].im);
		}
		Jxexp*=par.dv;
		Jxsqr*=par.dv;
		Jyexp*=par.dv;
		Jysqr*=par.dv;
		Jzexp*=par.dv;
		Jzsqr*=par.dv;
		Jsqrexp*=par.dv;
		norm*=par.dv;


		fprintf(pf, "%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
					ipsi,Jxexp, Jxsqr, Jyexp, Jysqr, Jzexp, Jzsqr, Jsqrexp, norm);
	}

	fclose(pf);
	free(Jxpsi); free(Jypsi); free(Jzpsi);


}

/************************************************************/
//Calculate the three vector components of the action of the J=L+S operator
void jOpp( zomplex* Jxpsi, zomplex* Jypsi, zomplex* Jzpsi,zomplex* psi, 
	double* vx, double* vy, double* vz,
 	fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,long_st ist, par_st par){
	long jz,jy,jyz,jx,jxyz;
	double density,x,y,z;

	/*** Find electron center of mass ***/
	x=y=z=0;
	omp_set_dynamic(0);
  	omp_set_num_threads(ist.nthreads);
  	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz) reduction(+:x,y,z)
  	for (jz = 0; jz < ist.nz; jz++) { 
  		for (jy = 0;jy<ist.ny;jy++){
  			jyz = ist.nx * (ist.ny * jz + jy);
  			for(jx = 0; jx<ist.nx;jx++){
				jxyz = jyz + jx;
				density = psi[jxyz].re *psi[jxyz].re + psi[jxyz].im * psi[jxyz].im 
				+ psi[jxyz+ist.ngrid].re *psi[jxyz+ist.ngrid].re + psi[jxyz+ist.ngrid].im * psi[jxyz+ist.ngrid].im;
				
				x += vx[jx] * density;
				y += vy[jy] * density;
				z += vz[jz] * density;
			}
		}
	}
	/*
	center.x = x*par.dv; center.y=y*par.dv;center.z=z*par.dv;
	printf("Center: <%f, %f, %f>\n", center.x, center.y , center.z);
	*/
	
	//spin up part
	lOpp(&Jxpsi[0],&Jypsi[0],&Jzpsi[0], &psi[0], 
		vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
	

	//spin dn part
	lOpp(&Jxpsi[ist.ngrid],&Jypsi[ist.ngrid],&Jzpsi[ist.ngrid], &psi[ist.ngrid], 
		vx,vy,vz,planfw,planbw, fftwpsi,ist, par);


	//add the spin part which acts locally
	omp_set_dynamic(0);
  	omp_set_num_threads(ist.nthreads);
	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < ist.nz; jz++) { 
  		for (jy = 0;jy<ist.ny;jy++){
  			jyz = ist.nx * (ist.ny * jz + jy);
  			for(jx = 0; jx<ist.nx;jx++){
				jxyz = jyz + jx;
				//Jx|r,up> = Lx|r,up> + 1/2 |r, dn>
				Jxpsi[jxyz].re += 0.5*psi[jxyz+ist.ngrid].re;
				Jxpsi[jxyz].im += 0.5*psi[jxyz+ist.ngrid].im;
				//Jx|r,dn> = Lx|r,dn> + 1/2 |r, up>
				Jxpsi[jxyz+ist.ngrid].re += 0.5*psi[jxyz].re;
				Jxpsi[jxyz+ist.ngrid].im += 0.5*psi[jxyz].im;

				//Jy|r,up> = Ly|r,up> - i/2 |r, dn>
				Jypsi[jxyz].re += 0.5* psi[jxyz+ist.ngrid].im;
				Jypsi[jxyz].im -= 0.5* psi[jxyz+ist.ngrid].re;
				//Jy|r,dn> = Ly|r,dn> + i/2 |r, up>
				Jypsi[jxyz+ist.ngrid].re -= 0.5* psi[jxyz].im;
				Jypsi[jxyz+ist.ngrid].im += 0.5* psi[jxyz].re;

				//Jz|r,up> = Lz|r,up> + 1/2 |r,up>
				Jzpsi[jxyz].re += 0.5*psi[jxyz].re;
				Jzpsi[jxyz].im += 0.5*psi[jxyz].im;
				//Jz|r,dn> = Lz|r,dn> - 1/2 |r,dn>
				Jzpsi[jxyz+ist.ngrid].re -= 0.5*psi[jxyz+ist.ngrid].re;
				Jzpsi[jxyz+ist.ngrid].im -= 0.5*psi[jxyz+ist.ngrid].im;

			}
		}
	}

}


/************************************************************/
//Calculate the three vector components of the action of the L operator on the 
//spatial part of the grid (no spin part)
void lOpp(zomplex* Lxpsi, zomplex* Lypsi, zomplex* Lzpsi, zomplex* psi, 
	double* vx, double* vy, double* vz,
 	fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,long_st ist, par_st par){
	
	double *kx, *ky, *kz,x, y, z;
	long jx,jy,jz, jyz, jxyz;

	xyz_st center;
	center.x = center.y = center.z = 0.00;


	/*** First use the fft to get the action of the p operator on each axis (-i d/dx) imag part from the fft definition***/


	//generate kx,ky,kz vectors
	if ((kx  = (double*)calloc(ist.nx,sizeof(double)))==NULL)nerror("kx");
  	if ((ky  = (double*)calloc(ist.ny,sizeof(double)))==NULL)nerror("ky");
  	if ((kz  = (double*)calloc(ist.nz,sizeof(double)))==NULL)nerror("kz");

  	/***initializing the k vectors ***/
  	//negative frequencies count backwards from end
  	//overall transform fw and bw aquires factor of nx*ny*nz so normalize here
  	for (kx[0] = 0.0, jx = 1; jx <= ist.nx / 2; jx++){
    	kx[ist.nx-jx] = -1.00 * (kx[jx] = (double)(jx) * par.dkx * 
    		ist.nx_1 * ist.ny_1 * ist.nz_1);
  	}

  	for (ky[0] = 0.0, jy = 1; jy <= ist.ny / 2; jy++){
    	ky[ist.ny-jy] = -1.00 * (ky[jy] = (double)(jy) * par.dky *
			ist.nx_1 * ist.ny_1 * ist.nz_1);
  	}

  	for (kz[0] = 0.0, jz = 1; jz <= ist.nz / 2; jz++){
    	kz[ist.nz-jz] = -1.00 * (kz[jz] = (double)(jz) * par.dkz *
			ist.nx_1 * ist.ny_1 * ist.nz_1);
  	}


	// Copy psi to fftwpsi
  	memcpy(&fftwpsi[0], &psi[0], ist.ngrid*sizeof(fftwpsi[0]));
  	// FT from r-space to k-space
  	fftw_execute(planfw);
  	
  	omp_set_dynamic(0);
  	omp_set_num_threads(ist.nthreads);
  	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < ist.nz; jz++) { 
  		for (jy = 0;jy<ist.ny;jy++){
  			jyz = ist.nx * (ist.ny * jz + jy);
  			for(jx = 0; jx<ist.nx;jx++){
				jxyz = jyz + jx;
				//multiply by kx to get partial along x-axis
				fftwpsi[jxyz][0] *=kx[jx];
				fftwpsi[jxyz][1] *=kx[jx];
			}
		}
	}

	// Inverse FT back to r-space
	fftw_execute(planbw);
	// Copy fftwpsi to psi to store Lx|psi> into |Lxpsi>
	memcpy(&Lxpsi[0], &fftwpsi[0], ist.ngrid*sizeof(Lxpsi[0]));
	

	// Copy psi to fftwpsi
  	memcpy(&fftwpsi[0], &psi[0], ist.ngrid*sizeof(fftwpsi[0]));

  	// FT from r-space to k-space
  	fftw_execute(planfw);
  	
  	omp_set_dynamic(0);
  	omp_set_num_threads(ist.nthreads);
 	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < ist.nz; jz++) { 
  		for (jy = 0;jy<ist.ny;jy++){
  			jyz = ist.nx * (ist.ny * jz + jy);
  			for(jx = 0; jx<ist.nx;jx++){
				jxyz = jyz + jx;
				//multiply by ky to get partial along y-axis
				fftwpsi[jxyz][0] *= ky[jy];
				fftwpsi[jxyz][1] *= ky[jy];
			}
		}
	}

	// Inverse FT back to r-space
	fftw_execute(planbw);
  
	// Copy fftwpsi to psi to store Ly|psi> into |Lypsi>
	memcpy(&Lypsi[0], &fftwpsi[0], ist.ngrid*sizeof(Lypsi[0]));
	
	

	// Copy psi to fftwpsi
  	memcpy(&fftwpsi[0], &psi[0], ist.ngrid*sizeof(fftwpsi[0]));

  	// FT from r-space to k-space
  	fftw_execute(planfw);
  	omp_set_dynamic(0);
  	omp_set_num_threads(ist.nthreads);
  	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < ist.nz; jz++) { 
  		for (jy = 0;jy<ist.ny;jy++){
  			jyz = ist.nx * (ist.ny * jz + jy);
  			for(jx = 0; jx<ist.nx;jx++){
				jxyz = jyz + jx;
				//multiply by kz to get partial along z-axis
				fftwpsi[jxyz][0] *= kz[jz];
				fftwpsi[jxyz][1] *= kz[jz];
			}
		}
	}

	// Inverse FT back to r-space
	fftw_execute(planbw);
  
	// Copy fftwpsi to psi to store Lz|psi> into |Lzpsi>
	memcpy(&Lzpsi[0], &fftwpsi[0], ist.ngrid*sizeof(Lzpsi[0]));

	free(kx); free(ky); free(kz);
	


	/*** Now do the cross product part at each grid point***/
	zomplex gradx,grady,gradz;
	omp_set_dynamic(0);
  	omp_set_num_threads(ist.nthreads);
	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz,gradx,grady,gradz,x,y,z)
  	for (jz = 0; jz < ist.nz; jz++) { 
  		z = vz[jz]-center.z;
  		for (jy = 0;jy<ist.ny;jy++){
  			jyz = ist.nx * (ist.ny * jz + jy);
  			y = vy[jy]- center.y;
  			for(jx = 0; jx<ist.nx;jx++){
  				x = vx[jx]-center.x;
				jxyz = jyz + jx;
				//copy over tmp varibales
				gradx.re = Lxpsi[jxyz].re; gradx.im = Lxpsi[jxyz].im;
				grady.re = Lypsi[jxyz].re; grady.im = Lypsi[jxyz].im;
				gradz.re = Lzpsi[jxyz].re; gradz.im = Lzpsi[jxyz].im;
				
				Lxpsi[jxyz].re = (y*gradz.re - z*grady.re);

				Lxpsi[jxyz].im = (y*gradz.im - z*grady.im);

				Lypsi[jxyz].re = (z*gradx.re - x*gradz.re);
		
				Lypsi[jxyz].im = (z*gradx.im - x*gradz.im);

				Lzpsi[jxyz].re = (x*grady.re - y*gradx.re);

				Lzpsi[jxyz].im = (x*grady.im - y*gradx.im);


			}
		}
	}




}
/************************************************************/
void lowPassFilter(zomplex* psi,fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,
	long_st ist, par_st par){
	double *kx,*ky,*kz;
	long jx,jy,jz,jyz,jxyz;
	//generate kx,ky,kz vectors
	if ((kx  = (double*)calloc(ist.nx,sizeof(double)))==NULL)nerror("kx");
  	if ((ky  = (double*)calloc(ist.ny,sizeof(double)))==NULL)nerror("ky");
  	if ((kz  = (double*)calloc(ist.nz,sizeof(double)))==NULL)nerror("kz");

  	/***initializing the k vectors ***/
  	//negative frequencies count backwards from end
  	//overall transform fw and bw aquires factor of nx*ny*nz so normalize here
  	for (kx[0] = 0.0, jx = 1; jx <= ist.nx / 2; jx++){
    	kx[ist.nx-jx] = -1.00 * (kx[jx] = (double)(jx) * par.dkx );
  	}

  	for (ky[0] = 0.0, jy = 1; jy <= ist.ny / 2; jy++){
    	ky[ist.ny-jy] = -1.00 * (ky[jy] = (double)(jy) * par.dky );
  	}

  	for (kz[0] = 0.0, jz = 1; jz <= ist.nz / 2; jz++){
    	kz[ist.nz-jz] = -1.00 * (kz[jz] = (double)(jz) * par.dkz );
  	}


	// Copy psi to fftwpsi
  	memcpy(&fftwpsi[0], &psi[0], ist.ngrid*sizeof(fftwpsi[0]));
  	// FT from r-space to k-space
  	fftw_execute(planfw);
  	
  	omp_set_dynamic(0);
  	omp_set_num_threads(ist.nthreads);
  	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz)
  	for (jz = 0; jz < ist.nz; jz++) { 
  		for (jy = 0;jy<ist.ny;jy++){
  			jyz = ist.nx * (ist.ny * jz + jy);
  			for(jx = 0; jx<ist.nx;jx++){
				jxyz = jyz + jx;
				
				if(sqr(kx[jx])+sqr(ky[jy])+sqr(kz[jz]) < sqr(3.0*par.dkz)){
					fftwpsi[jxyz][0] *=1.00;
					fftwpsi[jxyz][1] *=1.00;
				}
				else{
					fftwpsi[jxyz][0] =0;
					fftwpsi[jxyz][1] =0;
				}
			}
		}
	}

	// Inverse FT back to r-space
	fftw_execute(planbw);
	// Copy fftwpsi to psi to store Lx|psi> into |Lxpsi>
	memcpy(&psi[0], &fftwpsi[0], ist.ngrid*sizeof(psi[0]));
	
}

/************************************************************/
