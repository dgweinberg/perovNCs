#include "fd.h"
#define calcBessel(x,x1)   ((x) < EPSR ? 0 : (sin((x)) * ((x1)*(x1)) - cos((x)) * (x1)))

int main(int argc, char const *argv[]){
	zomplex* psitot;	
	long_st ist;
	par_st par;
	fftw_plan_loc planfw, planbw; fftw_complex *fftwpsi;
	int jx,jy,jz, jyz, jxyz, i;
	double dx, dy, dz, *vx, *vy, *vz, prefactor, r, r_scaled, x,y,z;
	long flags = 0;
	FILE *pf;
	

//initialize grid and par and ist structs and fftw vars
	//ist struct
	ist.nx = ist.ny = ist.nz = 128;
	ist.ngrid = ist.nx * ist.ny * ist.nz;
	ist.nspin = 2;
	ist.nspinngrid = ist.nspin * ist.ngrid;
	ist.mstot = 1;
	ist.nthreads = 6;
	ist.nx_1 = 1.0/(double)(ist.nx);
	ist.ny_1 = 1.0/(double)(ist.ny);
	ist.nz_1 = 1.0/(double)(ist.ny);
	//par struct
	par.xmin = -5.00;
	par.xmax =  5.00;
	par.dx = (par.xmax-par.xmin)/((double)ist.nx);
	par.ymin = -5.00;
	par.ymax =  5.00;
	par.dy = (par.ymax-par.ymin)/((double)ist.ny);
	par.zmin = -5.00;
	par.zmax =  5.00;
	par.dz = (par.zmax-par.zmin)/((double)ist.nz);

	par.dv = par.dx*par.dy*par.dz;


	par.dkx = TWOPI / ((double)ist.nx * par.dx);
	par.dky = TWOPI / ((double)ist.ny * par.dy);
	par.dkz = TWOPI / ((double)ist.nz * par.dz);



	/*** the grid in the x, y, and z directions ***/
  	if ((vx = (double *) calloc(ist.nx, sizeof(double))) == NULL) nerror("vx");
  	if ((vy = (double *) calloc(ist.ny, sizeof(double))) == NULL) nerror("vy");
  	if ((vz = (double *) calloc(ist.nz, sizeof(double))) == NULL) nerror("vz");
  	
  	/***initializing the grid vectors  ***/
  	for (jx = 0, dx = par.xmin+0.5*par.dx; jx < ist.nx; jx++, dx += par.dx) vx[jx] = dx;
  	for (jy = 0, dy = par.ymin+0.5*par.dy; jy < ist.ny; jy++, dy += par.dy) vy[jy] = dy;
  	for (jz = 0, dz = par.zmin+0.5*par.dz; jz < ist.nz; jz++, dz += par.dz) vz[jz] = dz;
	
	/*** initialization for the fast Fourier transform ***/
  	fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid);
  	planfw = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, fftwpsi, fftwpsi, FFTW_FORWARD, flags);
  	planbw = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, fftwpsi, fftwpsi, FFTW_BACKWARD, flags);


	//write a known wavefunction to the grid
  	if ((psitot = (zomplex *) calloc(ist.nspinngrid*ist.mstot, sizeof(zomplex))) == NULL) nerror("psitot");
	omp_set_dynamic(0);
  	omp_set_num_threads(ist.nthreads);
  	#pragma omp parallel for private (jz,jy,jyz,jx,jxyz,r,x,y,z)
  	for (jz = 0; jz < ist.nz; jz++) { 
  		for (jy = 0;jy<ist.ny;jy++){
  			jyz = ist.nx * (ist.ny * jz + jy);
  			for(jx = 0; jx<ist.nx;jx++){
				jxyz = jyz + jx;
				x = vx[jx];
				y = vy[jy];
				z = vz[jz];
				r = sqrt((x)*(x)+(y)*(y)+(z)*(z))+ EPSR;
				for( i = 0;i<ist.mstot;i++){
					switch (i){

						case(1000):
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r));
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].im = 0.00;
						break;




						case(0):
						//|l=1 ml=+1,ms=+1/2> = |j = 3/2, jz = 3/2>
						// <J^2> = 15/4 = 3.75
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r)) * (x/r);
						psitot[jxyz+i*ist.nspinngrid].im = exp(-0.5*(r*r)) * (y/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].im = 0.00;
						break;

						case(1):
						//|l=1 ml=+1,ms=-1/2> = sqrt(1/3)|j = 3/2, jz = 1/2> + sqrt(2/3)|j = 1/2, jz = 1/2>
						//<J^2> = 7/4 = 1.75
						psitot[jxyz+i*ist.nspinngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].re = exp(-0.5*(r*r)) * (x/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].im = exp(-0.5*(r*r)) * (y/r);
						break;

						case(2):
						//|l=1 ml=0,ms=+1/2> = sqrt(2/3)|j = 3/2, jz = 1/2> - sqrt(1/3)|j = 1/2, jz = 1/2>
						//<J^2> = 11/4 = 2.75
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r)) * (z/r);
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].im = 0.00;
						break;

						case(3):
						//|l=1 ml=0,ms=-1/2> = sqrt(2/3)|j = 3/2, jz = -1/2> + sqrt(1/3)|j = 1/2, jz = -1/2>
						//<J^2> = 11/4 = 2.75
						psitot[jxyz+i*ist.nspinngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].re = exp(-0.5*(r*r)) * (z/r);
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].im = 0.00;

						case(4):
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r)) * (x/r);
						psitot[jxyz+i*ist.nspinngrid].im = -1.00 * exp(-0.5*(r*r)) * (y/r);
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].im = 0.00;
						break;


						case(5):
						psitot[jxyz+i*ist.nspinngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].re = exp(-0.5*(r*r)) * (x/r);
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].im = -1.00 * exp(-0.5*(r*r)) * (y/r);
						break;

						case(6):
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r));
						psitot[jxyz+i*ist.nspinngrid].im = exp(-0.5*(r*r));
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].im = 0.00;
						break;


						case(7):
						psitot[jxyz+i*ist.nspinngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].re = exp(-0.5*(r*r));
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].im = exp(-0.5*(r*r));
						break;


						case(8):
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r)) * (y/r);
						psitot[jxyz+i*ist.nspinngrid].im = exp(-0.5*(r*r)) * (z/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].im = 0.00;
						break;

						case(9):
						psitot[jxyz+i*ist.nspinngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].re = exp(-0.5*(r*r)) * (y/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].im = exp(-0.5*(r*r)) * (z/r);
						break;

						case(10):
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r));
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].re = exp(-0.5*(r*r));
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].im = 0.00;
						break;

						case(11):
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r));
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].re = -1.00*exp(-0.5*(r*r));
						psitot[jxyz+i*ist.nspinngrid + ist.ngrid].im = 0.00;
						break;
						
						case(12):
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r)) * (y/r);
						psitot[jxyz+i*ist.nspinngrid].im = exp(-0.5*(r*r)) * (z/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].re = -1.00 * exp(-0.5*(r*r)) * (z/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].im = exp(-0.5*(r*r)) * (y/r);
						break;

						case(13):
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r)) * (y/r);
						psitot[jxyz+i*ist.nspinngrid].im = exp(-0.5*(r*r)) * (z/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].re = exp(-0.5*(r*r)) * (z/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].im = -1.00*exp(-0.5*(r*r)) * (y/r);
						break;


						case(14):
						psitot[jxyz+i*ist.nspinngrid].re = exp(-0.5*(r*r)) * (z/r);
						psitot[jxyz+i*ist.nspinngrid].im = exp(-0.5*(r*r)) * (x/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].im = 0.00;
						break;

						case(15):
						psitot[jxyz+i*ist.nspinngrid].re = 0.00;
						psitot[jxyz+i*ist.nspinngrid].im = 0.00;
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].re = exp(-0.5*(r*r)) * (z/r);
						psitot[jxyz+i*ist.nspinngrid+ist.ngrid].im = exp(-0.5*(r*r)) * (x/r);
						break;
					}
					

				}

			}
		}
	}
	//fclose(pf);
	normalize_all( psitot , par.dv,ist.mstot, ist.nspinngrid, ist.nthreads);

	calcAngularExp(psitot, vx, vy, vz, planfw, planbw, fftwpsi, ist, par);
	
	free(vx); free(vy); free(vz);
	free(psitot);
	fftw_destroy_plan(planfw);
  	fftw_destroy_plan(planbw);
  	fftw_free(fftwpsi);


 return 0;
}
