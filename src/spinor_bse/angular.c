/****************************************************************************/

#include "fd.h"

/****************************************************************************/

void spins(zomplex *sx, zomplex *sy, zomplex *sz,zomplex *psi,long_st ist,par_st par){
long a,i,b,j,jx, jy, jz, jgridup,jgriddn, jyz, index;
FILE *pfx, *pfy, *pfz;

pfx = fopen("sx.dat", "w"); pfy = fopen("sy.dat", "w"); pfz = fopen("sz.dat", "w");
//calculate spin matrix elements between all occupied (hole) orbitals
	for(i =0; i< ist.totalhomo;i++){
		for (j = 0; j< ist.totalhomo;j++){
			
			for (jz = 0; jz < ist.nz; jz++) {
        		for (jy = 0; jy < ist.ny; jy++) {
          			jyz = ist.nx * (ist.ny * jz + jy);
          			for (jx = 0; jx < ist.nx; jx++) {
          				jgridup = jyz + jx;
            			jgriddn = jgridup+ist.ngrid;
						
					//Spin x part
						//<j|r,dn> * <r,up|i>
						sx[i*ist.totalhomo+j].re += psi[j*ist.nspinngrid+jgriddn].re * psi[i*ist.nspinngrid+jgridup].re
												  + psi[j*ist.nspinngrid+jgriddn].im * psi[i*ist.nspinngrid+jgridup].im;
						sx[i*ist.totalhomo+j].im += psi[j*ist.nspinngrid+jgriddn].re * psi[i*ist.nspinngrid+jgridup].im
												  - psi[j*ist.nspinngrid+jgriddn].im * psi[i*ist.nspinngrid+jgridup].re;
						//<j|r,up> * <r,dn|i>
						sx[i*ist.totalhomo+j].re += psi[j*ist.nspinngrid+jgridup].re * psi[i*ist.nspinngrid+jgriddn].re
												  + psi[j*ist.nspinngrid+jgridup].im * psi[i*ist.nspinngrid+jgriddn].im;
						sx[i*ist.totalhomo+j].im += psi[j*ist.nspinngrid+jgridup].re * psi[i*ist.nspinngrid+jgriddn].im
												  - psi[j*ist.nspinngrid+jgridup].im * psi[i*ist.nspinngrid+jgriddn].re;
					//Spin y part
						//i*<j|r,dn> * <r,up|i>
						sy[i*ist.totalhomo+j].im += psi[j*ist.nspinngrid+jgriddn].re * psi[i*ist.nspinngrid+jgridup].re
												  + psi[j*ist.nspinngrid+jgriddn].im * psi[i*ist.nspinngrid+jgridup].im;
						sy[i*ist.totalhomo+j].re -= psi[j*ist.nspinngrid+jgriddn].re * psi[i*ist.nspinngrid+jgridup].im
												  - psi[j*ist.nspinngrid+jgriddn].im * psi[i*ist.nspinngrid+jgridup].re;
						//-i*<j|r,up> * <r,dn|i>
						sy[i*ist.totalhomo+j].im -= psi[j*ist.nspinngrid+jgridup].re * psi[i*ist.nspinngrid+jgriddn].re
												  + psi[j*ist.nspinngrid+jgridup].im * psi[i*ist.nspinngrid+jgriddn].im;
						sy[i*ist.totalhomo+j].re += psi[j*ist.nspinngrid+jgridup].re * psi[i*ist.nspinngrid+jgriddn].im
												  - psi[j*ist.nspinngrid+jgridup].im * psi[i*ist.nspinngrid+jgriddn].re;
					//Spin z part
						//<j|r,up> * <r,up|i>
						sz[i*ist.totalhomo+j].re += psi[j*ist.nspinngrid+jgridup].re * psi[i*ist.nspinngrid+jgridup].re
												  + psi[j*ist.nspinngrid+jgridup].im * psi[i*ist.nspinngrid+jgridup].im;
						sz[i*ist.totalhomo+j].im += psi[j*ist.nspinngrid+jgridup].re * psi[i*ist.nspinngrid+jgridup].im
												  - psi[j*ist.nspinngrid+jgridup].im * psi[i*ist.nspinngrid+jgridup].re;
						//<j|r,dn> * <r,dn|i>
						sz[i*ist.totalhomo+j].re -= psi[j*ist.nspinngrid+jgriddn].re * psi[i*ist.nspinngrid+jgriddn].re
												  + psi[j*ist.nspinngrid+jgriddn].im * psi[i*ist.nspinngrid+jgriddn].im;
						sz[i*ist.totalhomo+j].im -= psi[j*ist.nspinngrid+jgriddn].re * psi[i*ist.nspinngrid+jgriddn].im
												  - psi[j*ist.nspinngrid+jgriddn].im * psi[i*ist.nspinngrid+jgriddn].re;
				

					}
				}
			}
			//multiply all by (1/2)*dV and complex conjugate
			sx[i*ist.totalhomo+j].re*=-0.5*par.dv;
			sx[i*ist.totalhomo+j].im*=0.5*par.dv;
			sy[i*ist.totalhomo+j].re*=-0.5*par.dv;
			sy[i*ist.totalhomo+j].im*=0.5*par.dv;
			sz[i*ist.totalhomo+j].re*=-0.5*par.dv;
			sz[i*ist.totalhomo+j].im*=0.5*par.dv;

			fprintf (pfx,"%ld %ld %g %g\n",i,j,sx[i*ist.totalhomo+j].re, sx[i*ist.totalhomo+j].im);
      		fprintf (pfy,"%ld %ld %g %g\n",i,j,sy[i*ist.totalhomo+j].re, sy[i*ist.totalhomo+j].im);
      		fprintf (pfz,"%ld %ld %g %g\n",i,j,sz[i*ist.totalhomo+j].re, sz[i*ist.totalhomo+j].im);

		}
	}
//calculate spin matrix elements between all unoccupied (electron) orbitals
	for (a = ist.nlumo; a < ist.nlumo+ist.totallumo; a++) 
	{
		for (b = ist.nlumo; b < ist.nlumo+ist.totallumo; b++) 
		{
			index= sqr(ist.totalhomo)+(a-ist.nlumo)*ist.totallumo+(b-ist.nlumo);
			for (jz = 0; jz < ist.nz; jz++) {
        		for (jy = 0; jy < ist.ny; jy++) {
          			jyz = ist.nx * (ist.ny * jz + jy);
          			for (jx = 0; jx < ist.nx; jx++) {
          				jgridup = jyz + jx;
            			jgriddn = jgridup+ist.ngrid;
						
					//Spin x part
						//<b|r,dn> * <r,up|a>
						sx[index].re += psi[b*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgridup].re
									  + psi[b*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgridup].im;
						sx[index].im += psi[b*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgridup].im
									  - psi[b*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgridup].re;
						//<b|r,up> * <r,dn|a>
						sx[index].re += psi[b*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgriddn].re
									  + psi[b*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgriddn].im;
						sx[index].im += psi[b*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgriddn].im
									  - psi[b*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgriddn].re;
					//Spin y part
						//i*<b|r,dn> * <r,up|a>
						sy[index].im += psi[b*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgridup].re
									  + psi[b*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgridup].im;
						sy[index].re -= psi[b*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgridup].im
									  - psi[b*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgridup].re;
						//-i*<b|r,up> * <r,dn|a>
						sy[index].im -= psi[b*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgriddn].re
									  + psi[b*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgriddn].im;
						sy[index].re += psi[b*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgriddn].im
									  - psi[b*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgriddn].re;
					//Spin z part
						//<b|r,up> * <r,up|a>
						sz[index].re += psi[b*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgridup].re
									  + psi[b*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgridup].im;
						sz[index].im += psi[b*ist.nspinngrid+jgridup].re * psi[a*ist.nspinngrid+jgridup].im
									  - psi[b*ist.nspinngrid+jgridup].im * psi[a*ist.nspinngrid+jgridup].re;
						//<b|r,dn> * <r,dn|a>
						sz[index].re -= psi[b*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgriddn].re
									  + psi[b*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgriddn].im;
						sz[index].im -= psi[b*ist.nspinngrid+jgriddn].re * psi[a*ist.nspinngrid+jgriddn].im
									  - psi[b*ist.nspinngrid+jgriddn].im * psi[a*ist.nspinngrid+jgriddn].re;
					}
				}
			}
			//divide all by 2
			sx[index].re*=0.5*par.dv;
			sx[index].im*=0.5*par.dv;
			sy[index].re*=0.5*par.dv;
			sy[index].im*=0.5*par.dv;
			sz[index].re*=0.5*par.dv;
			sz[index].im*=0.5*par.dv;

			fprintf (pfx,"%ld %ld %g %g\n",a,b,sx[index].re, sx[index].im);
      		fprintf (pfy,"%ld %ld %g %g\n",a,b,sy[index].re, sy[index].im);
      		fprintf (pfz,"%ld %ld %g %g\n",a,b,sz[index].re, sz[index].im);

		}
	}
	fclose(pfx);
	fclose(pfy);
	fclose(pfz);
	return;

}

/****************************************************************************/

void angular(zomplex* lx, zomplex* ly, zomplex* lz, zomplex* lsqr, zomplex* ls, double *vx, double *vy, double *vz, zomplex *psi,
	fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi, long_st ist, par_st par){

	zomplex* Lxpsi = (zomplex*) calloc(ist.nspinngrid,sizeof(zomplex));
	zomplex* Lypsi = (zomplex*) calloc(ist.nspinngrid,sizeof(zomplex));
	zomplex* Lzpsi = (zomplex*) calloc(ist.nspinngrid,sizeof(zomplex));

	zomplex* Lxsqrpsi = (zomplex*) calloc(ist.nspinngrid,sizeof(zomplex));
	zomplex* Lysqrpsi = (zomplex*) calloc(ist.nspinngrid,sizeof(zomplex));
	zomplex* Lzsqrpsi = (zomplex*) calloc(ist.nspinngrid,sizeof(zomplex));
	zomplex* temp1 = (zomplex*) calloc(ist.ngrid,sizeof(zomplex));
	zomplex* temp2 = (zomplex*) calloc(ist.ngrid,sizeof(zomplex));
	long i,j,a,b,jgrid,index, jgridup, jgriddn;
	
	FILE* pfx = fopen("lx.dat", "w"); FILE* pfy = fopen("ly.dat", "w"); FILE* pfz = fopen("lz.dat", "w"); FILE* pfsqr = fopen("lsqr.dat", "w");
	FILE* pfls = fopen("ls.dat", "w");
	printf("Hole States:\n");
	for (i = 0;i<ist.totalhomo;i++){
		
		//spin up part
		lOpp(&Lxpsi[0],&Lypsi[0],&Lzpsi[0], &psi[i*ist.nspinngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		//spin dn part
		lOpp(&Lxpsi[ist.ngrid],&Lypsi[ist.ngrid],&Lzpsi[ist.ngrid], &psi[i*ist.nspinngrid+ist.ngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);


		//Lxsqr parts
		//spin up part
		lOpp(&Lxsqrpsi[0],&temp1[0],&temp2[0], &Lxpsi[0], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		//spin dn part
		lOpp(&Lxsqrpsi[ist.ngrid],&temp1[0],&temp2[0], &Lxpsi[ist.ngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);

		//Lysqr parts
		//spin up part
		lOpp(&temp1[0],&Lysqrpsi[0],&temp2[0], &Lypsi[0], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		//spin dn part
		lOpp(&temp1[0],&Lysqrpsi[ist.ngrid],&temp2[0], &Lypsi[ist.ngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);

		//Lzsqr parts
		//spin up part
		lOpp(&temp1[0],&temp2[0],&Lzsqrpsi[0], &Lzpsi[0], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		//spin dn part
		lOpp(&temp1[0],&temp2[0],&Lzsqrpsi[ist.ngrid], &Lzpsi[ist.ngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);



		

		for (j = 0;j<ist.totalhomo;j++){
			lx[i*ist.totalhomo+j].re=lx[i*ist.totalhomo+j].im=0.0;
			ly[i*ist.totalhomo+j].re=ly[i*ist.totalhomo+j].im=0.0;
			lz[i*ist.totalhomo+j].re=lz[i*ist.totalhomo+j].im=0.0;
			lsqr[i*ist.totalhomo+j].re = lsqr[i*ist.totalhomo+j].im=0.0;


			for(jgrid=0;jgrid<ist.nspinngrid;jgrid++){
				lx[i*ist.totalhomo+j].re -= Lxpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].re + Lxpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].im; // -<j|L_x|i>^*
				lx[i*ist.totalhomo+j].im += Lxpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].re - Lxpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].im;


				ly[i*ist.totalhomo+j].re -= Lypsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].re + Lypsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].im; // -<j|L_y|i>^*
				ly[i*ist.totalhomo+j].im += Lypsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].re - Lypsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].im;

				lz[i*ist.totalhomo+j].re -= Lzpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].re + Lzpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].im; // -<j|L_z|i>^*
				lz[i*ist.totalhomo+j].im += Lzpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].re - Lzpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].im;


				lsqr[i*ist.totalhomo+j].re += Lxsqrpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].re + Lxsqrpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].im; // <j|L_x^2|i>^*
				lsqr[i*ist.totalhomo+j].im -= Lxsqrpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].re - Lxsqrpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].im;

				lsqr[i*ist.totalhomo+j].re += Lysqrpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].re + Lysqrpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].im; // <j|L_y^2|i>^*
				lsqr[i*ist.totalhomo+j].im -= Lysqrpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].re - Lysqrpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].im;

				lsqr[i*ist.totalhomo+j].re += Lzsqrpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].re + Lzsqrpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].im; // <j|L_z^2|i>^*
				lsqr[i*ist.totalhomo+j].im -= Lzsqrpsi[jgrid].im*psi[j*ist.nspinngrid+jgrid].re - Lzsqrpsi[jgrid].re*psi[j*ist.nspinngrid+jgrid].im;



			}
			zomplex lsx,lsy,lsz;
			lsx.re=lsy.re=lsz.re=0.00;
			lsx.im=lsy.im=lsz.im=0.00;

			for(jgrid=0;jgrid<ist.ngrid;jgrid++){
				jgridup = jgrid;
				jgriddn = jgrid+ist.ngrid;


				ls[i*ist.totalhomo+j].re += psi[j*ist.nspinngrid+jgriddn].re * Lxpsi[jgridup].re
										  + psi[j*ist.nspinngrid+jgriddn].im * Lxpsi[jgridup].im;
				ls[i*ist.totalhomo+j].im += psi[j*ist.nspinngrid+jgriddn].re * Lxpsi[jgridup].im
										  - psi[j*ist.nspinngrid+jgriddn].im * Lxpsi[jgridup].re;
				//<j|r,up> * <r,dn|Lxi>
				ls[i*ist.totalhomo+j].re += psi[j*ist.nspinngrid+jgridup].re * Lxpsi[jgriddn].re
										  + psi[j*ist.nspinngrid+jgridup].im * Lxpsi[jgriddn].im;
				ls[i*ist.totalhomo+j].im += psi[j*ist.nspinngrid+jgridup].re * Lxpsi[jgriddn].im
										  - psi[j*ist.nspinngrid+jgridup].im * Lxpsi[jgriddn].re;

				//syly
				//i*<j|r,dn> * <r,up|Lyi>
				ls[i*ist.totalhomo+j].im += psi[j*ist.nspinngrid+jgriddn].re * Lypsi[jgridup].re
							  			  + psi[j*ist.nspinngrid+jgriddn].im * Lypsi[jgridup].im;
				ls[i*ist.totalhomo+j].re -= psi[j*ist.nspinngrid+jgriddn].re * Lypsi[jgridup].im
							  			  - psi[j*ist.nspinngrid+jgriddn].im * Lypsi[jgridup].re;
				//-i*<j|r,up> * <r,dn|Lyi>
				ls[i*ist.totalhomo+j].im -= psi[j*ist.nspinngrid+jgridup].re * Lypsi[jgriddn].re
							  			  + psi[j*ist.nspinngrid+jgridup].im * Lypsi[jgriddn].im;
				ls[i*ist.totalhomo+j].re += psi[j*ist.nspinngrid+jgridup].re * Lypsi[jgriddn].im
							  			  - psi[j*ist.nspinngrid+jgridup].im * Lypsi[jgriddn].re;

  			  	//szlz
  			  	//<j|r,up> * <r,up|Lzi>
				ls[i*ist.totalhomo+j].re += psi[j*ist.nspinngrid+jgridup].re * Lzpsi[jgridup].re
										  + psi[j*ist.nspinngrid+jgridup].im * Lzpsi[jgridup].im;
				ls[i*ist.totalhomo+j].im += psi[j*ist.nspinngrid+jgridup].re * Lzpsi[jgridup].im
										  - psi[j*ist.nspinngrid+jgridup].im * Lzpsi[jgridup].re;
				//<j|r,dn> * <r,dn|Lzi>
				ls[i*ist.totalhomo+j].re -= psi[j*ist.nspinngrid+jgriddn].re * Lzpsi[jgriddn].re
										  + psi[j*ist.nspinngrid+jgriddn].im * Lzpsi[jgriddn].im;
				ls[i*ist.totalhomo+j].im -= psi[j*ist.nspinngrid+jgriddn].re * Lzpsi[jgriddn].im
										  - psi[j*ist.nspinngrid+jgriddn].im * Lzpsi[jgriddn].re;
			  

			}

			//normalize
			lx[i*ist.totalhomo+j].re*=par.dv; ly[i*ist.totalhomo+j].re*=par.dv; lz[i*ist.totalhomo+j].re*=par.dv; lsqr[i*ist.totalhomo+j].re*=par.dv; 
			lx[i*ist.totalhomo+j].im*=par.dv; ly[i*ist.totalhomo+j].im*=par.dv; lz[i*ist.totalhomo+j].im*=par.dv; lsqr[i*ist.totalhomo+j].im*=par.dv;
			
			//normalize and complex conjugate
			ls[i*ist.totalhomo+j].re*=0.5*par.dv; ls[i*ist.totalhomo+j].im*= -0.5*par.dv;

			fprintf(pfx,"%ld\t%ld\t%lf\t%lf\n",i,j,
				lx[i*ist.totalhomo+j].re,lx[i*ist.totalhomo+j].im);
			fprintf(pfy,"%ld\t%ld\t%lf\t%lf\n",i,j,	
				ly[i*ist.totalhomo+j].re,ly[i*ist.totalhomo+j].im);
			fprintf(pfz,"%ld\t%ld\t%lf\t%lf\n",i,j,
				lz[i*ist.totalhomo+j].re,lz[i*ist.totalhomo+j].im);
			fprintf(pfsqr,"%ld\t%ld\t%lf\t%lf\n",i,j,
				lsqr[i*ist.totalhomo+j].re,lsqr[i*ist.totalhomo+j].im);
			fprintf(pfls,"%ld\t%ld\t%lf\t%lf\n",i,j,
				ls[i*ist.totalhomo+j].re,ls[i*ist.totalhomo+j].im);
		}
	}

	printf("Electron States:\n");
	for (a = ist.nlumo;a<ist.totallumo+ist.nlumo;a++){
		//spin up part
		lOpp(&Lxpsi[0],&Lypsi[0],&Lzpsi[0], &psi[a*ist.nspinngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		

		//spin dn part
		lOpp(&Lxpsi[ist.ngrid],&Lypsi[ist.ngrid],&Lzpsi[ist.ngrid], &psi[a*ist.nspinngrid+ist.ngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		

		//Lxsqr parts
		//spin up part
		lOpp(&Lxsqrpsi[0],&temp1[0],&temp2[0], &Lxpsi[0], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		//spin dn part
		lOpp(&Lxsqrpsi[ist.ngrid],&temp1[0],&temp2[0], &Lxpsi[ist.ngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);

		//Lysqr parts
		//spin up part
		lOpp(&temp1[0],&Lysqrpsi[0],&temp2[0], &Lypsi[0], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		//spin dn part
		lOpp(&temp1[0],&Lysqrpsi[ist.ngrid],&temp2[0], &Lypsi[ist.ngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);

		//Lzsqr parts
		//spin up part
		lOpp(&temp1[0],&temp2[0],&Lzsqrpsi[0], &Lzpsi[0], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		//spin dn part
		lOpp(&temp1[0],&temp2[0],&Lzsqrpsi[ist.ngrid], &Lzpsi[ist.ngrid], 
			vx,vy,vz,planfw,planbw, fftwpsi,ist, par);
		

		for (b = ist.nlumo;b<ist.totallumo+ist.nlumo;b++){
			index= sqr(ist.totalhomo)+(a-ist.nlumo)*ist.totallumo+(b-ist.nlumo);

			lx[index].re=lx[index].im=0.0;
			ly[index].re=ly[index].im=0.0;
			lz[index].re=lz[index].im=0.0;
			lsqr[index].re = lsqr[index].im=0.0;


			for(jgrid=0;jgrid<ist.nspinngrid;jgrid++){
				lx[index].re += Lxpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].re + Lxpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].im; // <b|L_x|a>
				lx[index].im += Lxpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].re - Lxpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].im;


				ly[index].re += Lypsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].re + Lypsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].im; // <b|L_y|a>
				ly[index].im += Lypsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].re - Lypsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].im;

				lz[index].re += Lzpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].re + Lzpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].im; // <b|L_z|a>
				lz[index].im += Lzpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].re - Lzpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].im;


				lsqr[index].re += Lxsqrpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].re + Lxsqrpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].im; // <b|L_x^2|a>
				lsqr[index].im += Lxsqrpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].re - Lxsqrpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].im;

				lsqr[index].re += Lysqrpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].re + Lysqrpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].im; // <b|L_y^2|a>
				lsqr[index].im += Lysqrpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].re - Lysqrpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].im;

				lsqr[index].re += Lzsqrpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].re + Lzsqrpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].im; // <b|L_z^2|a>
				lsqr[index].im += Lzsqrpsi[jgrid].im*psi[b*ist.nspinngrid+jgrid].re - Lzsqrpsi[jgrid].re*psi[b*ist.nspinngrid+jgrid].im;


			}

			zomplex lsx,lsy,lsz;
			lsx.re=lsy.re=lsz.re=0.00;
			lsx.im=lsy.im=lsz.im=0.00;

			for(jgrid=0;jgrid<ist.ngrid;jgrid++){
				jgridup = jgrid;
				jgriddn = jgrid+ist.ngrid;

				
				//sxlx
				//<b|r,dn> * <r,up|Lxa>
				ls[index].re += psi[b*ist.nspinngrid+jgriddn].re * Lxpsi[jgridup].re
							  + psi[b*ist.nspinngrid+jgriddn].im * Lxpsi[jgridup].im;
				ls[index].im += psi[b*ist.nspinngrid+jgriddn].re * Lxpsi[jgridup].im
							  - psi[b*ist.nspinngrid+jgriddn].im * Lxpsi[jgridup].re;
				//<b|r,up> * <r,dn|Lxa>
				ls[index].re += psi[b*ist.nspinngrid+jgridup].re * Lxpsi[jgriddn].re
							  + psi[b*ist.nspinngrid+jgridup].im * Lxpsi[jgriddn].im;
				ls[index].im += psi[b*ist.nspinngrid+jgridup].re * Lxpsi[jgriddn].im
							  - psi[b*ist.nspinngrid+jgridup].im * Lxpsi[jgriddn].re;

				//syly
				//i*<j|r,dn> * <r,up|Lyi>
				ls[index].im += psi[b*ist.nspinngrid+jgriddn].re * Lypsi[jgridup].re
				  			  + psi[b*ist.nspinngrid+jgriddn].im * Lypsi[jgridup].im;
				ls[index].re -= psi[b*ist.nspinngrid+jgriddn].re * Lypsi[jgridup].im
				  			  - psi[b*ist.nspinngrid+jgriddn].im * Lypsi[jgridup].re;
				//-i*<j|r,up> * <r,dn|Lyi>
				ls[index].im -= psi[b*ist.nspinngrid+jgridup].re * Lypsi[jgriddn].re
				  			  + psi[b*ist.nspinngrid+jgridup].im * Lypsi[jgriddn].im;
				ls[index].re += psi[b*ist.nspinngrid+jgridup].re * Lypsi[jgriddn].im
				  			  - psi[b*ist.nspinngrid+jgridup].im * Lypsi[jgriddn].re;

  			  	//szlz
  			  	//<j|r,up> * <r,up|Lzi>
				ls[index].re += psi[b*ist.nspinngrid+jgridup].re * Lzpsi[jgridup].re
							  + psi[b*ist.nspinngrid+jgridup].im * Lzpsi[jgridup].im;
				ls[index].im += psi[b*ist.nspinngrid+jgridup].re * Lzpsi[jgridup].im
							  - psi[b*ist.nspinngrid+jgridup].im * Lzpsi[jgridup].re;
				//<j|r,dn> * <r,dn|Lzi>
				ls[index].re -= psi[b*ist.nspinngrid+jgriddn].re * Lzpsi[jgriddn].re
							  + psi[b*ist.nspinngrid+jgriddn].im * Lzpsi[jgriddn].im;
				ls[index].im -= psi[b*ist.nspinngrid+jgriddn].re * Lzpsi[jgriddn].im
							  - psi[b*ist.nspinngrid+jgriddn].im * Lzpsi[jgriddn].re;
							  

			}
			lx[index].re*=par.dv; ly[index].re*=par.dv; lz[index].re*=par.dv; lsqr[index].re*=par.dv; 
			lx[index].im*=par.dv; ly[index].im*=par.dv; lz[index].im*=par.dv; lsqr[index].im*=par.dv;

			//normalize
			ls[index].re*=0.5*par.dv; ls[index].im*=0.5*par.dv;

			fprintf(pfx,"%ld\t%ld\t%lf\t%lf\n",a,b,
				lx[index].re,lx[index].im);
			fprintf(pfy,"%ld\t%ld\t%lf\t%lf\n",a,b,
				ly[index].re,ly[index].im);
			fprintf(pfz,"%ld\t%ld\t%lf\t%lf\n",a,b,
				lz[index].re,lz[index].im);
			fprintf(pfsqr,"%ld\t%ld\t%lf\t%lf\n",a,b,
				lsqr[index].re,lsqr[index].im);
			fprintf(pfls,"%ld\t%ld\t%lf\t%lf\n",a,b,
				ls[index].re,ls[index].im);
		}
	}
	
	free(Lxpsi); free(Lypsi); free(Lzpsi);
	free(Lxsqrpsi); free(Lysqrpsi); free(Lzsqrpsi);
	free(temp1); free(temp2);

}


/************************************************************/
//Calculate the three vector components of the action of the L operator on the 
//spatial part of the grid (no spin part)
void lOpp(zomplex* Lxpsi, zomplex* Lypsi, zomplex* Lzpsi, zomplex* psi, 
	double* vx, double* vy, double* vz,
 	fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,long_st ist, par_st par){
	
	double *kx, *ky, *kz, density, x, y, z;
	long jx,jy,jz, jyz, jxyz;



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
  		z = vz[jz];
  		for (jy = 0;jy<ist.ny;jy++){
  			jyz = ist.nx * (ist.ny * jz + jy);
  			y = vy[jy];
  			for(jx = 0; jx<ist.nx;jx++){
  				x = vx[jx];
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