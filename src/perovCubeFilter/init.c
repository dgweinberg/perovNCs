#include "fd.h"

/*****************************************************************************/

void init_size(long argc, char *argv[],par_st *par,long_st *ist)
{
  FILE* pf;
  if((pf = fopen("input.par" , "r"))==NULL){
    printf("can't open input.par to read... ");
    nerror("input");
  }
  printf("input as read from input.par\n");
  fscanf (pf,"%ld",&ist->nx);  printf("nx = %ld\t", ist->nx); /*** number of grid point in x ***/
  fscanf (pf,"%ld",&ist->ny);  printf("ny = %ld\t", ist->ny);/*** number of grid point in y ***/
  fscanf (pf,"%ld",&ist->nz);  printf("nz = %ld\n", ist->nz);/*** number of grid point in z ***/
  fscanf (pf,"%lg",&par->dx);  printf("dx = %lg\n", par->dx);/*** dx = minimum grid density ***/
  fscanf (pf,"%ld",&ist->ms);  printf("ms = %ld\t", ist->ms);/*** number of states per filter (energy windows) ***/
  fscanf (pf,"%ld",&ist->ns);  printf("ns = %ld\t", ist->ns);/*** number of filter cycles ***/
  fscanf (pf,"%ld",&ist->nc);    printf("nc = %ld\n", ist->nc);/*** length of newton interpolation ***/
  fscanf (pf,"%lg %lg",&par->VBmin, &par->VBmax); printf("VBmin = %lg\t, VBmzx=%lg\n", par->VBmin, par->VBmax);/***VB energy range***/
  fscanf (pf,"%lg %lg",&par->CBmin, &par->CBmax); printf("CBmin = %lg\t, CBmzx=%lg\n", par->CBmin, par->CBmax);/***CB energy range***/
  fscanf (pf,"%ld",&ist->nthreads);  printf("nthreads = %ld\n", ist->nthreads);
  fscanf (pf,"%ld",&ist->flagSO);  printf("flagSO = %ld\n", ist->flagSO);
  fscanf (pf,"%ld",&ist->flagCenter);  printf("flagCenter = %ld\n", ist->flagCenter);
  fclose(pf);

  // Paramters that rarely, if ever, change - John
  ist->npot = 8192;
  par->Ekinmax = 10.0; 
  ist->nproj = 5;
  // Get the number of atoms 
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ist->natom);
  fclose(pf);
  
  par->dy = par->dz = par->dx;
 
  return;
}

/*****************************************************************************/

void init_conf(double *rx,double *ry,double *rz,atm_st *atm,par_st *par,long_st *ist)
{
  long ntot, ntmp; FILE *pf;
  double xd, yd, zd;
  
  /*** read the pasivated nanocrystal configuration ***/
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ntot);
  assert(fabs((double)(ntot - ist->natom)) < 1.0e-15);
  
  read_conf(rx,ry,rz,atm,ntot,pf, ist);
  fclose(pf);
  printf("nnonlocal = %ld\n",ist->nnonlocal);
  
  xd = rint(0.5 * get_dot_ligand_size_z(rx, ntot) + 5.0);
  yd = rint(0.5 * get_dot_ligand_size_z(ry, ntot) + 5.0);
  zd = rint(0.5 * get_dot_ligand_size_z(rz, ntot) + 5.0);
  printf("xd = %g yd = %g zd = %g\n",xd,yd,zd);

  /***initial parameters for the pot reduce mass, etc. in the x direction ***/
  par->xmin = -xd;
  par->xmax = xd;
  ntmp  = (long)((par->xmax - par->xmin) / par->dx);
  if (ntmp > ist->nx) ist->nx = ntmp;
  par->xmin = -((double)(ist->nx) * par->dx) / 2.0;
  par->xmax = ((double)(ist->nx) * par->dx) / 2.0;
  /*par->dx  = (par->xmax - par->xmin) / (double)(ist->nx);*/
  par->dkx = TWOPI / ((double)ist->nx * par->dx);
  
  /***initial parameters for the pot reduce mass, etc. in the y direction ***/
  par->ymin = -yd;
  par->ymax = yd;
  ntmp  = (long)((par->ymax - par->ymin) / par->dy);
  if (ntmp > ist->ny) ist->ny = ntmp;
  /*par->dy  = (par->ymax - par->ymin) / (double)(ist->ny);*/
  par->ymin = -((double)(ist->ny) * par->dy) / 2.0;
  par->ymax = ((double)(ist->ny) * par->dy) / 2.0;
  par->dky = TWOPI / ((double)ist->ny * par->dy);

  /***initial parameters for the pot reduce mass, etc. in the z direction ***/
  par->zmin = -zd;
  par->zmax = zd;
  ntmp  = (long)((par->zmax - par->zmin) / par->dz);
  if (ntmp > ist->nz) ist->nz = ntmp;
  /*par->dz  = (par->zmax - par->zmin) / (double)(ist->nz);*/
  par->zmin = -((double)(ist->nz) * par->dz) / 2.0;
  par->zmax = ((double)(ist->nz) * par->dz) / 2.0;
  par->dkz = TWOPI / ((double)ist->nz * par->dz);
  printf("new xd = %g yd = %g zd = %g\n",par->xmax,par->ymax,par->zmax);

  ist->nx_1 = 1.0 / (double)(ist->nx);
  ist->ny_1 = 1.0 / (double)(ist->ny);
  ist->nz_1 = 1.0 / (double)(ist->nz);
  ist->ngrid = ist->nx * ist->ny * ist->nz;

  ist->nspin = 2;
  ist->mstot = ist->ms * ist->ns;
  ist->mstotngrid = ist->mstot * ist->ngrid;
  ist->nspinngrid = ist->nspin * ist->ngrid;
  par->sigma = par->dx;
  par->sigma_1 = sqrt(0.5) / par->sigma;
  //par->Rnlcut2 = 0.49 * 6.0 * log(10.0) + 3.0 * par->sigma;
  par->Rnlcut2 = 1.5 +  6.0 * log(10.0) + 3.0 * par->sigma;
  printf("Final Calculation parameters:\n");
  printf("nx = %ld  ny = %ld  nz = %ld npot = %ld\n", ist->nx, ist->ny, ist->nz, ist->npot);
  printf("ms = %ld  ns = %ld natom = %ld nc = %ld\n", ist->ms, ist->ns, ist->natom, ist->nc);
  printf("VBmin = %lg\t, VBmzx=%lg\n", par->VBmin, par->VBmax);
  printf("CBmin = %lg\t, CBmzx=%lg\n", par->CBmin, par->CBmax);
  printf("ngrid = %ld nspin = %ld\n", ist->ngrid, ist->nspin);
  printf("threads = %ld\n", ist->nthreads);
  printf("Rnlcut = %g\n", sqrt(par->Rnlcut2));
  fflush(stdout);

  return;
}

/*****************************************************************************/
void init_list(nlc_st *nlc,long *nl,double *vx,double *vy,double *vz,double *rx,double *ry,double *rz,atm_st *atm,par_st par,long_st ist)
{
  FILE *pf;
  long jatom, jx, jy, jz, jyz, jxyz;
  int iproj, *sgnProj;
  double dx, dy, dz, dxeps, dyeps, dzeps, dr_1, dr2;
  double *SOprojectors, *vr, dr_proj;
  double *nlcprojectors;
  long N = 1024;

  //gen projectors on the fly
  if ((SOprojectors = (double*) calloc(N * ist.nproj, sizeof(double)))==NULL){nerror("mem_projector");}
  if ((vr = (double*) calloc(N, sizeof(double)))==NULL){nerror("mem_vr");}
  gen_SO_projectors(par.dx, sqrt(par.Rnlcut2), ist.nproj, SOprojectors, vr); 
  
  if ((nlcprojectors = (double*) calloc(N * ist.nproj, sizeof(double)))==NULL){nerror("mem_projector");}
  if ((sgnProj = (int*) calloc(ist.nproj, sizeof(int)))==NULL){nerror("mem_projector");}
  
  dr_proj = vr[1];
  printf("Projector dr = %f\n",dr_proj); fflush(0);

  // Useful constants
  double tmp1 = 0.5 * sqrt(3.0 / PIE);
  double tmp2 = 0.5 * sqrt(3.0 / TWOPI);

  //pf = fopen("proj.dat", "w");
  // Find all the grid points within par.Rnlcut of each atom and calculate
  // r, r2, y1[1+m], proj(r), etc. at the grid points and store the results in nlc
  for (jatom = 0; jatom < ist.nnonlocal; jatom++) {
    
    //generate the nonlocal part for each atom
    gen_nlc_projectors(par.dx, sqrt(par.Rnlcut2), ist.nproj, nlcprojectors, sgnProj, vr, atm, jatom, ist);

    nl[jatom] = 0;
    for (jz = 0; jz < ist.nz; jz++) {
      dz = vz[jz] - rz[jatom];
      dzeps = dz + EPSDX;
      for (jy = 0; jy < ist.ny; jy++) {
      	jyz = ist.nx * (ist.ny * jz + jy);
      	dy = vy[jy] - ry[jatom];
      	dyeps = dy + EPSDX;
      	for (jx = 0; jx < ist.nx; jx++) {
      	  jxyz = jyz + jx;
      	  dx = vx[jx] - rx[jatom];
      	  dxeps = dx + EPSDX;
      	  dr2 = dx * dx + dy * dy + dz * dz;
      	  if (dr2 < par.Rnlcut2) {
      	    nlc[jatom*ist.nnlc + nl[jatom]].jxyz = jxyz;

      	    nlc[jatom*ist.nnlc + nl[jatom]].r  = sqrt(dr2);

      	    dr_1 = 1.0 / sqrt(dx * dx + dy * dy + dzeps * dzeps);
      	    nlc[jatom*ist.nnlc + nl[jatom]].y1[1].re = tmp1 * dzeps * dr_1;
      	    nlc[jatom*ist.nnlc + nl[jatom]].y1[1].im = 0.0;

      	    dr_1 = 1.0 / sqrt(dxeps * dxeps + dy * dy + dz * dz);
      	    nlc[jatom*ist.nnlc + nl[jatom]].y1[2].re  = -tmp2 * dxeps * dr_1;
      	    nlc[jatom*ist.nnlc + nl[jatom]].y1[0].re = tmp2 * dxeps * dr_1;

      	    dr_1 = 1.0 / sqrt(dx * dx + dyeps * dyeps + dz * dz);
      	    nlc[jatom*ist.nnlc + nl[jatom]].y1[2].im  = -tmp2 * dyeps * dr_1;
      	    nlc[jatom*ist.nnlc + nl[jatom]].y1[0].im = -tmp2 * dyeps * dr_1;


            //write projectors to nlc struct and scale projectors by the SO scaling for this atom
            for (iproj = 0; iproj< ist.nproj; iproj++){ 
                nlc[jatom*ist.nnlc + nl[jatom]].proj[iproj] = 
                  interpolate(sqrt(dr2),dr_proj,vr, &SOprojectors[N*iproj],0, N,0);
                
                nlc[jatom*ist.nnlc + nl[jatom]].proj[iproj] *= sqrt(atm[jatom].Vso);
                
                nlc[jatom*ist.nnlc + nl[jatom]].nlProj[iproj] =
                  interpolate(sqrt(dr2),dr_proj,vr, &nlcprojectors[N*iproj],0, N,0);
                nlc[jatom*ist.nnlc + nl[jatom]].nlProjSgn[iproj] = sgnProj[iproj];
            }

      	    if (dr2 > EPSDX) {
      	      nlc[jatom*ist.nnlc + nl[jatom]].r2_1 = sqr(dr_1);
      	      nlc[jatom*ist.nnlc + nl[jatom]].r2 = dr2;
      	    }
      	    else {
      	      nlc[jatom*ist.nnlc + nl[jatom]].r2_1 = 0.0;
      	      nlc[jatom*ist.nnlc + nl[jatom]].r2 = 1.0 / EPSDX;
      	    }
      	    nl[jatom]++;
      	  }
      	}
      }
    }
  }
  //fclose(pf);

  free(SOprojectors);
  free(vr);
  free(nlcprojectors);
  free(sgnProj);

  pf = fopen("list.dat" , "w");
  for (jatom = 0; jatom < ist.nnonlocal; jatom++) {
    fprintf(pf, "%ld %ld\n", jatom, nl[jatom]);
  }
  fclose(pf);
  
  return;
}


/*****************************************************************************/

void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,atm_st *atm,par_st *par,double *eval,long_st *ist,double *dr,double *vr,double *potatom,long *npot)
{
  FILE *pf;
  long jx, jy, jz, jyz, jxyz, jatom, jtmp, nn, msVB, msCB;
  double range,del, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
  double sum;

  if ((ksqrx  = (double*)calloc(ist->nx,sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry  = (double*)calloc(ist->ny,sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz  = (double*)calloc(ist->nz,sizeof(double)))==NULL)nerror("ksqrz");
  
  par->dv = par->dx * par->dy * par->dz;
  par->dr = sqrt(sqr(par->dx) + sqr(par->dy) + sqr(par->dz));
  printf("dx = %g dy = %g dz = %g dv = %g dr = %g\n", par->dx, par->dy, par->dz, par->dv, par->dr);

  /***initializing the ksqr vectors ***/
  //hold extra factor of 0.5 for the kinetic energy operation
  for (ksqrx[0] = 0.0, jx = 1; jx <= ist->nx / 2; jx++)
    ksqrx[jx] = (ksqrx[ist->nx-jx] = 0.5 * sqr((double)(jx) * par->dkx) *
		ist->nx_1 * ist->ny_1 * ist->nz_1);

  for (ksqry[0] = 0.0, jy = 1; jy <= ist->ny / 2; jy++)
    ksqry[jy] = (ksqry[ist->ny-jy] = 0.5 * sqr((double)(jy) * par->dky) *
		ist->ny_1 * ist->nx_1 * ist->nz_1);

  for (ksqrz[0] = 0.0, jz = 1; jz <= ist->nz / 2; jz++)
    ksqrz[jz] = (ksqrz[ist->nz-jz] = 0.5 * sqr((double)(jz) * par->dkz) *
		ist->nz_1 * ist->nx_1 * ist->ny_1);

  par->Ekinmax *= (ist->ny_1 * ist->nx_1 * ist->nz_1);
  for (jz = 0; jz < ist->nz; jz++) for (jy = 0; jy < ist->ny; jy++){
    for (jyz = ist->nx * (ist->ny * jz + jy), jx = 0; jx < ist->nx; jx++){
      jxyz = jyz + jx;
      ksqr[jxyz] = ksqrx[jx] + ksqry[jy] + ksqrz[jz];
      if (ksqr[jxyz] > par->Ekinmax) ksqr[jxyz] = par->Ekinmax;
    }
  }
  free(ksqrx); free(ksqry);  free(ksqrz);
  


  printf("read pot.."); fflush(0);
  /*** read pseudopotentials ***/


  read_pot(vr, potatom, npot, dr, atm, ist);
  printf("done\n"); fflush(0);
  


  /***initializing the potential vector  ***/
  for (jx = 0, dx = par->xmin; jx < ist->nx; jx++, dx += par->dx) vx[jx] = dx;
  for (jy = 0, dy = par->ymin; jy < ist->ny; jy++, dy += par->dy) vy[jy] = dy;
  for (jz = 0, dz = par->zmin; jz < ist->nz; jz++, dz += par->dz) vz[jz] = dz;

  omp_set_dynamic(0);
  omp_set_num_threads(ist->nthreads);
#pragma omp parallel for private(dx,dy,dz,del,jy,jx,jyz,jxyz,sum,jatom)
  for (jz = 0; jz < ist->nz; jz++) {
    for (jy = 0; jy < ist->ny; jy++) {
      jyz = ist->nx * (ist->ny * jz + jy);
      for (jx = 0; jx < ist->nx; jx++) {
      	jxyz = jyz + jx;
      	for (sum = 0.0, jatom = 0; jatom < ist->natom; jatom++){
      	  dx = vx[jx] - rx[jatom];
      	  dy = vy[jy] - ry[jatom];
      	  dz = vz[jz] - rz[jatom];
      	  del = sqrt(dx * dx + dy * dy + dz * dz);

          //cubic part of the function
          sum += (1.0-atm[jatom].strPar)*interpolate(del,dr[2*atm[jatom].natyp],vr,potatom,ist->npot,npot[2*atm[jatom].natyp],2*atm[jatom].natyp);

          //ortho part of the function
          sum += (atm[jatom].strPar)*interpolate(del,dr[2*atm[jatom].natyp+1],vr,potatom,ist->npot,npot[2*atm[jatom].natyp+1],2*atm[jatom].natyp+1);

      	}
      	potl[jxyz] = sum;
      }
    }
  }

  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;
  for (jxyz= 0; jxyz < ist->ngrid; jxyz++){
    if (par->Vmax < potl[jxyz]) par->Vmax = potl[jxyz];
    if (par->Vmin > potl[jxyz]) par->Vmin = potl[jxyz];
  }

  printf("dV = %g Vmin = %g Vmax = %g\n", par->Vmax-par->Vmin, par->Vmin, par->Vmax);

  par->dE = 0.5 * sqr(PIE) / (par->dx*par->dx) + 0.5 * sqr(PIE) / (par->dy*par->dy) + 0.5 * sqr(PIE) / (par->dz*par->dz);
  printf("dT = %g\n", par->dE);
  
  /*** setting the energy grid El ***/
  range = (par->CBmax - par->CBmin)+(par->VBmax - par->VBmin);
  msCB = (long) ((double)ist->ms / 2.0);
  msVB = ist->ms-msCB;
  if(msCB<1 || msVB<1){printf("error: init egrid\n"); fflush(0);exit(EXIT_FAILURE);}
  del = (par->VBmax - par->VBmin)/(double)msVB;
  printf("Del VB: %lg\n", del);
  for (jx = 0; jx < msVB; jx++) {
    eval[jx] = par->VBmax - (double)(jx) * del;
  }
  del = (par->CBmax - par->CBmin)/(double)msCB;
  printf("Del CB: %lg\n", del);
  for (jx = msVB; jx < ist->ms; jx++) {
    eval[jx] = par->CBmin + (double)(jx-msVB) * del;
  }

  for (pf = fopen("Egrid.dat","w"),jx = 0; jx < ist->ms; jx++) {
    fprintf(pf, "%g\n", eval[jx]);
  }
  fclose(pf);
  
  for (nn = jatom = 0; jatom < ist->natom; jatom++)
    if (((atm[jatom].natyp < 8) || (atm[jatom].natyp > 11)) && (atm[jatom].natyp % 2)) nn++;
  printf("nn = %ld\n", nn);
  ist->homo = 4*nn-1;
  ist->lumo = ist->homo+1;
  printf("homo = %ld lumo = %ld\n", ist->homo, ist->lumo);

  
  /*** for the nonlocal potential ***/
  /*** count the max number of grid points in Rnlcut of an atom***/
  for (ist->nnlc = 0, jatom = 0; jatom < ist->nnonlocal; jatom++) {
    for (jtmp =0, jz = 0; jz < ist->nz; jz++) {
      for (jy = 0; jy < ist->ny; jy++) {
      	for (jx = 0; jx < ist->nx; jx++) {
      	  dx = vx[jx] - rx[jatom];
      	  dy = vy[jy] - ry[jatom];
      	  dz = vz[jz] - rz[jatom];
      	  if (dx*dx+dy*dy+dz*dz < par->Rnlcut2) {
            jtmp++; 
          }
      	}
      }
    }
    if (jtmp > ist->nnlc) ist->nnlc = jtmp;
  }
  printf("nnlc = %ld\n", ist->nnlc);
  fflush(stdout);

  return;
}
/*****************************************************************************/

void init_psi(zomplex *psi,long_st ist,par_st par,long *idum)
{
  long jx, jy, jz, jzy, jxyz;
  long tidum = (*idum);

  for (jz = 0; jz < ist.nz; jz++) for (jy = 0; jy < ist.ny; jy++) {
    for (jzy = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++) {
      jxyz = jzy + jx;
      psi[jxyz].re = (-1.0 + 2.0 * ran_nrc(&tidum));
      psi[jxyz].im = (-1.0 + 2.0 * ran_nrc(&tidum));
    }
  }

  normalize(psi,par.dv,ist.ngrid,ist.nthreads);
  (*idum) = tidum;

  return;
}

/*****************************************************************************/
