/*****************************************************************************/

#include "fd.h"
#include "unistd.h"

/*****************************************************************************/

void init_size(long argc, char *argv[],par_st *par,long_st *ist) {
  long ieof, i, a, nval; 
  double evalloc, deloc, *eval, *de;
  char field[100], tmp[100];
  FILE *pf = fopen("input.par" , "r");
  if (pf ==NULL) {
    printf("no input in cwd\n");
    exit(EXIT_FAILURE);
  }
  fscanf (pf,"%ld %lf",&ist->nx, &par->dx);  /*** number of grid point in x and dx ***/
  printf("ist->nx = %ld\tist->dx = %lf\n", ist->nx, par->dx);fflush(0);
  fscanf (pf,"%ld %lf",&ist->ny, &par->dy);  /*** number of grid polong in y ***/
  printf("ist->ny = %ld\tist->dy = %lf\n", ist->ny, par->dy); fflush(0);
  fscanf (pf,"%ld %lf",&ist->nz, &par->dz);  /*** number of grid polong in z ***/
  printf("ist->nz = %ld\tist->dz = %lf\n", ist->nz, par->dz);fflush(0);
  fscanf (pf,"%lg %lg %lg",&par->epsx,&par->epsy,&par->epsz);  /*** dielectric constant for the x,y,z directions ***/
  printf("par->espz = %lg par->espz = %lg par->espz = %lg\n", par->epsx, par->epsy, par->epsz);fflush(0);
  fscanf (pf,"%lg %lg %lg",&par->deltae,&par->deltah,&par->deps);  /*** energy windows for h+ and e- states allowed and sige check ***/
  fscanf (pf,"%ld %ld", &ist->maxElecStates, &ist->maxHoleStates);  /*** maximum number of single-particle states used ***/
  fscanf (pf,"%ld",&ist->nthreads);  /*** number of grid polong in z ***/
  fclose(pf);

  // Set Defaults
  par->Ekinmax = 10.0; // kinetic energy maximum 
  ist->npot = 8192;    // length of pseudopotential files
  ist->natomtype = 12;
  ist->printFPDensity = 0; // do not print fixed point quasiparticle densities
  ist->calcDarkStates = 0; // calculate the bright excitonic states (with exchange-like term)

  // Optional functionality of the code, optional inputs
  if( access( "optionalInput.par", F_OK) != -1 ) {
    pf = fopen("optionalInput.par", "r");
    while (fscanf(pf, "%s", field) != EOF) {
      if (! strcmp(field, "printFPDensity")) {
  	    fscanf(pf, "%s %ld", tmp, &(ist->printFPDensity));
  	  }
      else if (! strcmp(field, "calcDarkStates")) {
        fscanf(pf, "%s %ld", tmp, &(ist->calcDarkStates));
      }
      else {
        printf("Invalid input field and/ or format - equal sign required after each field\n");
        printf("Only allowed fields are (case-sensitive): printFPDensity, calcDarkStates\n");
        printf("printFPDensity = 1 (= 0 is default)\n");
        printf("calcDarkStates = 1 (= 0 is default)\n");
        exit(EXIT_FAILURE);
      }
    }
    fclose(pf);
  }

  
  // Error checking related to the input.par file
  if (par->epsx < 0.0) nerror("error in epsilonx");
  if (par->epsy < 0.0) nerror("error in epsilony");
  if (par->epsz < 0.0) nerror("error in epsilonz");


  // Get the total number of atoms from the conf.par
  pf = fopen("conf.par" , "r");
  if (pf ==NULL) {
    printf("no conf in cwd\n");
    exit(EXIT_FAILURE);
  }
  fscanf(pf, "%ld", &ist->natom);
  fclose(pf);

  // Calculate the parameters that depend on the parameters from input.par 
  ist->nx_1 = 1.0 / (double)(ist->nx);
  ist->ny_1 = 1.0 / (double)(ist->ny);
  ist->nz_1 = 1.0 / (double)(ist->nz);
  ist->ngrid = ist->nx * ist->ny * ist->nz;
  ist->nspinngrid = 2*ist->nx * ist->ny * ist->nz;
  ist->ngrid_1 = 1.0 / (double)(ist->ngrid);

  //
  ist->nhomo = ist->nlumo = ist->totalhomo = ist->totallumo = 0;
  pf = fopen("eval.par" , "r");
  if (pf ==NULL) {
    printf("no eval in cwd\n");
    exit(EXIT_FAILURE);
  }
  
  //this is where we set which is eleectron and what is hole
  //TODO: Have the fermi energy setable!!

  for (i = ieof = 0; ieof != EOF; i++){
    ieof = fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
    if (deloc < par->deps && evalloc < -0.19) ist->nhomo = i;
    if (deloc < par->deps && evalloc < -0.19) ist->totalhomo++; 
    if (deloc < par->deps && evalloc > -0.19) ist->totallumo++; 
  }
  fclose(pf);

  nval = i - 1;
  pf = fopen("eval.par" , "r");
  for (i = 0; i <= ist->nhomo; i++) fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
  for (i = ist->nhomo+1; i < nval; i++) {
    fscanf(pf, "%ld %lg %lg", &a, &evalloc, &deloc);
    if (deloc < par->deps) {
      ist->nlumo = i;
      break;
    }
  }
  fclose(pf);

  // 
  printf("The number of openMP threads used = %ld\n", ist->nthreads);
  printf("Total # of filtered hole eigenstates = %ld\n", ist->totalhomo);
  printf("Total # of filtered electron eigenstates = %ld\n", ist->totallumo);
  printf("The eval.par index of the HOMO state = %ld\n", ist->nhomo);
  printf("The eval.par index of the LUMO state = %ld\n", ist->nlumo);

  // 
  eval = (double *) calloc(nval, sizeof(double)); 
  de = (double *) calloc(nval, sizeof(double)); 

  pf = fopen("eval.par" , "r");
  for (i = 0; i < nval; i++) fscanf(pf, "%ld %lg %lg", &a, &eval[i], &de[i]);
  fclose(pf);

  printf("The HOMO energy = % .8f % .8f\n", eval[ist->nhomo], eval[ist->nhomo]*AUTOEV);
  printf("The LUMO energy = % .8f % .8f\n", eval[ist->nlumo], eval[ist->nlumo]*AUTOEV);
  printf("Fundamental gap = %.10f %.10f\n", eval[ist->nlumo]-eval[ist->nhomo], 
    (eval[ist->nlumo]-eval[ist->nhomo])*AUTOEV);

  // 
  for (ist->totalhomo = 0, i = ist->nhomo; i >= 0; i--)
    if ((eval[ist->nhomo] - eval[i] <= par->deltah) && de[i] < par->deps) ist->totalhomo++;

  for (ist->totallumo = 0, i = ist->nlumo; i < nval; i++)
    if ((eval[i] - eval[ist->nlumo] <= par->deltae) && de[i] < par->deps) ist->totallumo++; 
  ist->ms = ist->totalhomo + ist->totallumo;

  // 
  printf("The total number of atoms = %ld\n", ist->natom);
  printf("Number of grid points used: nx = %ld  ny = %ld  nz = %ld\n", ist->nx, ist->ny, ist->nz);
  printf("Length of pseudopotential files, npot = %ld\n", ist->npot);
  printf("Kinetic energy maximum = %.4f\n", par->Ekinmax);
  printf("Dielectric Constant: epsx = %.4f epsy = %.4f epsz = %.4f\n", par->epsx, par->epsy, par->epsz);
  printf("Maximum sigma for an eigenstate, deps = %.4f\n", par->deps);
  printf("Energy windows, electrons = %.4f holes = %.4f\n", par->deltae, par->deltah);
  printf("The number of hole eigenstates within %.4f of the HOMO energy = %ld\n", par->deltah, ist->totalhomo);
  printf("The number of electron eigenstates within %.4f of the LUMO energy = %ld\n", par->deltae, ist->totallumo);
  printf("Total number of carrier states, ms = %ld\n", ist->ms);
  
  // Free the memory that was dynamically allocated within this function
  free(eval);  free(de);

  return;
}

/*****************************************************************************/

void init(double *potl, double *vx, double *vy, double *vz, double *ksqr, double *rx, double *ry, double *rz, par_st *par, long_st *ist)
{
  FILE *pf; 
  long ntmp, jx, jy, jz, jyz, jxyz, ie, ntot, jp, *npot, nn, flags=0;
  double del, mx, my, mz, xd, yd, zd, dx, dy, dz, *ksqrx, *ksqry, *ksqrz;
  double *vr, *potatom, *dr;
  atm_st *atm;

  // Allocate memory 
  if ((ksqrx = (double *) calloc(ist->nx, sizeof(double)))==NULL)nerror("ksqrx");
  if ((ksqry = (double *) calloc(ist->ny, sizeof(double)))==NULL)nerror("ksqry");
  if ((ksqrz = (double *) calloc(ist->nz, sizeof(double)))==NULL)nerror("ksqrz");

  // Read the passivated nanocrystal configuration from conf.par
  pf = fopen("conf.par" , "r");
  fscanf(pf,"%ld",&ntot);
  assert(fabs((double)(ntot - ist->natom)) < 1.0e-15);
  if ((atm = (atm_st *) calloc(ist->natom, sizeof(atm_st)))==NULL)nerror("atm");
  read_conf(rx,ry,rz,atm,ntot,pf);
  fclose (pf);

  // Set the box size  
  xd = rlong(0.5 * get_dot_ligand_size_z(rx,ntot) + 5.0);
  yd = rlong(0.5 * get_dot_ligand_size_z(ry,ntot) + 5.0);
  zd = rlong(0.5 * get_dot_ligand_size_z(rz,ntot) + 5.0);
  printf("QD (minimum) dimensions: xd = %.2f yd = %.2f zd = %.2f\n", xd, yd, zd);
  
  /// Initial parameters for the pot reduce mass, etc. in the x direction 
  par->xmin = -xd;
  par->xmax = xd;
  mx = 1.0;
  
  ntmp  = (long)((par->xmax - par->xmin) / par->dx);
  if (ntmp > ist->nx) ist->nx = ntmp;
  par->xmin = -((double)(ist->nx) * par->dx) / 2.0;
  par->xmax = ((double)(ist->nx) * par->dx) / 2.0;
  par->dkx = TWOPI / ((double)ist->nx * par->dx);
  
  // Initial parameters for the pot reduce mass, etc. in the y direction 
  par->ymin = -yd;
  par->ymax = yd;
  my = 1.0;
  
  ntmp  = (long)((par->ymax - par->ymin) / par->dy);
  if (ntmp > ist->ny) ist->ny = ntmp;
  par->ymin = -((double)(ist->ny) * par->dy) / 2.0;
  par->ymax = ((double)(ist->ny) * par->dy) / 2.0;
  par->dky = TWOPI / ((double)ist->ny * par->dy);

  // Initial parameters for the pot reduce mass, etc. in the z direction 
  par->zmin = -zd;
  par->zmax = zd;
  mz = 1.0;
  ntmp  = (long)((par->zmax - par->zmin) / par->dz);
  if (ntmp > ist->nz) ist->nz = ntmp;
  /*par->dz  = (par->zmax - par->zmin) / (double)(ist->nz);*/
  par->zmin = -((double)(ist->nz) * par->dz) / 2.0;
  par->zmax = ((double)(ist->nz) * par->dz) / 2.0;
  par->dkz = TWOPI / ((double)ist->nz * par->dz);
  
  printf("Final Box Dimensions: xd = %g yd = %g zd = %g\n",par->xmax,par->ymax,par->zmax);

  if ((xd < yd) && (xd < zd))  {
    par->minr = xd;
    par->boxl = (double)(ist->nx) * par->dx;
  }
  else if ((yd < xd) && (yd < zd))  {
    par->minr = yd;
    par->boxl = (double)(ist->ny) * par->dy;
  }
  else {
    par->minr = zd;
    par->boxl = (double)(ist->nz) * par->dz;
  }    
  
  par->gamma = 7.0 / (2.0 * par->minr);
  par->gamma2 = sqr(par->gamma);

  printf("Gamma = %.6f\n", par->gamma);
  
  par->Vmin = 1.0e10;
  par->Vmax = -1.0e10;

  par->dv = par->dx * par->dy * par->dz;
  par->dr = sqrt(sqr(par->dx) + sqr(par->dy) + sqr(par->dz));
  printf("Grid point spacing: dx = %.4f dy = %.4f dz = %.4f dv = %.6f dr = %.6f\n",
	  par->dx, par->dy, par->dz, par->dv, par->dr);


  /***initializing the ksqr vectors ***/
  for (ksqrx[0] = 0.0, jx = 1; jx <= ist->nx / 2; jx++)
    ksqrx[jx] = (ksqrx[ist->nx-jx] = 0.5 * sqr((double)(jx) * par->dkx) *
		ist->nx_1 * ist->ny_1 * ist->nz_1 / mx);

  for (ksqry[0] = 0.0, jy = 1; jy <= ist->ny / 2; jy++)
    ksqry[jy] = (ksqry[ist->ny-jy] = 0.5 * sqr((double)(jy) * par->dky) *
		ist->ny_1 * ist->nx_1 * ist->nz_1 / my);

  for (ksqrz[0] = 0.0, jz = 1; jz <= ist->nz / 2; jz++)
    ksqrz[jz] = (ksqrz[ist->nz-jz] = 0.5 * sqr((double)(jz) * par->dkz) *
		ist->nz_1 * ist->nx_1 * ist->ny_1 / mz);

  par->Ekinmax *= (ist->ny_1 * ist->nx_1 * ist->nz_1);
  for (jz = 0; jz < ist->nz; jz++) for (jy = 0; jy < ist->ny; jy++){
    for (jyz = ist->nx * (ist->ny * jz + jy), jx = 0; jx < ist->nx; jx++){
      jxyz = jyz + jx;
      ksqr[jxyz] = ksqrx[jx] + ksqry[jy] + ksqrz[jz];
      if (ksqr[jxyz] > par->Ekinmax) ksqr[jxyz] = par->Ekinmax;
    }
  }
  free(ksqrx); free(ksqry);  free(ksqrz); free(atm);

  /***initializing the potential vector  ***/
  for (jx = 0, dx = par->xmin; jx < ist->nx; jx++, dx += par->dx) vx[jx] = dx;
  for (jy = 0, dy = par->ymin; jy < ist->ny; jy++, dy += par->dy) vy[jy] = dy;
  for (jz = 0, dz = par->zmin; jz < ist->nz; jz++, dz += par->dz) vz[jz] = dz;

  
  return;
}


/****************************************************************************/

void init_pot(double *vx,double *vy,double *vz,zomplex *potq,zomplex *potqx,par_st par,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi)
{
  long jx, jy, jz, jyz, jxyz, sx, sy, sz;
  double dr, dr2, x2, y2, z2, *kx2, *ky2, *kz2, alpha, cosa, sina;
  double ex_1, ey_1, ez_1, sqrtexeyez_1, sqrtaveps, boxl2 = sqr(par.boxl), sqrk0;
  double gammaeps;
  zomplex *potr, *potrx, tmp;

  ex_1 = 1.0 / par.epsx;
  ey_1 = 1.0 / par.epsy;
  ez_1 = 1.0 / par.epsz;
  sqrtexeyez_1 = 1.0 / sqrt(par.epsx * par.epsy * par.epsz);
  sqrtaveps = sqrt((par.epsx + par.epsy + par.epsz) / 3.0);
  gammaeps = par.gamma * sqrtaveps;

  /*** no yukawa screening for the exchange ***/
  sqrk0 = par.gamma2 * sqr(sqrtaveps);

  if ((kx2  = (double*)calloc(ist.nx,sizeof(double)))==NULL)nerror("kx2");
  if ((ky2  = (double*)calloc(ist.ny,sizeof(double)))==NULL)nerror("ky2");
  if ((kz2  = (double*)calloc(ist.nz,sizeof(double)))==NULL)nerror("kz2");
  if ((potr = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("potr");
  if ((potrx = (zomplex*)calloc(ist.ngrid,sizeof(zomplex)))==NULL)nerror("potrx");
  
  for (jxyz = 0; jxyz < ist.ngrid; jxyz++) potr[jxyz].re = potr[jxyz].im = 0.0;
  for (jxyz = 0; jxyz < ist.ngrid; jxyz++) potrx[jxyz].re = potrx[jxyz].im = 0.0;

  for (jz = 0; jz < ist.nz; jz++) {
    z2 =sqr(vz[jz]);
    for (jy = 0; jy < ist.ny; jy++) {
      y2 = sqr(vy[jy]);
      jyz = ist.nx * (ist.ny * jz + jy);
      for (jx = 0; jx < ist.nx; jx++) {
      	x2 = sqr(vx[jx]);
      	jxyz = jyz + jx;
      	dr2 = (x2 + y2 + z2);
      	if (dr2 < boxl2) {
      	  dr = sqrt(x2 + y2 + z2);
      	  potr[jxyz].re = screenedcoulomb(dr, par.gamma);
      	  dr = sqrt(ex_1 * x2 + ey_1 * y2 + ez_1 * z2);
      	  potrx[jxyz].re = sqrtexeyez_1 * screenedcoulomb(dr, gammaeps);
      	}
      }
    }
  }

  for (jxyz = 0; jxyz < ist.ngrid; jxyz++)
    potq[jxyz].re = potq[jxyz].im = potqx[jxyz].re = potqx[jxyz].im = 0.0;
  
  memcpy(&fftwpsi[0], &potr[0], ist.ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&potq[0], &fftwpsi[0], ist.ngrid*sizeof(potq[0]));

  memcpy(&fftwpsi[0], &potrx[0], ist.ngrid*sizeof(fftwpsi[0]));
  fftw_execute(planfw);
  memcpy(&potqx[0], &fftwpsi[0], ist.ngrid*sizeof(potqx[0]));
  
  for (kx2[0] = 0.0, jx = 1; jx <= ist.nx / 2; jx++)
    kx2[jx] = (kx2[ist.nx-jx] = sqr((double)(jx) * par.dkx));
  for (ky2[0] = 0.0, jy = 1; jy <= ist.ny / 2; jy++)
    ky2[jy] = (ky2[ist.ny-jy] = sqr((double)(jy) * par.dky));
  for (kz2[0] = 0.0, jz = 1; jz <= ist.nz / 2; jz++)
    kz2[jz] = (kz2[ist.nz-jz] = sqr((double)(jz) * par.dkz));

  for (sz = 1.0, jz = 0; jz < ist.nz; jz++, sz = -sz) {
    z2 = kz2[jz];
    for (sy = 1.0, jy = 0; jy < ist.ny; jy++, sy = -sy) {
      y2 = ky2[jy];
      jyz = ist.nx * (ist.ny * jz + jy);
      for (sx = 1.0, jx = 0; jx < ist.nx; jx++, sx = -sx) {
      	x2 = kx2[jx];
      	jxyz = jyz + jx;
      	alpha = PIE * (double)(jx + jy + jz + ist.ngrid / 2);
      	cosa = cos(alpha);
      	sina = sin(alpha);

      	/*** hartree term ***/
      	tmp.re = potq[jxyz].re;
      	tmp.im = potq[jxyz].im;
      	
      	potq[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par.dv;
      	potq[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par.dv;
      	potq[jxyz].re += (FOURPI / (x2 + y2 + z2 + par.gamma2));
      	potq[jxyz].re *= ist.ngrid_1;
      	potq[jxyz].im *= ist.ngrid_1;

      	/*** screened exchange term ***/
      	tmp.re = potqx[jxyz].re;
      	tmp.im = potqx[jxyz].im;

      	potqx[jxyz].re = (tmp.re * cosa - tmp.im * sina) * par.dv;
      	potqx[jxyz].im = (tmp.re * sina + tmp.im * cosa) * par.dv;
      	//potqx[jxyz].re += FOURPI / (par.epsx * x2 + par.epsy * y2 + par.epsz * z2 + sqrk0);
      	potqx[jxyz].re += (FOURPI * (1.0 - exp(-0.25* (par.epsx * x2 + par.epsy * y2 + par.epsz * z2) / sqrk0)) / (par.epsx * x2 + par.epsy * y2 + par.epsz * z2 + EPSR));
        //printf("denominator = % .12f\n", par.epsx * x2 + par.epsy * y2 + par.epsz * z2);
		    potqx[jxyz].re *= ist.ngrid_1;
      	potqx[jxyz].im *= ist.ngrid_1;
      }
    }
  }

  free(potr);  free(potrx); free(kx2); free(ky2); free(kz2);

  return;
}

/************************************************************************/

#define SQRTPI       (sqrt(3.14159265358979323846))

double screenedcoulomb(double dr, double gamma)
{
  //if (dr < EPSR) return (gamma);
  if (dr < EPSR) return (2.0*gamma/SQRTPI);
  return (erf(gamma * dr) / dr);
  //return ((1.0 - exp(-gamma * dr)) / dr);
}

/************************************************************************/

void init_psi(zomplex *psi,double *vx,double *vy,double *vz,long_st ist,par_st par,long *idum)
{
  long jx, jy, jz, jzy, jxyz;
  long tidum = (*idum);

  for (jz = 0; jz < ist.nz; jz++) for (jy = 0; jy < ist.ny; jy++){
    for (jzy = ist.nx * (ist.ny * jz + jy), jx = 0; jx < ist.nx; jx++){
      jxyz = jzy + jx;
      psi[jxyz].re = (-1.0 + 2.0 * ran_nrc(&tidum));
      psi[jxyz].im = 0.0;
    }
  }
  normalize_zomplex(psi, par.dv, ist.ngrid);
  (*idum) = tidum;
  return;
}

/************************************************************************/
