/*****************************************************************************/
// Required libraries

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <fftw3.h>
#include <omp.h>
#include <mkl.h>

/*****************************************************************************/
// Application specific structures

typedef struct st0 {
  double re, im;
} zomplex;

typedef fftw_plan fftw_plan_loc;

typedef struct st1 {
  double dx, dy, dz, dr, dkx, dky, dkz, dv, dt;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  double Vmin, Vmax, dE, dE_1, VBmin, VBmax, CBmin, CBmax, Ekinmax;
  double Rnlcut2, sigma, sigma_1;
} par_st;

typedef struct st4 {
  long natom, lumo, homo, flaghomo, nnonlocal;
  long ms, ns, nc, npot, mstot, mstotngrid, nspinngrid;
  long nx, ny, nz, ngrid, nspin, nnlc, nproj;
  long natomtype;
  long atomspresent[20];
  long nthreads;
  long flagSO, flagCenter;
  double nx_1, ny_1, nz_1;
} long_st;

typedef struct st5 {
  double x, y, z;
} xyz_st;

typedef struct st9 {
  long natyp;
  int Zval;
  char atyp[3];
  double Vso, strPar, nlcPar[2];
} atm_st;

typedef struct st11 {
  long jxyz;
  zomplex y1[3];
  double proj[5];
  double nlProj[5];
  int nlProjSgn[5];
  double r, r2_1, r2, Vr;
} nlc_st;


/*****************************************************************************/
// Macro definitions

// double macros
#define sqr(x)       ((x) * (x))
#define cube(x)      ((x) * (x) * (x))
#define forth(x)     ((x) * (x) * (x) * (x))
#define ISWAP(a,b)   {double temp = (a); (a) = (b); (b) = -temp;}
#define ISWAP2(a,b)  {double temp = (a); (a) = 2.0 * (b); (b) = -2.0 * temp;}

// complex number macros
#define cplus(a,b,c)  {zomplex tmp; (tmp).re=(a).re+(b).re; (tmp).im=(a).im+(b).im; c=tmp;}
#define cminus(a,b,c) {zomplex tmp; (tmp).re=(a).re-(b).re; (tmp).im=(a).im-(b).im; c=tmp;}
#define cmul(a,b,c)   {zomplex tmp; (tmp).re=(a).re*(b).re-(a).im*(b).im; (tmp).im=(a).re*(b).im+(a).im*(b).re; c=tmp;}
#define cmuls(a,b,c)  {zomplex tmp; (tmp).re=(a).re*(b).re+(a).im*(b).im; (tmp).im=(a).re*(b).im-(a).im*(b).re; c=tmp;}
#define cdev(a,b,c)   {zomplex tmp; double mechane; mechane = (1.0 / ((b).re*(b).re+(b).im*(b).im)); (tmp).re=((a).re*(b).re+(a).im*(b).im) * mechane; (tmp).im = ((a).im * (b).re - (a).re*(b).im)*mechane; c=tmp;}
#define cexp(a,c)     {double texp = exp((a).re); zomplex tmp; (tmp).re=texp * cos((a).im); (tmp).im = texp * sin((a).im); c= tmp;}
#define cexpminx(a,c) {double texp = exp(-(a).re); (c).re=texp * cos((a).im); (c).im = -texp * sin((a).im);}

// constant macros
#define AUTOEV    27.211385
#define PIE       3.14159265358979323846
#define TWOPI     6.28318530717958647692
#define SVDEPS    1.0e-10
#define EPSR      1.0e-10
#define EPSE      1.0e-10
#define EPSR02	  1.0e-2
#define EPSR0     1.0e-4
#define EPSCHI    1.0e-8
#define DENE      1.0e-10
#define EPSDX     1.0e-20
#define ANGTOBOHR 1.889726125

/*****************************************************************************/
// Function declarations

//init.c
void init(double *potl,double *vx,double *vy,double *vz,double *ksqr,double *rx,double *ry,double *rz,atm_st *atm,par_st *par,double *eval,long_st *ist,double *dr,double *vr,double *potatom,long *npot);
void init_psi(zomplex *psi,long_st ist,par_st par,long *idum);
void init_size(long, char *argv[],par_st *,long_st *);
void init_conf(double *rx,double *ry,double *rz,atm_st *atm,par_st *par,long_st *ist);
void init_list(nlc_st *nlc,long *nl,double *vx,double *vy,double *vz,double *rx,double *ry,double *rz,atm_st *atm,par_st par,long_st ist);


//read.c
long assign_atom_number(char atyp[2]);
void assign_atom_type(char *atype,long j);
long get_number_of_atom_types(atm_st *atm,long_st ist,long *list);
void read_conf(double *rx,double *ry,double *rz,atm_st *atm,long n,FILE * pf, long_st *ist);
void read_pot(double *vr,double *pot,long *npot,double *dr,atm_st *atm,long_st *ist);
void setAtmStr(double *rx,double *ry,double *rz,atm_st *atm,long ntot,long_st *ist );
double getBondAngle(long index1,long index2,long index3,double *rx,double *ry,double *rz);

//size.c
double get_dot_ligand_size_z(double *rz, long n);

//interpolate.c
double interpolate(double r,double dr,double *vr,double *pot,long npot,long n,long j);

//rand.c
double ran_nrc(long *idum);
void Randomize();

//hamiltonian.c
void kinetic(zomplex *psi,double *ksqr,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long_st ist);
void spin_orbit_proj_pot(zomplex *phi, zomplex *psi, long_st ist, par_st par, nlc_st *nlc, long* nl );
void nonlocal_proj_pot(zomplex *phi, zomplex *psi, long_st ist, par_st par, nlc_st *nlc, long* nl );
void hamiltonian(zomplex *phi,zomplex *psi,double *potl,double *ksqr,long_st ist,par_st par,nlc_st *nlc,long *nl,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void timeReverseAll(zomplex *psitot, zomplex *dest,long_st ist);
void hamiltonian_t(zomplex *phi, zomplex *psi, double *potl, double *ksqr, long_st ist, par_st par, nlc_st *nlc,
                long *nl, fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,clock_t *timers );

//filter.c
void filtering(zomplex *psi,double *potl,double *ksqr,zomplex *an,double *zn,double *el,long_st ist,par_st par,nlc_st *nlc,long *nl,long tid,long jns);
void filter(zomplex *psin,zomplex *psim1,zomplex *psims,double *potl,double *ksqr,par_st par,nlc_st *nlc,long *nl,zomplex *an,double *zn,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi,long tid,long jns);

//hnorm.c
void hnorm(zomplex *psim1,zomplex *psin,double *potl,double *ksqr,par_st par,nlc_st *nlc,long *nl,double zm1,long_st ist,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);

//norm.c
double norm(zomplex *, double,long,long);
double normalize(zomplex *,double,long,long);
void normalize_all(zomplex *psi,double dr,long ms,long ngrid,long nthreads);

//energy.c
double energy(zomplex *psi,zomplex *phi,double *potl,double *ksqr,long_st ist,par_st par,nlc_st *nlc,long *nl,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void energy_all(zomplex *psi,zomplex *phi,zomplex *psims,double *potl,double *ksqr,double *ene,long_st ist,par_st par,nlc_st *nlc,long *nl,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void get_energy_range(zomplex *psi,zomplex *phi,double *potl,double *vx,double *vy,double *vz,double *ksqr,long_st ist,par_st *par,nlc_st *nlc,long *nl,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
void calc_sigma_E(zomplex *psi,zomplex *phi,zomplex *psitot,double *potl,double *ksqr,double *eval,long_st ist,par_st par,nlc_st *nlc,long *nl,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);


//coeff.c
void coefficient(zomplex *an,double *samp, double* el,par_st par,long_st ist);
void chebyshev_reordered(double *,double,double,long);
double samp_points_ashkenazy(zomplex *point,double min,double max,long nc);
void check_function(zomplex *an,zomplex *samp,long_st ist,par_st par, double el);


//Hmat.c
void Hmat(zomplex *psi,zomplex *phi,MKL_Complex16 *psi0,double *potl,double *ksqr,double *eval,long_st ist,par_st par,nlc_st *nlc,long *nl,fftw_plan_loc planfw,fftw_plan_loc planbw,fftw_complex *fftwpsi);
MKL_Complex16 dotp(zomplex *psi,MKL_Complex16 *phi,long m,long ngrid,double dv);

//nerror.c
void nerror(char *);

//ortho.c
long portho(MKL_Complex16 *,double,long_st);

//projectors.c
void gen_SO_projectors(double dx, double rcut, long nproj, double*  projectors, double* vr);
void gen_nlc_projectors(double dx, double rcut, long nproj, double*  projectors,int* sgnProj, double* vr, atm_st *atm,long jatom,  long_st ist);
//angular.c
void calcAngularExp(zomplex* psitot, double* vx, double* vy, double* vz,
  fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,long_st ist, par_st par, int start, int stop);
void jOpp( zomplex* Jxpsi, zomplex* Jypsi, zomplex* Jzpsi,zomplex* psi, 
  double* vx, double* vy, double* vz,
  fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,long_st ist, par_st par);
void lOpp(zomplex* Lxpsi, zomplex* Lypsi, zomplex* Lzpsi, zomplex* psi, 
  double* vx, double* vy, double* vz,
  fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,long_st ist, par_st par);
void lowPassFilter(zomplex* psi,fftw_plan_loc planfw, fftw_plan_loc planbw, fftw_complex *fftwpsi,
  long_st ist, par_st par);


//write.c
void writeCubeFile(double *rho, par_st par, long_st ist, char *fileName);
/*****************************************************************************/
