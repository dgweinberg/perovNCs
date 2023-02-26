/*****************************************************************************/
// Main file for cube printing utility.
#include "fd.h"

/*****************************************************************************/
int countlines(char *filename);


/*****************************************************************************/
int main(int argc, char *argv[])
{
  FILE *ppsi; zomplex *psitot;
  par_st par;   long_st  ist; atm_st *atm;
  double *rx, *ry, *rz;
  double *rho; 
  long  jms;
  int start, end;
  time_t currentTime = time(NULL);

  //command line input parsing
  if (argc!=3){
    printf("Usage: makecube start end");
    exit(EXIT_FAILURE);
  }
  start = atoi(argv[1]);
  end = atoi(argv[2]);
  if (start>end){
    printf("Invaid start (%d), end(%d): start > end\n", start,end);
    exit(EXIT_FAILURE);
  }
  if (start<0){
    printf("Invaid start (%d): start < 0\n", start);
    exit(EXIT_FAILURE);
  }

  printf("This calculation began at: %s", ctime(&currentTime)); 
  fflush(stdout);



  /*** read initial setup from input.par ***/
  init_size(argc, argv, &par, &ist);

  /*** allocating memory ***/
  /*** the positions of the atoms in the x, y, and z directions ***/
  if ((rx = (double *) calloc(ist.natom, sizeof(double))) == NULL) nerror("rx");
  if ((ry = (double *) calloc(ist.natom, sizeof(double))) == NULL) nerror("ry");
  if ((rz = (double *) calloc(ist.natom, sizeof(double))) == NULL) nerror("rz");
  if ((atm = (atm_st *) calloc(ist.natom, sizeof(atm_st))) == NULL) nerror("atm");
  init_conf(rx, ry, rz, atm, &par, &ist);
	
  
  if ((rho = (double *)calloc(ist.ngrid,sizeof(double)))==NULL) nerror("rho");


  //count number of states found
  jms = countlines("eval.dat");
  printf("%ld total states in psi.dat\n", jms);
  
  //allocate memory for psitot
  if ((psitot = (zomplex *) calloc(ist.nspinngrid, sizeof(zomplex))) == NULL) nerror("psitot");


  //read psi from file
	ppsi = fopen("psi.dat" , "r");

	
  
  char filename[20];
  for (int j = start; j<=end; j++){ 
    printf("Reading state %d from psi.dat\n", j);
    if(fseek(ppsi,j*ist.nspinngrid*sizeof(zomplex),SEEK_SET)!=0)
    {
      printf("Error reading from psi.dat!\n"); exit(EXIT_FAILURE);
    }
    fread (&psitot[0],sizeof(zomplex),ist.nspinngrid,ppsi);

    for (long i = 0;i<ist.ngrid; i++){
      rho[i] = sqr(psitot[i].re)+sqr(psitot[i].im);
              
    }
    sprintf(filename, "rhoUp%i.cube", j);
    writeCubeFile(rho, par,ist, filename);
    
    for (long i = 0;i<ist.ngrid; i++){
      rho[i]=sqr(psitot[ist.ngrid+i].re)+sqr(psitot[ist.ngrid+i].im);
    }
    sprintf(filename, "rhoDn%i.cube", j);
    writeCubeFile(rho, par,ist, filename);

    for (long i = 0;i<ist.ngrid; i++){
      rho[i]= sqr(psitot[i].re)+sqr(psitot[i].im)+
              sqr(psitot[ist.ngrid+i].re)+sqr(psitot[ist.ngrid+i].im);
    }
    sprintf(filename, "rhoTot%i.cube", j);
    writeCubeFile(rho, par,ist, filename);

  }
  fclose(ppsi);  

  return 0;


}


/*****************************************************************************/
int countlines(char *filename){
  FILE* fp = fopen(filename,"r");
  int lines = 0;
  int ch;
  while(1){
    ch = fgetc(fp);
    if (feof(fp)){ break; }
    if(ch == '\n')
    {
      lines++;
    }
  }
  fclose(fp);
  return lines;
}


