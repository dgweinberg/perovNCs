#include "fd.h"

/*****************************************************************************/

void read_conf(double *rx,double *ry,double *rz,atm_st *atm,long ntot,FILE *pf, long_st *ist)
{
  long i,j, nnonlocal; FILE *pw;
  double xd, yd, zd;

  for (j=0;j<20;j++){
    ist->atomspresent[j] = -1;
  }
  ist->natomtype = 0;
  for (nnonlocal = 0, xd = yd = zd = 0.0, i = 0; i < ntot; i++){
    fscanf (pf,"%s %lf %lf %lf",atm[i].atyp,&rx[i],&ry[i],&rz[i]);
    atm[i].Zval = assign_atom_number(atm[i].atyp);
    //update list of atom numbers in ist->atoms
    for (j=0;j<=ist->natomtype;j++){
      //if atom already in list break
      if (atm[i].Zval==ist->atomspresent[j]) {
        atm[i].natyp = j;
        break; 
      }
      //else check if we can write this position
      else if (j==ist->natomtype){
        ist->atomspresent[j] = atm[i].Zval;
        ist->natomtype++;
        printf("natomtype = %ld\n", ist->natomtype);
        atm[i].natyp = j;
        printf("Found new atom with Zval: %d\n",atm[i].Zval);
        printf("The atoms so far found are [ ");
        for (int k =0; k<ist->natomtype;k++){
          printf("%ld ",ist->atomspresent[k]);
        }
        printf("]\n");
        break;
      }
      else if (j==19){
        printf("Too mandy distinct atom types! Exiting..."); fflush(0);
        exit(EXIT_FAILURE);
      }
    }
    
    /*
      Cd = 0
      Se = 1
      In = 2
      As = 3
      Si = 4
      H  = 5
      Zn = 6
      S  = 7
      P1 = 8
      P2 = 9
      P3 = 10
      P4 = 11
      Te = 12
      Cdz = 13
      Sez = 14
      Ga  = 15
  */

    atm[i].Vso = 0.0;
    if (atm[i].Zval == 0) atm[i].Vso = 0.0;
    if (atm[i].Zval == 1) atm[i].Vso = 0.0;
    
    //In
    //if (atm[i].natyp == 2) atm[i].Vso = 0.48 * cube(6.057/0.529177) / 8.0 / 27.2114;
    //if (atm[i].natyp == 2) atm[i].Vso = 0.00876822726*2.0;
    if (atm[i].Zval == 2) atm[i].Vso = 3.42974738;
    
    //As
    //if (atm[i].natyp == 3) atm[i].Vso = 0.0976 * cube(6.057/0.529177) / 8.0 / 27.2114;
    //if (atm[i].natyp == 3) atm[i].Vso = 0.00223915662*2.0;
    if (atm[i].Zval == 3) atm[i].Vso = 0.78183193;
    
    //Si
    if (atm[i].Zval == 4) atm[i].Vso = 0.0;
    //H
    if (atm[i].Zval == 5) atm[i].Vso = 0.0;
    if (atm[i].Zval == 6) atm[i].Vso = 0.0;
    if (atm[i].Zval == 7) atm[i].Vso = 0.0;
    if (atm[i].Zval == 12) atm[i].Vso = 0.0;
    if (atm[i].Zval == 13) atm[i].Vso = 0.0;
    if (atm[i].Zval == 14) atm[i].Vso = 0.0;
    
    
    /*
    if (atm[i].Zval == 15) atm[i].Vso = 1.15454347; // Ga
    if (atm[i].Zval == 3) atm[i].Vso = 3.15164354; // As
    */

    




    if (atm[i].Zval < 8 || atm[i].Zval > 11) nnonlocal++;
    xd += rx[i];
    yd += ry[i];
    zd += rz[i];
  }
  xd /= (double)(ntot);
  yd /= (double)(ntot);
  zd /= (double)(ntot);
  if (ist->flagCenter){
    for (i = 0; i < ntot; i++){
      rx[i] -= xd;
      ry[i] -= yd;
      rz[i] -= zd;
    }
  }



  setAtmStr(rx, ry, rz, atm, ntot, ist);
  


  //CsPbI specific hard coded :'( 
  double CsSO[2], ISO[2], PbSO[2];
  double CsNL1[2], INL1[2], PbNL1[2];
  double CsNL2[2], INL2[2], PbNL2[2];
  if((pw = fopen("Cs_ortho_SO.par","r"))==NULL) nerror("CsOSO");
  fscanf(pw,"%lg %*g", &CsSO[0]);
  fscanf(pw,"%lg %*g", &CsNL1[0]);
  fscanf(pw,"%lg %*g", &CsNL2[0]);
  fclose(pw);

  if((pw = fopen("I_ortho_SO.par","r"))==NULL) nerror("IOSO");
  fscanf(pw,"%lg %*g", &ISO[0]);
  fscanf(pw,"%lg %*g", &INL1[0]);
  fscanf(pw,"%lg %*g", &INL2[0]);
  fclose(pw);

  if((pw = fopen("Pb_ortho_SO.par","r"))==NULL) nerror("PbOSO");
  fscanf(pw,"%lg %*g", &PbSO[0]);
  fscanf(pw,"%lg %*g", &PbNL1[0]);
  fscanf(pw,"%lg %*g", &PbNL2[0]);
  fclose(pw);


  if((pw = fopen("Cs_cubic_SO.par","r"))==NULL) nerror("CsCSO");
  fscanf(pw,"%lg %*g", &CsSO[1]);
  fscanf(pw,"%lg %*g", &CsNL1[1]);
  fscanf(pw,"%lg %*g", &CsNL2[1]);
  fclose(pw);

  if((pw = fopen("I_cubic_SO.par","r"))==NULL) nerror("ICSO");
  fscanf(pw,"%lg %*g", &ISO[1]);
  fscanf(pw,"%lg %*g", &INL1[1]);
  fscanf(pw,"%lg %*g", &INL2[1]);
  fclose(pw);

  if((pw = fopen("Pb_cubic_SO.par","r"))==NULL) nerror("PbCSO");
  fscanf(pw,"%lg %*g", &PbSO[1]);
  fscanf(pw,"%lg %*g", &PbNL1[1]);
  fscanf(pw,"%lg %*g", &PbNL2[1]);
  fclose(pw);


  for(i = 0; i < ntot; i++) {
    if (atm[i].Zval == 53){//I
      atm[i].Vso = atm[i].strPar*ISO[0]+(1.0-atm[i].strPar)*ISO[1];
      atm[i].nlcPar[0] = atm[i].strPar*INL1[0]+(1.0-atm[i].strPar)*INL1[1];
      atm[i].nlcPar[1] = atm[i].strPar*INL2[0]+(1.0-atm[i].strPar)*INL2[1];
    } 
    if (atm[i].Zval == 55) {  //Cs
      atm[i].Vso = atm[i].strPar*CsSO[0]+(1.0-atm[i].strPar)*CsSO[1];
      atm[i].nlcPar[0] = atm[i].strPar*CsNL1[0]+(1.0-atm[i].strPar)*CsNL1[1];
      atm[i].nlcPar[1] = atm[i].strPar*CsNL2[0]+(1.0-atm[i].strPar)*CsNL2[1];
    }
    if (atm[i].Zval == 82) {//Pb
      atm[i].Vso = atm[i].strPar*PbSO[0]+(1.0-atm[i].strPar)*PbSO[1];
      atm[i].nlcPar[0] = atm[i].strPar*PbNL1[0]+(1.0-atm[i].strPar)*PbNL1[1];
      atm[i].nlcPar[1] = atm[i].strPar*PbNL2[0]+(1.0-atm[i].strPar)*PbNL2[1];
    }
  } 




  pw = fopen("conf.dat" , "w");
  fprintf(pw,"%ld\n", ntot);
  for (i = 0; i < ntot; i++) {
    fprintf(pw, "%s %g %g %g %d %g %g\n", atm[i].atyp, rx[i], ry[i], rz[i], atm[i].Zval, atm[i].Vso,atm[i].strPar);
  }
  fclose(pw);
  
  
  ist->nnonlocal = nnonlocal;

  return;
}

/*****************************************************************************/
void setAtmStr(double *rx,double *ry,double *rz,atm_st *atm,long ntot,long_st *ist ){
  long j;
  double distSqr,bondAngle, aveIpar;
  long bonded[10];
  int nbonded;

  double PbIBondMax = 3.3*ANGTOBOHR; 
  //double orthoBondAngle =150.8431;
  double orthoBondAngle = 160.6344; //using the larger of the two ortho angles

  

  //first sweep only look for I and Cs
  for (long i = 0;i<ntot;i++ ){
  //check what type of atom this is. (Currently only do Pb, I, Cs)
    
    // for an I atom
    if (atm[i].Zval == 53) { 
      nbonded = 0;

      //get list of bonded lead atoms
      for(j=0; j<ntot; j++){
        if (atm[j].Zval == 82){
          distSqr = sqr(rx[i]-rx[j])+ sqr(ry[i]-ry[j])+ sqr(rz[i]-rz[j]);
          if (distSqr<sqr(PbIBondMax)){
            bonded[nbonded]=j;
            nbonded++;
          }
        }
      }


      //do something to calculate the bond angle....
      if(nbonded!=2){
        atm[i].strPar=0.0; //default to cubic if edge atom
        printf("atom %ld (I) is an edge atom (nbonded=%d)\n",i, nbonded);
      } 
      else{
        bondAngle = getBondAngle(bonded[0],i,bonded[1],rx,ry,rz);
        printf("atom %ld (I) has bond angle %g\n",i, bondAngle);

        if (bondAngle < orthoBondAngle) {
          atm[i].strPar=1.0;
        }
        else{
          atm[i].strPar = 1.0 - ((bondAngle - orthoBondAngle) / (180.0 - orthoBondAngle));
        }
      }

    } 

    

    //for a Cs
    else if (atm[i].Zval == 55) { 
      atm[i].strPar=0.0; //just ignore the Cs ... 
      printf("atom %ld (Cs)\n",i);
    }
    else if (atm[i].Zval == 82) continue;
    else{
      printf("Unknown atom with Zval %d!\nExiting...", atm[i].Zval);
      fflush(0);
      exit(EXIT_FAILURE);
    } 

  }
    

  //second sweep only look for Pb
  for (long i = 0;i<ntot;i++ ){
    
    if (atm[i].Zval == 53) continue;
    else if (atm[i].Zval == 55) continue;

    //for a Pb
    else if (atm[i].Zval == 82) {
      nbonded = 0;

      //get list of bonded I atoms
      for(j=0; j<ntot; j++){
        if (atm[j].Zval == 53){
          distSqr = sqr(rx[i]-rx[j])+ sqr(ry[i]-ry[j])+ sqr(rz[i]-rz[j]);
          if (distSqr<sqr(PbIBondMax)){
            bonded[nbonded]=j;
            nbonded++;
          }
        }
      }

      //go based on average bond length
      if(nbonded!=6){
        atm[i].strPar=0.0; //default to cubic if edge atom /
        printf("atom %ld (Pb) is an edge atom (nbonded=%d)\n",i, nbonded);
      } 
      else{
        aveIpar=0.0;
        for(j=0;j<6;j++){
          aveIpar+=atm[bonded[j]].strPar;
        }
        aveIpar/=6.0;

        printf("atom %ld (Pb) has neighbors with pars {%g, %g, %g, %g, %g, %g} (ave: %g) \n",i, atm[bonded[0]].strPar, atm[bonded[1]].strPar,atm[bonded[2]].strPar,
                                                                                        atm[bonded[3]].strPar, atm[bonded[4]].strPar, atm[bonded[5]].strPar,
                                                                                        aveIpar);

        atm[i].strPar = aveIpar;
      }
    }


    //for an unknown atom
    else{
      printf("Unknown atom with Zval %d!\nExiting...", atm[i].Zval);
      fflush(0);
      exit(EXIT_FAILURE);
    }
  } 
  return;

}

/*****************************************************************************/
double getBondAngle(long index1,long index2,long index3,double *rx,double *ry,double *rz){
  double dx1 = rx[index1]-rx[index2];
  double dy1 = ry[index1]-ry[index2];
  double dz1 = rz[index1]-rz[index2]; 
  double l1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);

  double dx2 = rx[index3]-rx[index2];
  double dy2 = ry[index3]-ry[index2];
  double dz2 = rz[index3]-rz[index2];
  double l2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2); 

  double dot = dx1*dx2 + dy1*dy2 + dz1*dz2;
  if (-1.0>dot/(l1*l2)||1.0<dot/(l1*l2)){return 180.0;}  
  if (-1.0>dot/(l1*l2)||1.0<dot/(l1*l2)){
    printf("Error in Bond Angle (%lf)!!\n", dot/(l1*l2));
    printf("1: %f %f %f\n", rx[index1], ry[index1],rz[index1] );
    printf("2: %f %f %f\n", rx[index2], ry[index2],rz[index2] );
    printf("3: %f %f %f\n", rx[index3], ry[index3],rz[index3] );
  }

  return acos(dot/(l1*l2))*180.0/PIE;

}


/*****************************************************************************/

void read_pot(double *vr,double *pot,long *npot,double *dr,atm_st *atm,long_st *ist)
{
  FILE *pf;  long i, j, iscan; char str[100], atype[3];
  long n  = ist->npot;
  long ntype = ist->natomtype;
  /*double *a, *b;

  if ((a = (double*)calloc(ntype,sizeof(double)))==NULL)nerror("a");
  if ((b = (double*)calloc(ntype,sizeof(double)))==NULL)nerror("b");

  a[8] = 0.64;
  a[9] = -0.384;
  a[10] = 0.04;
    a[11] = -0.684;
  a[10] = 0.64;
  a[11] = -0.384;
  b[8] = b[9] = 2.2287033;
  b[10] = 2.2287033;
  b[11] = 2.2287033;
  */
  //double a[4] = { 0.64, -0.384, 0.04, -0.684};
  //double b[4] = {2.2287033, 2.2287033, 2.2287033, 2.2287033};

  /*
      Cd = 0
      Se = 1
      In = 2
      As = 3
      Si = 4
      H  = 5
      Zn = 6
      S  = 7
      P1 = 8
      P2 = 9
      P3 = 10
      P4 = 11
      Te = 12
      Cdz = 13
      Sez = 14
      Ga  = 15
      
      I = 53
      Cs = 55
      Pb = 82



  */

  for (j = 0; j < 2*ntype*n; j++) pot[j] = vr[j] = 0;
  for (j = 0; j < 2*ntype; j++) npot[j] = 0;
  
  for (j = 0; j < ntype; j++){
    long iatm  = ist->atomspresent[j];
    assign_atom_type(atype,iatm);
    

  //read cubic Pot
    sprintf (str, "Cpot");
    for (int strnum = 0;strnum<3; strnum++){  
      if(atype[strnum]=='\0') break;
      strncat(str,&atype[strnum],1);
    }
    strcat(str, ".par");
 
    pf = fopen(str , "r");
    if (pf != NULL) {
      for (npot[2*j] = iscan = 0; iscan != EOF; npot[2*j]++) {
        iscan = fscanf (pf,"%lg %lg",&vr[j*2*n+npot[2*j]],&pot[j*2*n+npot[2*j]]);
      }
      fclose(pf);
      npot[2*j]--;
    }      
    else {
      printf("Unable to open pot file %s to read! Exiting...\n", str); fflush(0);
      exit(EXIT_FAILURE);
    }   
    
    /*** shift the potentials and get the r-spacing ***/
    for (i = 0; i < npot[2*j]; i++) {
      pot[j*2*n+i] -= pot[j*2*n+npot[2*j]-1];
      dr[2*j] = vr[j*2*n+1] - vr[j*2*n+0];
    }
      
    /*** print shifted pot ***/
  
    sprintf (str, "Cpot");
    for (int strnum = 0;strnum<3; strnum++){
      if(atype[strnum]=='\0') break;
      strncat(str,&atype[strnum],1);
    }
    strcat(str, ".dat");
    printf("%s",str);

    pf = fopen(str , "w");
    for (i = 0; i < npot[2*j]; i++) fprintf (pf,"%g %g\n",vr[j*2*n+i],pot[j*2*n+i]);
    fclose(pf);

  //read ortho Pot
    sprintf (str, "Opot");
    for (int strnum = 0;strnum<3; strnum++){  
      if(atype[strnum]=='\0') break;
      strncat(str,&atype[strnum],1);
    }
    strcat(str, ".par");
 
    pf = fopen(str , "r");
    if (pf != NULL) {
      for (npot[2*j+1] = iscan = 0; iscan != EOF; npot[2*j+1]++) {
        iscan = fscanf (pf,"%lg %lg",&vr[j*2*n+n+npot[2*j+1]],&pot[j*2*n+n+npot[2*j+1]]);
      }
      fclose(pf);
      npot[2*j+1]--;
    }      
    else {
      printf("Unable to open pot file %s to read! Exiting...\n", str); fflush(0);
      exit(EXIT_FAILURE);
    }   
    
    /*** shift the potentials and get the r-spacing ***/
    for (i = 0; i < npot[2*j+1]; i++) {
      pot[j*2*n+n+i] -= pot[j*2*n+n+npot[2*j+1]-1];
      dr[2*j+1] = vr[j*2*n+n+1] - vr[j*2*n+n+0];
    }
      
    /*** print shifted pot ***/
  
    sprintf (str, "Opot");
    for (int strnum = 0;strnum<3; strnum++){
      if(atype[strnum]=='\0') break;
      strncat(str,&atype[strnum],1);
    }
    strcat(str, ".dat");
    printf("%s",str);

    pf = fopen(str , "w");
    for (i = 0; i < npot[2*j+1]; i++) fprintf (pf,"%g %g\n",vr[j*2*n+n+i],pot[j*2*n+n+i]);
    fclose(pf);


  }
  //free(a); free(b);
  return;
}

/*****************************************************************************/

long assign_atom_number(char atyp[3])
{
  char strerror[100];
  
  if ((atyp[0] == 'C') && (atyp[1] == 'd')  && (atyp[2] == '\0')) return(0);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(1);
  else if ((atyp[0] == 'I') && (atyp[1] == 'n') && (atyp[2] == '\0')) return(2);
  else if ((atyp[0] == 'A') && (atyp[1] == 's') && (atyp[2] == '\0')) return(3);
  else if ((atyp[0] == 'S') && (atyp[1] == 'i') && (atyp[2] == '\0')) return(4);
  else if ((atyp[0] == 'H') && (atyp[1] == '\0') && (atyp[2] == '\0'))  return(5);
  else if ((atyp[0] == 'Z') && (atyp[1] == 'n') && (atyp[2] == '\0'))  return(6);
  else if ((atyp[0] == 'S') && (atyp[1] == '\0') && (atyp[2] == '\0'))  return(7);
  else if ((atyp[0] == 'P') && (atyp[1] == '1') && (atyp[2] == '\0'))  return(8);
  else if ((atyp[0] == 'P') && (atyp[1] == '2') && (atyp[2] == '\0'))  return(9);
  else if ((atyp[0] == 'P') && (atyp[1] == '3') && (atyp[2] == '\0'))  return(10);
  else if ((atyp[0] == 'P') && (atyp[1] == '4') && (atyp[2] == '\0'))  return(11);
  else if ((atyp[0] == 'T') && (atyp[1] == 'e') && (atyp[2] == '\0')) return(12);
  else if ((atyp[0] == 'C') && (atyp[1] == 'd') && (atyp[2] == 'z')) return(13);
  else if ((atyp[0] == 'S') && (atyp[1] == 'e') && (atyp[2] == 'z')) return(14);
  else if ((atyp[0] == 'G') && (atyp[1] == 'a') && (atyp[2] == '\0')) return(15);
  else if ((atyp[0] == 'I') && (atyp[1] == '\0') && (atyp[2] == '\0')) return(53);
  else if ((atyp[0] == 'C') && (atyp[1] == 's') && (atyp[2] == '\0')) return(55);
  else if ((atyp[0] == 'P') && (atyp[1] == 'b') && (atyp[2] == '\0')) return(82);
  else if ((atyp[0] == 'X') && (atyp[1] == '\0') && (atyp[2] == '\0')) return(-1);
  
  else {
    sprintf (strerror,"atom type %s not in current list",atyp);
    nerror (strerror);
  }
  return(0);
}

/*****************************************************************************/

void assign_atom_type(char *atyp,long j)
{
  if (j == 0) {atyp[0] = 'C'; atyp[1] = 'd'; atyp[2] = '\0';}
  else if (j == 1) {atyp[0] = 'S'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 2) {atyp[0] = 'I'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 3) {atyp[0] = 'A'; atyp[1] = 's'; atyp[2] = '\0';}
  else if (j == 4) {atyp[0] = 'S'; atyp[1] = 'i'; atyp[2] = '\0';}
  else if (j == 5) {atyp[0] = 'H'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 6) {atyp[0] = 'Z'; atyp[1] = 'n'; atyp[2] = '\0';}
  else if (j == 7) {atyp[0] = 'S'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 8) {atyp[0] = 'P'; atyp[1] = '1'; atyp[2] = '\0';}
  else if (j == 9) {atyp[0] = 'P'; atyp[1] = '2'; atyp[2] = '\0';}
  else if (j == 10) {atyp[0] = 'P'; atyp[1] = '3'; atyp[2] = '\0';}
  else if (j == 11) {atyp[0] = 'P'; atyp[1] = '4'; atyp[2] = '\0';}
  else if (j == 12) {atyp[0] = 'T'; atyp[1] = 'e'; atyp[2] = '\0';}
  else if (j == 13) {atyp[0] = 'C'; atyp[1] = 'd'; atyp[2] = 'z';}
  else if (j == 14) {atyp[0] = 'S'; atyp[1] = 'e'; atyp[2] = 'z';}
  else if (j == 15) {atyp[0] = 'G'; atyp[1] = 'a'; atyp[2] = '\0';}
  else if (j == 53) {atyp[0] = 'I'; atyp[1] = '\0'; atyp[2] = '\0';}
  else if (j == 55) {atyp[0] = 'C'; atyp[1] = 's'; atyp[2] = '\0';}
  else if (j == 82) {atyp[0] = 'P'; atyp[1] = 'b'; atyp[2] = '\0';}
  else if (j == -1) {atyp[0] = 'X'; atyp[1] = '\0'; atyp[2] = '\0';}

  return;
}

/*****************************************************************************/

long get_number_of_atom_types(atm_st *atm,long_st ist,long *list)
{
  long jatom, jlist, nlist, flag;

  for (jlist = 0; jlist < 200; jlist++) list[jlist] = -1;
  list[0] = 8; list[1] = 9; list[2] = 10; list[3] = 11;
  for (nlist = 4, jatom = 0; jatom < ist.natom; jatom++){
    for (flag = 0, jlist = 0; jlist < nlist; jlist++)
      if (atm[jatom].natyp == list[jlist]) flag = 1;
    if (flag != 1) {list[nlist] = atm[jatom].natyp; nlist++;}
  }
  //for (jlist = 0; jlist < nlist-4; jlist++) printf ("atom in list %ld\n",list[jlist+4]);
  return(nlist-4);
}

/*****************************************************************************/
