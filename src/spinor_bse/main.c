/*****************************************************************************/

#include "fd.h"
#include <float.h>

/*****************************************************************************/

int main(int argc, char *argv[]) {
    FILE *ppsi, *peval;  
    long i, a, j, thomo, tlumo, indexfirsthomo, flags=0;
    double *eval, *de, *ksqr, *vx, *vy, *vz,  *rx, *ry, *rz;
    double *potl, *h0mat, *mx, *my, *mz, *rs;
    zomplex *mux, *muy, *muz, *potq, *potqx, *poth,*psidummy, *psi, *bsmat, *direct, *exchange;
    zomplex *sx, *sy, *sz;
    zomplex *lx, *ly, *lz, *lsqr, *ls;
    fftw_plan_loc *planfw, *planbw; fftw_complex *fftwpsi;
    par_st par;
    long_st ist;

    /*************************************************************************/
    writeCurrentTime(stdout);
    fflush(0);
    writeSeparation(stdout);
    fflush(0);
    init_size(argc, argv, &par, &ist);
  
    /*************************************************************************/
    fftwpsi = fftw_malloc(sizeof (fftw_complex )*ist.ngrid*ist.nthreads);
    potl  = (double *) calloc(ist.ngrid, sizeof(double));
    potq  = (zomplex *) calloc(ist.ngrid, sizeof(zomplex));
    potqx  = (zomplex *) calloc(ist.ngrid, sizeof(zomplex));
    poth = (zomplex *) calloc(ist.ngrid*ist.nthreads, sizeof(zomplex));
    ksqr = (double *) calloc(ist.ngrid, sizeof(double));
    vx = (double *) calloc(ist.nx, sizeof(double));
    vy = (double *) calloc(ist.ny, sizeof(double));
    vz = (double *) calloc(ist.nz, sizeof(double));
    rx = (double *) calloc(ist.natom, sizeof(double));
    ry = (double *) calloc(ist.natom, sizeof(double));
    rz = (double *) calloc(ist.natom, sizeof(double));
  
    /**************************************************************************/
    //Initilize con figuration, grid, and ksqr grid
    init(potl, vx, vy, vz, ksqr, rx, ry, rz, &par, &ist);

    /*** initialization for the fast Fourier transform ***/
    fftw_plan_with_nthreads(ist.nthreads);
  
    planfw = (fftw_plan_loc *) calloc(ist.nthreads, sizeof(fftw_plan_loc));
    planbw = (fftw_plan_loc *) calloc(ist.nthreads, sizeof(fftw_plan_loc));
    for (i = 0; i < ist.nthreads; i++) { 
        planfw[i] = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, &fftwpsi[i*ist.ngrid], 
                                     &fftwpsi[i*ist.ngrid], FFTW_FORWARD, flags);
        planbw[i] = fftw_plan_dft_3d(ist.nz, ist.ny, ist.nx, &fftwpsi[i*ist.ngrid],
                                     &fftwpsi[i*ist.ngrid], FFTW_BACKWARD, flags);
    }
    init_pot(vx, vy, vz, potq, potqx, par, ist, planfw[0], planbw[0], &fftwpsi[0]);

    /*************************************************************************/
    eval = (double *) calloc(ist.nhomo+1, sizeof(double)); 
    de = (double *) calloc(ist.nhomo+1, sizeof(double)); 

    peval = fopen("eval.par" , "r");
    for (i = 0; i < ist.nhomo+1; i++)
        fscanf (peval,"%ld %lg %lg",&a,&eval[i],&de[i]);
    fclose (peval);

    peval = fopen("eval.dat" , "w");
    for (i = 0; i < ist.nhomo+1; i++){
        fprintf (peval,"%ld %g %g\n",i,eval[i],de[i]);
    }
    fclose(peval);
  
    for (thomo = 0, i = ist.nhomo; i >= 0; i--){
        if (de[i] < par.deps) thomo++;
        if (thomo == ist.totalhomo) {
            indexfirsthomo = i;
            break;
        }
    }

    printf("The index of lowest energy occupied level used = %ld\n", indexfirsthomo); 

    free(eval);
    free(de);

    psidummy = (zomplex*)calloc(ist.nspinngrid,sizeof(zomplex));
    eval = (double*)calloc(ist.ms,sizeof(double)); 
    de = (double*)calloc(ist.ms,sizeof(double)); 


    /**********************************************************************/
    /*               Read in eval.par and store it in memory              */
    /**********************************************************************/
    struct eindex { // Struct I used to store index evals and deps in memory
        long index;
        double evalue;
        double deps;
    }; 
    
    struct eindex *evalindex;
    int size = 1024; // Initial size of array in bytes, is allowed to resize
    long counter = 0;

    if ((evalindex = calloc(size, sizeof(struct eindex))) == NULL) {
        fprintf(stderr, "Memory error in main\n");
        exit(EXIT_FAILURE);
    }

    peval = fopen("eval.par" , "r");
    if (peval) {
         while (fscanf(peval, "%ld %lg %lg", &evalindex[counter].index, 
                &evalindex[counter].evalue, &evalindex[counter].deps) == 3) {
            counter++;
            if (counter == size - 1) {
                size *= 2;
                // TODO: This does not check if realloc succeeded or not
                evalindex = realloc(evalindex, size * sizeof(struct eindex));
            }
        }
    }
    fclose(peval);

    // Read ead in psi.par as we have already read in psi.par
    ppsi = fopen("psi.par" , "r");
    if (ppsi ==NULL) {
        printf("no psi.par in cwd\n");
        exit(EXIT_FAILURE);
    }
    printf("allocating memory for %ld hole and %ld electron states\n",(ist.totalhomo),(ist.totallumo));
    zomplex *psihomo = calloc((ist.totalhomo)*ist.nspinngrid, sizeof(zomplex));
    zomplex *psilumo = calloc((ist.totallumo)*ist.nspinngrid, sizeof(zomplex));
    if (!psihomo || !psilumo) {printf("Failed to allocate memory for psihomo/psilumo\n");exit(EXIT_FAILURE);}

    long foffset = ist.nspinngrid * sizeof(zomplex);  // for random access 
    char fname[80] = {0};
    long nstates = 0;  // Total number of states
    fseek(ppsi, foffset * ist.nhomo, SEEK_SET); 
    counter = ist.nhomo; // Set loop counter to HOMO index
    thomo = 0;           // Counter for hole states used
    while (counter >= indexfirsthomo && thomo < ist.maxHoleStates) {
        fread(psidummy, sizeof(zomplex), ist.nspinngrid, ppsi);
        if (evalindex[counter].deps < par.deps) {    
            normalize_zomplex(psidummy, par.dv, ist.nspinngrid);
            //sprintf(fname, "pzv%ld.dat", counter);
            eval[nstates] = evalindex[counter].evalue;
            de[nstates] = evalindex[counter].deps;
            for (j = 0; j < ist.nspinngrid; j++) {
                psihomo[thomo * ist.nspinngrid + j] = psidummy[j];
            }
            nstates++;
            thomo++;
			
        }
        counter--;
        
        int ret =0;
        if (counter >=0) {
            ret= fseek(ppsi, -2 * foffset, SEEK_CUR);
        }
        if(ret) printf("nonzero exit in fseek!\n"); fflush(0); 
        
    }

    fseek(ppsi, foffset * (ist.nhomo + 1), SEEK_SET);
    counter = ist.nhomo + 1;
    tlumo = 0;

    while (tlumo < ist.maxElecStates && tlumo < ist.totallumo) {
        fread(psidummy, sizeof(zomplex), ist.nspinngrid, ppsi);
        if (evalindex[counter].deps < par.deps) {
            normalize_zomplex(psidummy, par.dv, ist.nspinngrid);
            //sprintf(fname, "pzc%ld.dat", counter);
            
		    eval[nstates] = evalindex[counter].evalue;
            de[nstates] = evalindex[counter].deps;
            for (j = 0; j < ist.nspinngrid; j++) {      
                psilumo[tlumo * ist.nspinngrid + j] = psidummy[j];
            }
        	nstates++;
        	tlumo++;
			
        }
        counter++;

    }
    free(psidummy);  printf("freeing psidummy\n");fflush(0);
    psi = (zomplex *) calloc(nstates, ist.nspinngrid*sizeof(zomplex));
    if (!psi) {printf("Failed to allocate memory for psi"); fflush(0); exit(EXIT_FAILURE);}
    fclose(ppsi); printf("closed ppsi");fflush(0);
    
    /**********************************************************************/

    for (i = thomo - 1, a = 0; i >= 0 && a < thomo; i--, a++) {
        for (j = 0; j < ist.nspinngrid; j++) {
            psi[a * ist.nspinngrid + j].re = psihomo[i * ist.nspinngrid + j].re;
            psi[a * ist.nspinngrid + j].im = psihomo[i * ist.nspinngrid + j].im;
        }
    }

    for (i = thomo, a = 0; i < tlumo + thomo && a < tlumo; i++, a++) {
        for (j = 0; j < ist.nspinngrid; j++) {
            psi[i * ist.nspinngrid + j].re = psilumo[a * ist.nspinngrid + j].re;
            psi[i * ist.nspinngrid + j].im = psilumo[a * ist.nspinngrid + j].im;
        }
    }

    double tmp1, tmp2;
    /* Rearrange eval */
    for (i = 0, j = thomo - 1; i < j; i++, j--) {
        tmp1 = eval[i];    tmp2 = de[i];
        eval[i] = eval[j]; de[i] = de[j];
        eval[j] = tmp1;    de[j] = tmp2;
    }
    /* Write smaller psi.par and eval.par for augerBSE */
    FILE *new_eval = fopen("BSEeval.par", "w");
    FILE *new_psi = fopen("BSEpsi.par", "w");
    if (!new_eval || !new_psi) {printf("Failed opening BSEeval.par");fflush(0); exit(EXIT_FAILURE);}

    fwrite(psi, sizeof(zomplex) * ist.nspinngrid, nstates, new_psi);

    for (i = 0; i < nstates; i++) {
        fprintf(new_eval, "%ld %.*g %.*g\n", i, DBL_DIG, eval[i], DBL_DIG, de[i]); 
    }
    fclose(new_eval);
    fclose(new_psi);
    /**********************************************************************/
    
    
    free(evalindex);
    free(psihomo);
    free(psilumo);
    

    printf("Number of hole eigenstates used in the BSE calculation = %ld\n", thomo);
    printf("Number of electron eigenstates used in the BSE calculation =  %ld\n", tlumo);
    printf("Total number of eigenstates used in the BSE calculation = %ld\n", nstates);
    fflush(0);
    for (int state  = 0; state<nstates; state++){
        printf("\nStats on state%d: (E=%lg)\n",state, eval[state]);
        double perUp = 0;
        double perDn = 0;
        for (long igrid = 0; igrid < ist.ngrid; igrid++){
            perUp+=sqr(psi[state*ist.nspinngrid+igrid].re)+sqr(psi[state*ist.nspinngrid+igrid].im);
            perDn+=sqr(psi[state*ist.nspinngrid+ist.ngrid+igrid].re)+sqr(psi[state*ist.nspinngrid+ist.ngrid+igrid].im);
        }
        printf(" Spin up fraction: %f\n", perUp*par.dv);
        printf(" Spin dn fraction: %f\n", perDn*par.dv);
    }

    /*************************************************************************/
    ist.ms = nstates;
    ist.nlumo = thomo;
    ist.nhomo = thomo - 1;
    ist.totallumo = tlumo;
    ist.totalhomo = thomo;
  
    
    /*** this routine computes the coulomb coupling between
         single excitons.  On input - it requires the eigenstates stored in psi,
         the eigenvalues stored in eval and poth computed in init_pot.
         On output it stores the coulomb matrix elements on the disk
         in the following format: a, i, b, j, ene_ai, ene_bj, vjbai, vabji.
         a - the index of the electron in exciton Sai.
         i - the index of the hole in exciton Sai.
         b - the index of the electron in exciton Sbj.
         j - the index of the hole in exciton Sbj.
         ene_ai - the energy of exciton Sai.
         ene_bj - the energy of exciton Sbj.
         vjbai and vabji are the coulomb matrix elements need to be used to
         generate the spin-depedent matrix elements as described by
         the last equation in our codument.  ***/

    ist.ms2 = ist.totallumo * ist.totalhomo;
    bsmat = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex));
    direct = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex)); 
    exchange = (zomplex *) calloc(ist.ms2*ist.ms2, sizeof(zomplex)); 
    h0mat = (double *) calloc(ist.ms2*ist.ms2, sizeof(double)); 
    
    single_coulomb_openmp(psi, potq, potqx, poth, eval, ist, par, planfw, planbw, fftwpsi, bsmat,direct,exchange, h0mat);

    ppsi = fopen("bsRE.dat", "w");
    for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
        for (j = 0; j < ist.ms2; j++)
            fprintf(ppsi,"%.*g\t", DBL_DIG, bsmat[i*ist.ms2+j].re);
    fclose(ppsi);

    ppsi = fopen("bsIM.dat", "w");
    for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
        for (j = 0; j < ist.ms2; j++)
            fprintf(ppsi,"%.*g\t", DBL_DIG, bsmat[i*ist.ms2+j].im);
    fclose(ppsi);

    ppsi = fopen("h0.dat", "w");
    for (i = 0; i < ist.ms2; i++, fprintf(ppsi,"\n"))
        for (j = 0; j < ist.ms2; j++)
             fprintf(ppsi,"%.*g\t", DBL_DIG, h0mat[i*ist.ms2+j]);
    fclose(ppsi);


    sx = (zomplex *) calloc(ist.totalhomo*ist.totalhomo+ist.totallumo*ist.totallumo, sizeof(zomplex)); //<psi_r|Sx|psi_s>
    sy = (zomplex *) calloc(ist.totalhomo*ist.totalhomo+ist.totallumo*ist.totallumo, sizeof(zomplex)); //<psi_r|Sy|psi_s>
    sz = (zomplex *) calloc(ist.totalhomo*ist.totalhomo+ist.totallumo*ist.totallumo, sizeof(zomplex)); //<psi_r|Sz|psi_s>
    lx = (zomplex *) calloc(ist.totalhomo*ist.totalhomo+ist.totallumo*ist.totallumo, sizeof(zomplex)); //<psi_r|Lx|psi_s>
    ly = (zomplex *) calloc(ist.totalhomo*ist.totalhomo+ist.totallumo*ist.totallumo, sizeof(zomplex)); //<psi_r|Ly|psi_s>
    lz = (zomplex *) calloc(ist.totalhomo*ist.totalhomo+ist.totallumo*ist.totallumo, sizeof(zomplex)); //<psi_r|Lz|psi_s>
    lsqr = (zomplex *) calloc(ist.totalhomo*ist.totalhomo+ist.totallumo*ist.totallumo, sizeof(zomplex)); //<psi_r|L^2|psi_s>
    ls  = (zomplex *) calloc(ist.totalhomo*ist.totalhomo+ist.totallumo*ist.totallumo, sizeof(zomplex)); //<psi_r|L*S|psi_s>
    mux = (zomplex *) calloc(ist.totallumo*ist.totalhomo, sizeof(zomplex)); // <psi_i|ux|psi_a>
    muy = (zomplex *) calloc(ist.totallumo*ist.totalhomo, sizeof(zomplex)); // <psi_i|uy|psi_a>
    muz = (zomplex *) calloc(ist.totallumo*ist.totalhomo, sizeof(zomplex)); // <psi_i|uz|psi_a>
    
    spins(sx,sy,sz,psi,ist,par);
    angular(lx,ly,lz,lsqr,ls,vx,vy,vz,psi,planfw[0], planbw[0], &fftwpsi[0],ist,par);
    dipole(vx, vy, vz, psi, mux, muy, muz, eval, ist, par);
    bethe_salpeter(bsmat, direct, exchange, h0mat, psi, vz, mux, muy, muz, mx, my, mz,sx,sy,sz,lx,ly,lz,lsqr,ls, ist, par);
  
    writeSeparation(stdout);
    writeCurrentTime(stdout);

    /***********************************************************************/
    free(psi); free(potq); free(potqx); free(eval);  free(de); 
    free(ksqr); free(vx); free(vy);  free(vz); 
    free(rx); free(ry); free(rz); free(poth);  free(potl);
    free(bsmat); free(h0mat); free(mux); free(muy); free(muz);
    //free(mx); free(my); free(mz); free(rs);
    free(sx); free(sy); free(sz);
    free(lx); free(ly); free(lz); free(lsqr);
    free(planfw); free(planbw);

    return 0;
}

/*****************************************************************************/
