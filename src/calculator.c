/* FILE: calculator.c
 * AUTH: Maxim Yurkin
 * DESCR: All the initialization is done here before actually
 *        calculating internal fields. Calculation of couple constants.
 *
 *        Previous versions were by Alfons Hoekstra
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include <time.h>
#include "vars.h"
#include "cmplx.h"
#include "Romberg.h"
#include "const.h"
#include "comm.h"
#include "debug.h"
#include "memory.h"
#include "crosssec.h"
#include "io.h"
#include "fft.h"

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in crosssec.c */
extern Parms_1D parms[2],parms_alpha;
extern angle_set beta_int,gamma_int,theta_int,phi_int;
/* defined and initialized in io.c */
extern int PolRelation,avg_inc_pol;
extern char alldir_parms[],scat_grid_parms[];
/* defined and initialized in timing.c */
extern clock_t Timing_Init, tstart_main;
extern unsigned long TotalEval;

/* used in CalculateE.c */
double *muel_phi;  /* used to store values of Mueller matrix for different phi (to integrate) */
double *muel_phi1; /* additional for integrating with different multipliers */
doublecomplex *EplaneX, *EplaneY; /* scattered E (for scattering in one plane)
                                     for two incident polarizations */
double *Eplane_buffer;            /* buffer to accumulate Eplane */
double dtheta_deg,dtheta_rad;     /* delta theta in deg and radians */
doublecomplex *ampl_alphaX,*ampl_alphaY; /* storing amplitude matrix for
                                            different values of alpha */
double *muel_alpha;  /* storing mueller matrix for different values of alpha */
/* used in fft.c */
double *tab1,*tab2,*tab3,*tab4,*tab5,*tab6,*tab7,*tab8,*tab9,*tab10; /* tables of integrals */
int **tab_index;      /* matrix for indexation of table arrays */
/* used in crosssec.c */
double *E2_alldir;  /* square of E, calculated for alldir */
double *E2_alldir_buffer; /* buffer to accumulate E2_alldir */
doublecomplex cc[MAX_NMAT][3];         /* couple constants */
/* used in iterative.c */
doublecomplex *xvec,*rvec,*vec1,*vec2,*vec3,*Avecbuffer;  /* vectors for iterative solvers */

/* LOCAL VARIABLES */

int block_theta;       /* size of one block of mueller matrix - 16*nTheta */
int finish_avg;        /* whether to stop orientation averaging */
double *out;           /* used to collect both mueller matrix and integral
                          scattering quantities when orient_avg */
/* EXTERNAL FUNCTIONS */

/* CalculateE.c */
extern int CalculateE(char which,int type);
extern void MuellerMatrix(void);

/*============================================================*/

void CoupleConstant(doublecomplex mrel,char which,doublecomplex *res)
{
  doublecomplex CoupleConstant[3];
  doublecomplex tempa,tempb,cm,m2,t1;
  double temp,V,b1,b2,b3;
  int i,j;
  double S,prop2[3];
  int asym;           /* whether polarizability is asymmetric */
  double *incPol;

  asym = (PolRelation==POL_CLDR || PolRelation==POL_SO);
  if (asym) j=3;
  else j=1;
  if (PolRelation==POL_LDR || PolRelation==POL_CLDR) {
    b1=LDR_B1;
    b2=LDR_B2;
    b3=LDR_B3;
  }
  if (PolRelation==POL_SO) {
    b1=SO_B1;
    b2=SO_B2;
    b3=SO_B3;
  }
  /* calculate the CM couple constant CC=(3V/4pi)*(m^2-1)/(m^2+2) */
  V=gridspace*gridspace*gridspace;     /* volume of one dipole */
  temp = (3*V)/(4*PI);
  cSquare(mrel,m2);                    /* m2=m^2 */
  tempa[RE] = m2[RE] - 1.0;
  tempa[IM] = tempb[IM] = m2[IM];                      
  tempb[RE] = m2[RE] + 2.0;
  cDiv(tempa,tempb,CoupleConstant[0]);
  CoupleConstant[0][RE] *= temp;
  CoupleConstant[0][IM] *= temp;
  
  if (PolRelation!=POL_CM) {
    if (PolRelation!=POL_RR) {
        for (i=0;i<3;i++) prop2[i]=prop[i]*prop[i];
        if (PolRelation==POL_LDR) {
          if (avg_inc_pol) S=0.5*(1-DotProd(prop2,prop2));
          else {  
            if (which=='X') incPol=incPolX;
            else if (which=='Y') incPol=incPolY;
            S=prop2[0]*incPol[0]*incPol[0]+prop2[1]*incPol[1]*incPol[1]+prop2[2]*incPol[2]*incPol[2];
          }
       }
    }
    memcpy(cm,CoupleConstant[0],sizeof(doublecomplex));
    for (i=0;i<j;i++) {
      t1[RE]=0.0;                     
      t1[IM]=2*kd*kd*kd/3;                       /* t1=2/3*i*kd^3         */
      if (PolRelation!=POL_RR) {
        if (PolRelation!=POL_LDR) S=prop2[i];
        t1[RE]+=(b1+(b2+b3*S)*m2[RE])*kd*kd;   /* t1+=(b1+(b2+b3*S)*m^2)*kd^2  */
        t1[IM]+=(b2+b3*S)*m2[IM]*kd*kd;
      }
      t1[RE]/=V;     
      t1[IM]/=V;
      cMultSelf(t1,cm);
      t1[RE]=1-t1[RE];
      t1[IM]=-t1[IM];
      cDiv(cm,t1,CoupleConstant[i]);        /* CC[i]=cm/(1-(cm/V)*t1) */
    }
  }
  if (asym) {
    if (!orient_avg) {
      printz("CoupleConstant:{%.10g%+.10gi, %.10g%+.10gi, %.10g%+.10gi}\n",
        CoupleConstant[0][RE],CoupleConstant[0][IM],CoupleConstant[1][RE],
        CoupleConstant[1][IM],CoupleConstant[2][RE],CoupleConstant[2][IM]);
      fprintz(logfile, "CoupleConstant:{%.10g%+.10gi, %.10g%+.10gi, %.10g%+.10gi}\n",
        CoupleConstant[0][RE],CoupleConstant[0][IM],CoupleConstant[1][RE],
        CoupleConstant[1][IM],CoupleConstant[2][RE],CoupleConstant[2][IM]);
    }
  }
  else {
    memcpy(CoupleConstant[1],CoupleConstant[0],sizeof(doublecomplex));
    memcpy(CoupleConstant[2],CoupleConstant[0],sizeof(doublecomplex));
    if (!orient_avg) {
      printz("CoupleConstant:%.10g%+.10gi\n",CoupleConstant[0][RE],CoupleConstant[0][IM]);
      fprintz(logfile,"CoupleConstant:%.10g%+.10gi\n",CoupleConstant[0][RE],CoupleConstant[0][IM]);
    }
  }
  memcpy(res,CoupleConstant,3*sizeof(doublecomplex));
}

/*============================================================*/

void InitCC(char which)
    /* calculate cc and cc_sqrt */
{
  int i,j;

  for(i=0;i<Nmat;i++) {
    CoupleConstant(ref_index[i],which,cc[i]);
    for(j=0;j<3;j++) {
      cSqrt(cc[i][j],cc_sqrt[i][j]);
    }
  }
}

/*============================================================*/

double *ReadTableFile(char *sh_fname, int size_multiplier)
{
  FILE *ftab;
  double *tab_n;
  int size;
  char fname[200];
  int i;

  size=TAB_SIZE*size_multiplier;
  memory+=size*sizeof(double);
  if (!prognose) {
    /* allocate memory for tab_n */
    if ((tab_n = (double *) malloc(size*sizeof(double))) == NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc tab_n");
    /* open file */
    strcpy(fname,TAB_PATH);
    strcat(fname,sh_fname);
    if ((ftab=fopen(fname,"r"))==NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not open file '%s'",fname);
    /* scan file */
    for (i=0; i<size; i++) if (fscanf(ftab,"%lf\t",&(tab_n[i]))!=1)
      LogError(EC_ERROR,ALL,POSIT,"Scan error in file '%s'. Probably file is too small",fname);
    if (!feof(ftab)) LogError(EC_WARN,ONE,POSIT,"File '%s' is longer than specified size (%d)",fname,size);
    /* close file */
    fclose(ftab);
  }
  return tab_n;
}

/*============================================================*/

void ReadTables(void)
{
  int i, j, ymax, Rm2, Rm2x;

  tab1=ReadTableFile("/t1f.dat",1);
  tab2=ReadTableFile("/t2f.dat",6);
  tab3=ReadTableFile("/t3f.dat",3);
  tab4=ReadTableFile("/t4f.dat",18);
  tab5=ReadTableFile("/t5f.dat",6);
  tab6=ReadTableFile("/t6f.dat",36);
  tab7=ReadTableFile("/t7f.dat",1);
  tab8=ReadTableFile("/t8f.dat",6);
  tab9=ReadTableFile("/t9f.dat",1);
  tab10=ReadTableFile("/t10f.dat",6);

  if (!prognose) {
    /* allocate memory for tab_index */
    if ((tab_index = iMatrix(1,TAB_RMAX,0,TAB_RMAX)) == NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc tab_index");
    /* fill tab_index */
    Rm2=TAB_RMAX*TAB_RMAX;
    tab_index[1][0] = 0;
    for (i=1; i<=TAB_RMAX; i++) {
      Rm2x=Rm2-i*i;
      ymax = MIN(i,(int)floor(sqrt(Rm2x)));
      for (j=0; j<ymax; j++) {
        tab_index[i][j+1] = tab_index[i][j] + MIN(j,(int)floor(sqrt(Rm2x-j*j)))+1;
      }
      if (i<TAB_RMAX) tab_index[i+1][0] = tab_index[i][ymax] + MIN(ymax,(int)floor(sqrt(Rm2x-ymax*ymax)))+1;
    }
  }
  /* printz("P[5,3]=%d (should be 41)\n",tab_index[5][3]); */
}

/*============================================================*/

void FreeTables(void)
{
  Free_iMatrix(tab_index,1,TAB_RMAX,0);
  free(tab1);
  free(tab2);
  free(tab3);
  free(tab4);
  free(tab5);
  free(tab6);
  free(tab7);
  free(tab8);
  free(tab9); 
  free(tab10);
}

/*============================================================*/

void SaveMueller(double *muel)
   /* saves mueller matrix (averaged) to file */
{
  FILE *mueller;
  char stringbuffer[200];
  int i,j;
  double theta;
  clock_t tstart;

  tstart=clock();

  strcpy(stringbuffer,directory);
  strcat(stringbuffer,"/mueller");

  if ((mueller = fopen (stringbuffer, "w"))==NULL)
    LogError(EC_ERROR,ONE,POSIT,"File 'mueller' write error");

  fprintf(mueller,"theta s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44\n");

  for (i=0;i<nTheta;i++) {
    theta=i*dtheta_deg;
    fprintf(mueller,"%.2f",theta);
    for (j=0;j<16;j++) fprintf(mueller," %.10E",muel[16*i+j]);
    fprintf(mueller,"\n");
  }
  fclose(mueller);

  Timing_FileIO += clock() - tstart;
}

/*==============================================================*/

void SaveCS(double Cext, double Cabs)
   /* save calculated crossections (averaged) to file */
{        
  FILE *CCfile;
  char stringbuffer[200];
  clock_t tstart;

  tstart=clock();

  strcpy(stringbuffer,directory);
  strcat(stringbuffer,"/CrossSec");

  if ((CCfile=fopen(stringbuffer,"w"))==NULL)
    LogError(EC_ERROR,ONE,POSIT,"Could not write to file '%s'",stringbuffer);

  printf("Cext\t= %.10lg\nQext\t= %.10lg\n",Cext,Cext*inv_G);
  fprintf(CCfile,"Cext\t= %.10lg\nQext\t= %.10lg\n",Cext,Cext*inv_G);
  printf("Cabs\t= %.10lg\nQabs\t= %.10lg\n",Cabs,Cabs*inv_G);
  fprintf(CCfile,"Cabs\t= %.10lg\nQabs\t= %.10lg\n",Cabs,Cabs*inv_G);
  
  fclose(CCfile);

  Timing_FileIO += clock() - tstart;
}

/*============================================================*/

void calculate_one_orientation(double *res /* where to put results when doing orientation averaging */)
     /* performs calculation for one orientation */
{     
  clock_t tstart;
  
  if (orient_avg) {
    alph_deg=0;
    InitRotation();

    printz("\nORIENTATION STEP beta=%g gamma=%g\n",bet_deg,gam_deg);
    fprintz(logfile,"\nORIENTATION STEP beta=%g gamma=%g\n",bet_deg,gam_deg);
  }
  
  /* calculate scattered field for y - polarized incident light */
  printz("\nhere we go, calc Y\n\n");
  if (!orient_avg) fprintz(logfile,"\nhere we go, calc Y\n\n");
  InitCC('Y');
  if (symR && !all_dir && !scat_grid) {
    if (CalculateE('Y',CE_PARPER)==CHP_EXIT) return;
  }
  else { /* no rotational symmetry */
  /* in case of all_dir we run twice to get the full electric field */
  /* with incoming light in X and Y direction. In case of rotational */
  /* symmetry this is not needed but requires lots more programming */
  /* so we leave this optimization to a later time. */

    if(CalculateE('Y',CE_NORMAL)==CHP_EXIT) return;

    printz("\nhere we go, calc X\n\n");
    if (!orient_avg) fprintz(logfile,"\nhere we go, calc X\n\n");
    if(PolRelation==POL_LDR && !avg_inc_pol) InitCC('X');

    if(CalculateE('X',CE_NORMAL)==CHP_EXIT) return;
  }
  D("CalculateE complete");
  MuellerMatrix();
  D("MuellerMatrix complete");
  if (ringid==ROOT && orient_avg) {
    tstart=clock();
    printf("\nError of alpha integration (Mueller) is %g\n",
      Romberg1D(parms_alpha,block_theta,muel_alpha,&res[2]));
    memcpy(res,&muel_alpha[-2],2*sizeof(double));
    D("integration alpha complete");
    Timing_Integration += clock() - tstart;
  }
  TotalEval++;
}

/*============================================================*/

void orient_integrand(int beta_i, int gamma_i, double *res) 
   /* function that provides interface with Romberg integration */
{
  BcastOrient(&beta_i,&gamma_i,&finish_avg);
  if (finish_avg) return;
  
  bet_deg=beta_int.val[beta_i];
  gam_deg=gamma_int.val[gamma_i];
  calculate_one_orientation(res);
}

/*============================================================*/

void AllocateEverything(void)
  /* allocates a lot of arrays and performs memory analysis */
{
  int temp_int;
  double memtot,memmax;

  /* allocate all the memory */
  if (!prognose) if ((xvec = cVector(nlocalRows)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc xvec");
  memory+=nlocalRows*sizeof(doublecomplex);
  if (!prognose) if ((rvec = cVector(nlocalRows))== NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc rvec");
  memory+=nlocalRows*sizeof(doublecomplex);
  if (!prognose) if ((pvec = cVector(nlocalRows)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc pvec");
  memory+=nlocalRows*sizeof(doublecomplex);
  if (!prognose) if ((Einc = cVector(nlocalRows)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Einc");
  memory+=nlocalRows*sizeof(doublecomplex);
  if (!prognose) if ((Avecbuffer = cVector(nlocalRows)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Avecbuffer");
  memory+=nlocalRows*sizeof(doublecomplex);
  if (IterMethod==IT_BICGSTAB || IterMethod==IT_QMR_CS) {
    /* additional vectors for iterative methods */
    if (!prognose) if ((vec1 = cVector(nlocalRows)) == NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc vec1");
    memory+=nlocalRows*sizeof(doublecomplex);
    if (!prognose) if ((vec2 = cVector(nlocalRows)) == NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc vec2");
    memory+=nlocalRows*sizeof(doublecomplex);
    if (!prognose) if ((vec3 = cVector(nlocalRows)) == NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc vec3");
    memory+=nlocalRows*sizeof(doublecomplex);
  }
  if (yzplane) {
    if (!prognose) {
      if ((EplaneX = cVector(2*nTheta)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc EplaneX");
      if ((EplaneY = cVector(2*nTheta)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc EplaneY");
#ifdef PARALLEL
      /* buffer for accumulate operation */
      if ((Eplane_buffer = dVector(0,4*nTheta-1)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc Eplane_buffer");
#endif
    }
    memory+=4*nTheta*sizeof(doublecomplex);
#ifdef PARALLEL
    memory+=4*nTheta*sizeof(double);  /* extra memory for buffer */
#endif
  }
  if (all_dir) {
    ReadAlldirParms(alldir_parms);
    /* calculate size of vectors; 4 - because first it is used to store
       per and par components of the field, and only afterwards squares */
    temp_int=4*theta_int.N*phi_int.N;
    if (!prognose) {
      if ((E2_alldir = dVector(0,temp_int-1)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc E2_alldir");
#ifdef PARALLEL
      /* buffer for accumulate operation */
      if ((E2_alldir_buffer = dVector(0,temp_int-1)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc E2_alldir_buffer");
#endif
    }
    memory+=temp_int*sizeof(double);
#ifdef PARALLEL
    memory+=temp_int*sizeof(double);  /* extra memory for buffer */
#endif
  }
  if (scat_grid) {
    ReadScatGridParms(scat_grid_parms);
    /* calculate size of vectors - holds all per-par combinations*/
    temp_int=2*angles.N;
    if (!prognose) {
      if ((EgridX = cVector(temp_int)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc EgridX");
      if ((EgridY = cVector(temp_int)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc EgridY");
      if (phi_integr && ringid==ROOT)
        if ((muel_phi = dVector(0,16*angles.phi.N+1)) == NULL)
          LogError(EC_ERROR,ONE,POSIT,"Could not malloc muel_phi");
        if ((muel_phi1 = dVector(0,16*angles.phi.N+1)) == NULL)
          LogError(EC_ERROR,ONE,POSIT,"Could not malloc muel_phi1");
#ifdef PARALLEL
      /* buffer for accumulate operation */
      if ((Egrid_buffer = dVector(0,2*temp_int-1)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc Egrid_buffer");
#endif
    }
    memory+=2*temp_int*sizeof(doublecomplex);
    if (phi_integr && ringid==ROOT) memory+=2*16*angles.phi.N*sizeof(double);
#ifdef PARALLEL
    memory+=2*temp_int*sizeof(double);  /* extra memory for buffer */
#endif
  }
  if (orient_avg) {
    if (!prognose && ringid==ROOT) {
      if ((ampl_alphaX = cVector(2*nTheta*alpha_int.N)) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc ampl_alphaX");
      if ((ampl_alphaY = cVector(2*nTheta*alpha_int.N)) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc ampl_alphaY");
      if ((muel_alpha = dVector(0,block_theta*alpha_int.N+1)) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc muel_alpha");
      muel_alpha+=2;
      if ((out = dVector(0,block_theta+1)) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc out");
    }
    memory += 4*nTheta*alpha_int.N*sizeof(doublecomplex) +
              (block_theta*(alpha_int.N+1)+4)*sizeof(double);
  }
  /* estimate of the memory (only the fastest scaling part):
       MatVec - 288*Ndip (more exactly gridX*gridY*gridZ*72)
       others - nvoid_Ndip*271(+144 for BiCGStab and QMR_CS)
     PARALLEL: above is total; division over processors of MatVec is uniform,
               others - according to local_nvoid_Ndip  */
  memtot=memory/MBYTE;
  AccumulateMax(&memtot,&memmax);
  printz("Total memory usage: %.1f Mb\n",memtot);
  fprintz(logfile,"Total memory usage: %.1f Mb\n",memtot);
#ifdef PARALLEL
  fprintz(logfile,"Maximum memory usage of single processor: %.1f Mb\n",memmax);
#endif
  if (prognose) Stop(0);
}

/*============================================================*/

void FreeEverything(void)
  /* frees all allocated vectors; should not be called in prognose mode,
     since arrays are not actually allocated */
{
  if (IntRelation == G_SO) FreeTables();
  Free_FFT_Dmat();
  Free_cVector(xvec);
  Free_cVector(rvec);
  Free_cVector(pvec);
  Free_cVector(Einc);
  Free_cVector(Avecbuffer);
  if (IterMethod==IT_BICGSTAB || IterMethod==IT_QMR_CS) {
    Free_cVector(vec1);
    Free_cVector(vec2);
    Free_cVector(vec3);
  }
  if (yzplane) {
    Free_cVector(EplaneX);
    Free_cVector(EplaneY);
#ifdef PARALLEL
    Free_dVector(Eplane_buffer,0);
#endif
  }
  if (all_dir) {
    free(theta_int.val);
    free(phi_int.val);
    Free_dVector(E2_alldir,0);
#ifdef PARALLEL
    Free_dVector(E2_alldir_buffer,0);
#endif
  }
  if (scat_grid) {
    free(angles.theta.val);
    free(angles.phi.val);
    Free_cVector(EgridX);
    Free_cVector(EgridY);
    if (phi_integr && ringid==ROOT) {
      free(muel_phi);
      free(muel_phi1);
    }
#ifdef PARALLEL
    Free_dVector(Egrid_buffer,0);
#endif
  }
  /* these 3 were allocated in MakeParticle */
  Free_dVector(DipoleCoord,0);
  free(position);
  free(material);

  if (orient_avg) {
    if (ringid==ROOT) {
      Free_cVector(ampl_alphaX);
      Free_cVector(ampl_alphaY);
      Free_dVector(muel_alpha-2,0);
      Free_dVector(out,0);
    }
    free(alpha_int.val);
    free(beta_int.val);
    free(gamma_int.val);
  }
}

/*============================================================*/

void Calculator (void)
{
  char bufstr[200];

  /* initialize variables */
  dtheta_deg = 180.0 / ((double)(nTheta-1));
  dtheta_rad = Deg2Rad(dtheta_deg);
  block_theta=16*nTheta;
  /* if not enough symmetry, calculate for +- theta (for one plane) */
  if (!(symY || orient_avg)) nTheta=2*(nTheta-1);
  finish_avg=FALSE;
  /* read tables if needed */
  if (IntRelation == G_SO) ReadTables();
  /* initialize D matrix (for matvec) */
  D("InitDmatrix started");
  InitDmatrix();
  D("InitDmatrix complete");
  /* allocate most (that is not already allocated; perform memory analysis */
  AllocateEverything();
  /* finish init */
  if (!orient_avg) alpha_int.N=1;
  Timing_Init = clock() - tstart_main;
  /* main calculation part */
  if (orient_avg) {
    if (ringid==ROOT) {
      sprintf(bufstr,"%s/log_orient_avg",directory);
      Romberg2D(parms,orient_integrand,block_theta+2,out,bufstr);
      finish_avg=TRUE;
      BcastOrient(&finish_avg,&finish_avg,&finish_avg);
      SaveMueller(&out[2]);
      SaveCS(out[0],out[1]);
    }
    else while (!finish_avg) orient_integrand(0,0,NULL);
  }
  else calculate_one_orientation(NULL);
  /* cleaning */
  FreeEverything();
}



