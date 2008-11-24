/* FILE: calculator.c
 * AUTH: Maxim Yurkin
 * DESCR: All the initialization is done here before actually
 *        calculating internal fields. Calculation of couple constants.
 *
 *        Previous versions were by Alfons Hoekstra
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
#include "timing.h"

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in crosssec.c */
extern const Parms_1D parms[2],parms_alpha;
extern const angle_set beta_int,gamma_int,theta_int,phi_int;
/* defined and initialized in param.c */
extern const int avg_inc_pol;
extern const char alldir_parms[],scat_grid_parms[];
/* defined and initialized in timing.c */
extern TIME_TYPE Timing_Init;
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
doublecomplex *expsX,*expsY,*expsZ;   /* arrays of exponents along 3 axes (for calc_field) */
/* used in iterative.c */
doublecomplex *rvec,*vec1,*vec2,*vec3,*Avecbuffer;  /* vectors for iterative solvers */

/* LOCAL VARIABLES */

static size_t block_theta;    /* size of one block of mueller matrix - 16*nTheta */
static int finish_avg;        /* whether to stop orientation averaging */
static double *out;           /* used to collect both mueller matrix and integral
                                 scattering quantities when orient_avg */
/* EXTERNAL FUNCTIONS */

/* CalculateE.c */
extern int CalculateE(char which,int type);
extern void MuellerMatrix(void);

/*============================================================*/

static void CoupleConstant(doublecomplex *mrel,const char which,doublecomplex *res)
{
  doublecomplex coup_con[3];
  doublecomplex tempa,tempb,cm,m2,t1;
  double temp,V,b1,b2,b3;
  int i,imax,j,jmax;  /* counters: i is for 'asym', j is for 'anysotropy' */
  double S,prop2[3];
  int asym;           /* whether polarizability is asymmetric (for isotropic m) */
  const double *incPol;
  int pol_avg=TRUE;

  asym = (PolRelation==POL_CLDR || PolRelation==POL_SO);
  /* this should never happen !!! */
  if (asym && anisotropy) LogError(EC_ERROR,ONE_POS,"Incompatibility error in CoupleConstant");
  if (asym) imax=3;
  else imax=1;
  if (anisotropy) jmax=3;
  else jmax=1;
  if (PolRelation==POL_LDR || PolRelation==POL_CLDR) {
    b1=LDR_B1;
    b2=LDR_B2;
    b3=LDR_B3;
  }
  else if (PolRelation==POL_SO) {
    b1=SO_B1;
    b2=SO_B2;
    b3=SO_B3;
  }
  /* calculate the CM couple constant CC=(3V/4pi)*(m^2-1)/(m^2+2) */
  V=gridspace*gridspace*gridspace;     /* volume of one dipole */
  temp = (3*V)/(4*PI);
  for (j=0;j<jmax;j++) {
    cSquare(mrel[j],m2);                    /* m2=m^2 */
    tempa[RE] = m2[RE] - 1.0;
    tempa[IM] = tempb[IM] = m2[IM];
    tempb[RE] = m2[RE] + 2.0;
    cDiv(tempa,tempb,coup_con[j]);
    coup_con[j][RE] *= temp;
    coup_con[j][IM] *= temp;

    if (PolRelation!=POL_CM) {
      if (PolRelation==POL_LDR || PolRelation==POL_CLDR || PolRelation==POL_SO) {
        /* set prop_i^2 */
        for (i=0;i<3;i++) {
          if (pol_avg && PolRelation==POL_SO) prop2[i]=ONE_THIRD;
          else prop2[i]=prop[i]*prop[i];
        }
        /* determine S coefficient for LDR */
        if (PolRelation==POL_LDR) {
          if (avg_inc_pol) S=0.5*(1-DotProd(prop2,prop2));
          else {
            if (which=='X') incPol=incPolX;
            else if (which=='Y') incPol=incPolY;
            S=prop2[0]*incPol[0]*incPol[0]+prop2[1]*incPol[1]*incPol[1]
             +prop2[2]*incPol[2]*incPol[2];
          }
        }
      }
      cEqual(coup_con[j],cm);
      for (i=0;i<imax;i++) {
        /* RR correction */
        t1[RE]=0.0;
        t1[IM]=2*kd*kd*kd/3;                       /* t1=2/3*i*kd^3         */
        /* plus more advanced corrections */
        if (PolRelation==POL_FCD)  /* t1+=((4/3)kd^2+(2/3pi)*log((pi-kd)/(pi+kd))*kd^3) */
          t1[RE]+=2*ONE_THIRD*kd*kd*(2+kd*INV_PI*log((PI-kd)/(PI+kd)));
        else if (PolRelation==POL_LDR || PolRelation==POL_CLDR || PolRelation==POL_SO) {
          if (PolRelation!=POL_LDR) S=prop2[i];
          t1[RE]+=(b1+(b2+b3*S)*m2[RE])*kd*kd;   /* t1+=(b1+(b2+b3*S)*m^2)*kd^2  */
          t1[IM]+=(b2+b3*S)*m2[IM]*kd*kd;
        }
        /* CC[i]=cm/(1-(cm/V)*t1) */
        cMultReal(1.0/V,t1,t1);
        cMultSelf(t1,cm);
        t1[RE]=1-t1[RE];
        t1[IM]=-t1[IM];
        /* 'i+j' is not robust. It assumes that only one counter is used */
        cDiv(cm,t1,coup_con[i+j]);
      }
    }
  }
  if (asym || anisotropy) {
    if (!orient_avg) {
      PRINTBOTHZ(logfile, "CoupleConstant:(%.10g%+.10gi,%.10g%+.10gi,%.10g%+.10gi)\n",
                 coup_con[0][RE],coup_con[0][IM],coup_con[1][RE],
                 coup_con[1][IM],coup_con[2][RE],coup_con[2][IM]);
    }
  }
  else {
    cEqual(coup_con[0],coup_con[1]);
    cEqual(coup_con[0],coup_con[2]);
    if (!orient_avg) {
      PRINTBOTHZ(logfile,"CoupleConstant:%.10g%+.10gi\n",
                 coup_con[0][RE],coup_con[0][IM]);
    }
  }
  memcpy(res,coup_con,3*sizeof(doublecomplex));
}

/*============================================================*/

static void InitCC(const char which)
    /* calculate cc and cc_sqrt */
{
  int i,j;

  for(i=0;i<Nmat;i++) {
    CoupleConstant(ref_index+Ncomp*i,which,cc[i]);
    for(j=0;j<3;j++) cSqrt(cc[i][j],cc_sqrt[i][j]);
  }
}

/*============================================================*/

static double *ReadTableFile(const char *sh_fname,const int size_multiplier)
{
  FILE *ftab;
  double *tab_n;
  int size;
  char fname[MAX_FNAME];
  int i;

  size=TAB_SIZE*size_multiplier;
  memory+=size*sizeof(double);
  if (!prognose) {
    /* allocate memory for tab_n */
    MALLOC_VECTOR(tab_n,double,size,ALL);
    /* open file */
    strcpy(fname,TAB_PATH);
    strcat(fname,sh_fname);
    ftab=FOpenErr(fname,"r",ALL_POS);
    /* scan file */
    for (i=0; i<size; i++) if (fscanf(ftab,"%lf\t",&(tab_n[i]))!=1)
      LogError(EC_ERROR,ALL_POS,"Scan error in file '%s'. Probably file is too small",fname);
    if (!feof(ftab))
      LogError(EC_WARN,ONE_POS,"File '%s' is longer than specified size (%d)",fname,size);
    /* close file */
    FCloseErr(ftab,fname,ALL_POS);
  }
  return tab_n;
}

/*============================================================*/

static void ReadTables(void)
{
  int i, j, ymax, Rm2, Rm2x;

  tab1=ReadTableFile(TAB_FNAME(1),1);
  tab2=ReadTableFile(TAB_FNAME(2),6);
  tab3=ReadTableFile(TAB_FNAME(3),3);
  tab4=ReadTableFile(TAB_FNAME(4),18);
  tab5=ReadTableFile(TAB_FNAME(5),6);
  tab6=ReadTableFile(TAB_FNAME(6),36);
  tab7=ReadTableFile(TAB_FNAME(7),1);
  tab8=ReadTableFile(TAB_FNAME(8),6);
  tab9=ReadTableFile(TAB_FNAME(9),1);
  tab10=ReadTableFile(TAB_FNAME(10),6);

  if (!prognose) {
    /* allocate memory for tab_index */
    MALLOC_IMATRIX(tab_index,1,TAB_RMAX,0,TAB_RMAX,ALL);
    /* fill tab_index */
    Rm2=TAB_RMAX*TAB_RMAX;
    tab_index[1][0] = 0;
    for (i=1; i<=TAB_RMAX; i++) {
      Rm2x=Rm2-i*i;
      ymax = MIN(i,(int)floor(sqrt(Rm2x)));
      for (j=0; j<ymax; j++) {
        tab_index[i][j+1] = tab_index[i][j] + MIN(j,(int)floor(sqrt(Rm2x-j*j)))+1;
      }
      if (i<TAB_RMAX) tab_index[i+1][0] = tab_index[i][ymax]
      + MIN(ymax,(int)floor(sqrt(Rm2x-ymax*ymax)))+1;
    }
  }
  /* PRINTZ("P[5,3]=%d (should be 41)\n",tab_index[5][3]); */
}

/*============================================================*/

static void FreeTables(void)
{
  Free_iMatrix(tab_index,1,TAB_RMAX,0);
  Free_general(tab1);
  Free_general(tab2);
  Free_general(tab3);
  Free_general(tab4);
  Free_general(tab5);
  Free_general(tab6);
  Free_general(tab7);
  Free_general(tab8);
  Free_general(tab9);
  Free_general(tab10);
}

/*============================================================*/

static void SaveMueller(double *muel)
   /* saves mueller matrix (averaged) to file */
{
  FILE *mueller;
  char fname[MAX_FNAME];
  int i,j;
  double theta;
  TIME_TYPE tstart;

  tstart=GET_TIME();

  strcpy(fname,directory);
  strcat(fname,"/" F_MUEL);

  mueller=FOpenErr(fname,"w",ONE_POS);

  fprintf(mueller,"theta s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44\n");

  for (i=0;i<nTheta;i++) {
    theta=i*dtheta_deg;
    fprintf(mueller,"%.2f",theta);
    for (j=0;j<16;j++) fprintf(mueller," %.10E",muel[16*i+j]);
    fprintf(mueller,"\n");
  }
  FCloseErr(mueller,F_MUEL,ONE_POS);

  Timing_FileIO += GET_TIME() - tstart;
}

/*==============================================================*/

static void SaveCS(const double Cext,const double Cabs)
   /* save calculated crossections (averaged) to file */
{
  FILE *CCfile;
  char fname[MAX_FNAME];
  TIME_TYPE tstart;

  tstart=GET_TIME();

  strcpy(fname,directory);
  strcat(fname,"/" F_CS);

  CCfile=FOpenErr(fname,"w",ONE_POS);

  PrintBoth(CCfile,"Cext\t= %.10g\nQext\t= %.10g\n",Cext,Cext*inv_G);
  PrintBoth(CCfile,"Cabs\t= %.10g\nQabs\t= %.10g\n",Cabs,Cabs*inv_G);

  FCloseErr(CCfile,F_CS,ONE_POS);

  Timing_FileIO += GET_TIME() - tstart;
}

/*============================================================*/

static void calculate_one_orientation(double *res)
  /* performs calculation for one orientation; may do orientation averaging and put
     the result in res */
{
  TIME_TYPE tstart;

  if (orient_avg) {
    alph_deg=0;
    InitRotation();
    PRINTBOTHZ(logfile,"\nORIENTATION STEP beta=%g gamma=%g\n",bet_deg,gam_deg);
  }

  /* calculate scattered field for y - polarized incident light */
  PRINTZ("\nhere we go, calc Y\n\n");
  if (!orient_avg) FPRINTZ(logfile,"\nhere we go, calc Y\n\n");
  InitCC('Y');
  if (symR && !scat_grid) {
    if (CalculateE('Y',CE_PARPER)==CHP_EXIT) return;
  }
  else { /* no rotational symmetry */
  /* in case of scat_grid we run twice to get the full electric field */
  /* with incoming light in X and Y direction. In case of rotational */
  /* symmetry this is not needed but requires lots more programming */
  /* so we leave this optimization to a later time. */
    if(CalculateE('Y',CE_NORMAL)==CHP_EXIT) return;

    PRINTZ("\nhere we go, calc X\n\n");
    if (!orient_avg) FPRINTZ(logfile,"\nhere we go, calc X\n\n");
    if(PolRelation==POL_LDR && !avg_inc_pol) InitCC('X');

    if(CalculateE('X',CE_NORMAL)==CHP_EXIT) return;
  }
  D("CalculateE finished");
  MuellerMatrix();
  D("MuellerMatrix finished");
  if (ringid==ROOT && orient_avg) {
    tstart=GET_TIME();
    printf("\nError of alpha integration (Mueller) is %g\n",
      Romberg1D(parms_alpha,block_theta,muel_alpha,res+2));
    memcpy(res,muel_alpha-2,2*sizeof(double));
    D("Integration over alpha completed on ROOT");
    Timing_Integration += GET_TIME() - tstart;
  }
  TotalEval++;
}

/*============================================================*/

static double orient_integrand(int beta_i,int gamma_i, double *res)
   /* function that provides interface with Romberg integration */
{
  BcastOrient(&beta_i,&gamma_i,&finish_avg);
  if (finish_avg) return 0;

  bet_deg=beta_int.val[beta_i];
  gam_deg=gamma_int.val[gamma_i];
  calculate_one_orientation(res);
  return 0;
}

/*============================================================*/

static void AllocateEverything(void)
  /* allocates a lot of arrays and performs memory analysis */
{
  double tmp;
  size_t temp_int;
  double memmax;

  /* allocate all the memory */
  tmp=sizeof(doublecomplex)*(double)nlocalRows;
  if (!prognose) {
    MALLOC_VECTOR(xvec,complex,nlocalRows,ALL);
    MALLOC_VECTOR(rvec,complex,nlocalRows,ALL);
    MALLOC_VECTOR(pvec,complex,nlocalRows,ALL);
    MALLOC_VECTOR(Einc,complex,nlocalRows,ALL);
    MALLOC_VECTOR(Avecbuffer,complex,nlocalRows,ALL);
  }
  memory+=5*tmp;
  if (IterMethod==IT_BICGSTAB || IterMethod==IT_QMR_CS) {
    /* additional vectors for iterative methods */
    if (!prognose) {
      MALLOC_VECTOR(vec1,complex,nlocalRows,ALL);
      MALLOC_VECTOR(vec2,complex,nlocalRows,ALL);
      MALLOC_VECTOR(vec3,complex,nlocalRows,ALL);
    }
    memory+=3*tmp;
  }
  MALLOC_VECTOR(expsX,complex,boxX,ALL);
  MALLOC_VECTOR(expsY,complex,boxY,ALL);
  MALLOC_VECTOR(expsZ,complex,local_Nz_unif,ALL);
  if (yzplane) {
    tmp=2*(double)nTheta;
    if (!prognose) {
      CheckOverflow(2*tmp,ONE_POS,"AllocateEverything()");
      temp_int=tmp;
      MALLOC_VECTOR(EplaneX,complex,temp_int,ALL);
      MALLOC_VECTOR(EplaneY,complex,temp_int,ALL);
    }
    memory+=2*tmp*sizeof(doublecomplex);
#ifdef PARALLEL
    if (ringid==ROOT) { /* buffer for accumulate operation */
      if (!prognose) MALLOC_VECTOR(Eplane_buffer,double,2*temp_int,ONE);
      memory+=2*tmp*sizeof(double);
    }
#endif
  }
  if (all_dir) {
    ReadAlldirParms(alldir_parms);
    /* calculate size of vectors; 4 - because first it is used to store
       per and par components of the field, and only afterwards squares */
    tmp=4*((double)theta_int.N)*phi_int.N;
    if (!prognose) {
      CheckOverflow(tmp,ONE_POS,"AllocateEverything()");
      temp_int=tmp;
      MALLOC_VECTOR(E2_alldir,double,temp_int,ALL);
    }
    memory+=tmp*sizeof(double);
#ifdef PARALLEL
    if (ringid==ROOT) { /* buffer for accumulate operation */
      if (!prognose) MALLOC_VECTOR(E2_alldir_buffer,double,temp_int,ONE);
      memory+=tmp*sizeof(double);
    }
#endif
  }
  if (scat_grid) {
    ReadScatGridParms(scat_grid_parms);
    /* calculate size of vectors - holds all per-par combinations*/
    tmp=2*(double)angles.N;
    if (!prognose) {
      CheckOverflow(2*tmp,ONE_POS,"AllocateEverything()");
      temp_int=tmp;
      MALLOC_VECTOR(EgridX,complex,temp_int,ALL);
      MALLOC_VECTOR(EgridY,complex,temp_int,ALL);
    }
    memory+=2*tmp*sizeof(doublecomplex);
#ifdef PARALLEL
    if (ringid==ROOT) {  /* buffer for accumulate operation */
      if (!prognose) MALLOC_VECTOR(Egrid_buffer,double,2*temp_int,ONE);
      memory+=2*tmp*sizeof(double);
    }
#endif
    if (phi_integr && ringid==ROOT) {
      tmp=16*(double)angles.phi.N;
      if (!prognose) {
        CheckOverflow(tmp,ONE_POS,"AllocateEverything()");
        temp_int=tmp;
        MALLOC_VECTOR(muel_phi,double,temp_int,ONE);
        MALLOC_VECTOR(muel_phi1,double,temp_int,ONE);
      }
      memory+=2*tmp*sizeof(double);
    }
  }
  if (orient_avg) {
    tmp=2*((double)nTheta)*alpha_int.N;
    if (!prognose) {
      /* this covers these 2 and next 2 mallocs */
      CheckOverflow(8*tmp+2,ONE_POS,"AllocateEverything()");
      temp_int=tmp;
      MALLOC_VECTOR(ampl_alphaX,complex,temp_int,ONE);
      MALLOC_VECTOR(ampl_alphaY,complex,temp_int,ONE);
    }
    memory += 2*tmp*sizeof(doublecomplex);
    if (ringid==ROOT) {
      if (!prognose) {
        MALLOC_VECTOR(muel_alpha,double,block_theta*alpha_int.N+2,ONE);
        muel_alpha+=2;
        MALLOC_VECTOR(out,double,block_theta+2,ONE);
      }
      memory += (8*tmp*(1+1.0/alpha_int.N)+4)*sizeof(double);
    }
  }
  /* estimate of the memory (only the fastest scaling part):
       MatVec - (288+384nprocs/boxX [+192/nprocs])*Ndip
                more exactly: gridX*gridY*gridZ*(36+48nprocs/boxX [+24/nprocs])
                value in [] is only for parallel mode
       others - nvoid_Ndip*271(+144 for BiCGStab and QMR_CS)
     PARALLEL: above is total; division over processors of MatVec is uniform,
               others - according to local_nvoid_Ndip  */
  memory/=MBYTE;
  AccumulateMax(&memory,&memmax);
  PRINTBOTHZ(logfile,"Total memory usage: %.1f Mb\n",memory);
#ifdef PARALLEL
  PRINTBOTHZ(logfile,"Maximum memory usage of single processor: %.1f Mb\n",memmax);
#endif
}

/*============================================================*/

static void FreeEverything(void)
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
  Free_cVector(expsX);
  Free_cVector(expsY);
  Free_cVector(expsZ);
  if (yzplane) {
    Free_cVector(EplaneX);
    Free_cVector(EplaneY);
#ifdef PARALLEL
    Free_general(Eplane_buffer);
#endif
  }
  if (all_dir) {
    Free_general(theta_int.val);
    Free_general(phi_int.val);
    Free_general(E2_alldir);
#ifdef PARALLEL
    Free_general(E2_alldir_buffer);
#endif
  }
  if (scat_grid) {
    Free_general(angles.theta.val);
    Free_general(angles.phi.val);
    Free_cVector(EgridX);
    Free_cVector(EgridY);
    if (phi_integr && ringid==ROOT) {
      Free_general(muel_phi);
      Free_general(muel_phi1);
    }
#ifdef PARALLEL
    Free_general(Egrid_buffer);
#endif
  }
  /* these 3 were allocated in MakeParticle */
  Free_general(DipoleCoord);
  Free_general(position);
  Free_general(material);

  if (orient_avg) {
    if (ringid==ROOT) {
      Free_cVector(ampl_alphaX);
      Free_cVector(ampl_alphaY);
      Free_general(muel_alpha-2);
      Free_general(out);
    }
    Free_general(alpha_int.val);
    Free_general(beta_int.val);
    Free_general(gamma_int.val);
  }
}

/*============================================================*/

void Calculator (void)
{
  char fname[MAX_FNAME];

  /* initialize variables */
  dtheta_deg = 180.0 / ((double)(nTheta-1));
  dtheta_rad = Deg2Rad(dtheta_deg);
  block_theta= 16*(size_t)nTheta;
  /* if not enough symmetry, calculate for +- theta (for one plane) */
  if (!(symY || orient_avg)) nTheta=2*(nTheta-1);
  finish_avg=FALSE;
  /* read tables if needed */
  if (IntRelation == G_SO) ReadTables();
  /* initialize D matrix (for matvec) */
  D("InitDmatrix started");
  InitDmatrix();
  D("InitDmatrix finished");
  /* allocate most (that is not already allocated; perform memory analysis */
  AllocateEverything();
  /* finish init */
  if (!orient_avg) alpha_int.N=1;
  Timing_Init = GET_TIME() - tstart_main;
  /* prognose stops here */
  if (prognose) return;
  /* main calculation part */
  if (orient_avg) {
    if (ringid==ROOT) {
      sprintf(fname,"%s/" F_LOG_ORAVG,directory);
      D("Romberg2D started on ROOT");
      Romberg2D(parms,orient_integrand,block_theta+2,out,fname);
      D("Romberg2D finished on ROOT");
      finish_avg=TRUE;
      /* first two are dummy variables; this call corresponds to similar call in orient_integrand
         by other processors */
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
