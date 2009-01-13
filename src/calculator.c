/* FILE: calculator.c
 * AUTH: Alfons Hoekstra
 * DESCR: All the initialization is done here before actually
 *        calculating internal fields
 *
 *        Currently is developed by Maxim Yurkin
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "cmplx.h"
#include "Romberg.h"
#include "const.h"
#include "comm.h"
#include "debug.h"
#include "memory.h"
#include "crosssec.h"

double dpl;
double lambda;
doublecomplex *x,*p,*r,*vec1,*vec2,*vec3,*Einc;
doublecomplex *Eper, *Epar; /* per and par components of scattered E (for scattering in one plane) */
doublecomplex *EgridX,*EgridY; /* E calculated on a grid for many different directions (holds Eper and Epar)
                                  for two incident fields */
double *E2_alldir;  /* square of E, calculated for alldir */
double *Eparper_buffer; /* buffer to accumulate Eper and Epar */
double *Egrid_buffer,*E2_alldir_buffer; /* buffers to accumulate Egrid and E2_alldir */

extern FILE *logfile;
extern int prognose;
extern doublecomplex ref_index[MAXNMAT];
extern Parms_1D parms[];
extern angle_set theta_int,phi_int;
extern int yzplane;
extern char directory[];
extern double prop[3];          /* propagation vector */
extern double incPolX[3],incPolY[3]; /* incident polarizations corresponding
				        to X and Y for default incidence */
extern int orient_avg;
extern angle_set alpha_int, beta_int, gamma_int;
extern double alph, bet, gam, alph_deg, bet_deg, gam_deg;

extern int all_dir,scat_grid,phi_integr;
extern scat_grid_angles angles;
extern char alldir_parms[],scat_grid_parms[];
double *muel_phi;   /* used to store values of Mueller matrix for different phi (to integrate) */
double *muel_phi1;  /* additional for integrating with different multipliers */

extern int shape;             /* shape of the particle */
extern int Nmat,mat_count[MAXNMAT+1];
extern double gridspace;

extern int symR;
extern int PolRelation;
extern int IterMethod;

int nDip;             /* number of dipoles */
int nTheta;	      /* number of angles in scattering profile */
int nlocalRows;       /* number of local rows of decomposition (only real dipoles) */
doublecomplex *Avecbuffer; /* buffer to hold result of local Matvec */
double *DipoleCoord;   /* vector to hold the coordinates of the dipoles */
short int *position;  /* position of the dipoles */
char *material;        /* material: index for cc */
double WaveNum;       /* wavenumber of incident light */
int memory=0;         /* total memory usage in bytes */
double *incPol;		/* used polarization */
doublecomplex cc[MAXNMAT][3];         /* couple constants */
doublecomplex cc_sqrt[MAXNMAT][3];   /* sqrt of couple constants */
double kd;             /* k*d=2*PI/dpl */
double *tab1,*tab2,*tab3,*tab4,*tab5,*tab6,*tab7,*tab8,*tab9,*tab10; /* tables of integrals */
int **tab_index;      /* matrix for indexation of table arrays */

extern int avg_inc_pol;

doublecomplex *ampl_alpha;  /* storing amplitude matrix for different values of alpha */
double *muel_alpha;        /* stroring mueller matrix for different values of alpha */
Parms_1D parms_alpha;
int block_theta;       /* size of one block of mueller matrix - 16*nTheta */

extern unsigned long TotalEval;
extern clock_t Timing_FileIO, Timing_Integration, Timing_Init;
extern clock_t tstart_main;

int finish_avg;

extern void init_rotation(void);
extern void init_Dmatrix(void);
extern void CalculateE(char which,int type);
extern void mueller_matrix(void);
extern void free_FFT_Dmat(void);

/*============================================================*/
                                         
void save_in_mie(void)
{
  FILE *out;
  double r,rin;
  char buffer[1024];
  int i,Nmie,inner,outer;
  extern int beamtype;
  extern double beam_w0,beam_x0,beam_y0,beam_z0,coat_ratio;
  clock_t tstart;

  inner=0; outer=0;
  if (shape==COATED) {
    inner=1;
    outer=0;
  }
  if (beamtype==PLANE) Nmie=1; else Nmie=2;
  r=(2.0*PI/dpl)*pow((0.75/PI)*nvoid_Ndip,0.333333333333);

  inner=0; outer=0; rin=r;
  if (shape==COATED) {
    inner=1;
    outer=0;
    rin=coat_ratio*r;
    rin=2.0*PI*pow((0.75/PI)*mat_count[1]/(dpl*dpl*dpl),0.333333333333);
  }

  tstart=clock();

  for(i=0;i<Nmie;i++) {
    sprintz(buffer,"%s/in_mie%i",directory,i);
    
    if ((out=fopen (buffer, "w"))==NULL)
      LogError(EC_ERROR,ONE,POSIT,"File 'in_mie' write error");
    
    if (beamtype==PLANE || (beam_x0==0.0 && beam_y0==0.0 && beam_z0==0.0)) {
      fprintz(out,"%i\n",nTheta);
      fprintz(out,"%.8g\n%.8g\n",rin,r);
      fprintz(out,"%.9g\n%.8g\n%.8g\n%.8g\n",
	      ref_index[inner][re],ref_index[inner][im],
	      ref_index[outer][re],ref_index[outer][im]);
      fprintz(out,"1.0\n%.7E\n",lambda);
      if (beamtype==PLANE) {
	fprintz(out,"%.7E\n",1e6*lambda);
	fprintz(out,"0\n");
      }
      else {
	fprintz(out,"%.7E\n",beam_w0);
	fprintz(out,"99\n");
      }
    } 
    else { /* focussed beam */
      fprintz(out,"%.7E\n",lambda);
      fprintz(out,"%.7E\n",beam_w0);
      fprintz(out,"%.7E\n",r/PI*lambda);
      fprintz(out,"%.7E\n",ref_index[0][re]);
      fprintz(out,"%.7E\n1.0\n",ref_index[0][im]);
      if (i==1) fprintz(out,"%.7E\n%.7E\n%.7E\n",beam_x0,beam_y0,beam_z0);
      else fprintz(out,"%.7E\n%.7E\n%.7E\n",beam_y0,-beam_x0,beam_z0);
      fprintz(out,"%.1f\n0.0\n",90.0*i);
      fprintz(out,"%.15E\n",180.0/(nTheta-1));
      fprintz(out,"%.7E\n",(double)(nTheta-1));
      fprintz(out,"12.0\n0\n");
    }
    fclose(out);
  }
  Timing_FileIO += clock() - tstart;
}

/*============================================================*/

void coupleconstant(doublecomplex mrel,char which,doublecomplex *res)
{
  doublecomplex CoupleConstant[3];
  doublecomplex tempa,tempb,cm,m2,t1;
  double temp,V,b1,b2,b3;
  int i,j;
  double S,prop2[3];
  int asym;           /* whether polarizability is asymmetric */
 
  asym = (PolRelation==CLDR || PolRelation==SOrd);
  if (asym) j=3;
  else j=1;
  if (PolRelation==LDR || PolRelation==CLDR) {
    b1=LDR_B1;
    b2=LDR_B2;
    b3=LDR_B3;
  }
  if (PolRelation==SOrd) {
    b1=SO_B1;
    b2=SO_B2;
    b3=SO_B3;
  }
  /* calculate the CM couple constant CC=(3V/4pi)*(m^2-1)/(m^2+2) */
  V=gridspace*gridspace*gridspace;     /* volume of one dipole */
  temp = (3*V)/(4*PI);
  cSquare(mrel,m2);                    /* m2=m^2 */
  tempa[re] = m2[re] - 1.0;
  tempa[im] = tempb[im] = m2[im];                      
  tempb[re] = m2[re] + 2.0;
  cDiv(tempa,tempb,CoupleConstant[0]);
  CoupleConstant[0][re] *= temp;
  CoupleConstant[0][im] *= temp;
  
  if (PolRelation!=CM) {
    if (PolRelation!=RADCOR) {
        for (i=0;i<3;i++) prop2[i]=prop[i]*prop[i];
        if (PolRelation==LDR) {
          if (avg_inc_pol) S=0.5*(1-DotProd(prop2,prop2));
          else {  
            if (which=='X') incPol=incPolX; 
            else incPol=incPolY;
            S=prop2[0]*incPol[0]*incPol[0]+prop2[1]*incPol[1]*incPol[1]+prop2[2]*incPol[2]*incPol[2];
          }
       }
    }
    memcpy(cm,CoupleConstant[0],sizeof(doublecomplex));
    for (i=0;i<j;i++) {
      t1[re]=0.0;                     
      t1[im]=2*kd*kd*kd/3;                       /* t1=2/3*i*kd^3         */
      if (PolRelation!=RADCOR) {
        if (PolRelation!=LDR) S=prop2[i];
        t1[re]+=(b1+(b2+b3*S)*m2[re])*kd*kd;   /* t1+=(b1+(b2+b3*S)*m^2)*kd^2  */
        t1[im]+=(b2+b3*S)*m2[im]*kd*kd;
      }
      t1[re]/=V;     
      t1[im]/=V;
      cMultSelf(t1,cm);
      t1[re]=1-t1[re];
      t1[im]=-t1[im];
      cDiv(cm,t1,CoupleConstant[i]);        /* CC[i]=cm/(1-(cm/V)*t1) */
    }
  }
  if (asym) {
    if (!orient_avg) {
      printz("CoupleConstant:{%.8g%+.8gi, %.8g%+.8gi, %.8g%+.8gi}\n",CoupleConstant[0][re],CoupleConstant[0][im],
    	CoupleConstant[1][re],CoupleConstant[1][im],CoupleConstant[2][re],CoupleConstant[2][im]);
      fprintz(logfile, "CoupleConstant:{%.8g%+.8gi, %.8g%+.8gi, %.8g%+.8gi}\n",CoupleConstant[0][re],CoupleConstant[0][im],
    	CoupleConstant[1][re],CoupleConstant[1][im],CoupleConstant[2][re],CoupleConstant[2][im]);
    }
  }
  else {
    memcpy(CoupleConstant[1],CoupleConstant[0],sizeof(doublecomplex));
    memcpy(CoupleConstant[2],CoupleConstant[0],sizeof(doublecomplex));
    if (!orient_avg) {
      printz("CoupleConstant:%.8g%+.8gi\n",CoupleConstant[0][re],CoupleConstant[0][im]);
      fprintz(logfile,"CoupleConstant:%.8g%+.8gi\n",CoupleConstant[0][re],CoupleConstant[0][im]);
    }
  }
  memcpy(res,CoupleConstant,3*sizeof(doublecomplex));
}

/*============================================================*/

void init_cc(char which)
    /* calculate cc and cc_sqrt */
{
  int i,j;

  for(i=0;i<Nmat;i++) {
    coupleconstant(ref_index[i],which,cc[i]);
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
    if ((tab_index = imatrix(1,TAB_RMAX,0,TAB_RMAX)) == NULL)
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
  free_imatrix(tab_index,1,TAB_RMAX,0);
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

void save_mueller(double *muel)
   /* saves mueller matrix (averaged) to file */
{
  FILE *mueller;
  char stringbuffer[200];
  int i,j;
  double dtheta,theta;
  clock_t tstart;

  tstart=clock();

  strcpy(stringbuffer,directory);
  strcat(stringbuffer,"/mueller");
  
  if ((mueller = fopen (stringbuffer, "w"))==NULL)
    LogError(EC_ERROR,ONE,POSIT,"File 'mueller' write error");

  fprintf(mueller,"theta s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44\n");
  
  dtheta = 180.0 / ((double)(nTheta-1));
  for (i=0;i<nTheta;i++) {
    theta=i*dtheta;
    fprintf(mueller,"%.2f",theta);
    for (j=0;j<16;j++) fprintf(mueller," %.7E",muel[16*i+j]);
    fprintf(mueller,"\n");
  }
  fclose(mueller);

  Timing_FileIO += clock() - tstart;
}

/*==============================================================*/

void save_CS(double Cext, double Cabs)
   /* save calculated crossections (averaged) to file */
{        
  FILE *CCfile;
  char stringbuffer[200];
  double a_eq,G;
  clock_t tstart;

  tstart=clock();
  
  strcpy(stringbuffer,directory);
  strcat(stringbuffer,"/CrossSec");

  if ((CCfile=fopen(stringbuffer,"w"))==NULL)
    LogError(EC_ERROR,ONE,POSIT,"Could not write to file '%s'",stringbuffer);

  a_eq = pow((0.75/PI)*nvoid_Ndip,0.333333333333)*gridspace;
  G = PI*a_eq*a_eq;
  fprintf(CCfile,"x=%.8lg\n\n",WaveNum*a_eq);

  printf("Cext\t= %12.8lg\nQext\t= %12.8lg\n",Cext,Cext/G);
  fprintf(CCfile,"Cext\t= %12.8lg\nQext\t= %12.8lg\n",Cext,Cext/G);
  printf("Cabs\t= %12.8lg\nQabs\t= %12.8lg\n",Cabs,Cabs/G);
  fprintf(CCfile,"Cabs\t= %12.8lg\nQabs\t= %12.8lg\n",Cabs,Cabs/G);
  
  fclose(CCfile);

  Timing_FileIO += clock() - tstart;
}

/*============================================================*/

void calculate_one_orientation(double *res /* where to put results when doing orientation averaging */)
     /* performs calculation for one orientation */
{     
  int i;
  clock_t tstart;
  
  if (orient_avg) {
    alph_deg=0;
    init_rotation();

    printz("\nORIENTATION STEP beta=%g gamma=%g\n",bet_deg,gam_deg);
    fprintz(logfile,"\nORIENTATION STEP beta=%g gamma=%g\n",bet_deg,gam_deg);
  }
  
  /* calculate scattered field for y - polarized incident light */
  printz("\nhere we go, calc Y\n\n");
  if (!orient_avg) fprintz(logfile,"\nhere we go, calc Y\n\n");
  init_cc('Y');
  if (symR==true && all_dir==false && scat_grid==false) {
    CalculateE('Y',PAR_AND_PER);
  }
  else { /* no rotational symmetry */
  /* in case of all_dir we run twice to get the full electric field */
  /* with incoming light in X and Y direction. In case of rotational */
  /* symmetry this is not needed but requires lots more programming */
  /* so we leave this optimization to a later time. */

    CalculateE('Y',NORMAL);

    printz("\nhere we go, calc X\n\n");
    if (!orient_avg) fprintz(logfile,"\nhere we go, calc X\n\n");
    if(PolRelation==LDR && !avg_inc_pol) init_cc('X');

    CalculateE('X',NORMAL);
  }
  D("CalculateE complete");
  mueller_matrix();
  D("mueller_matrix complete");
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
  int i;
  double si;

  Bcast_orient(&beta_i,&gamma_i,&finish_avg);
  if (finish_avg) return;
  
  bet_deg=beta_int.val[beta_i];
  gam_deg=gamma_int.val[gamma_i];
  calculate_one_orientation(res);
}

/*============================================================*/

void Calculator (void)
{
  int temp_int,ntheta_mod;
  doublecomplex mrel;		/* relative refractive index of particle */
  double temp, memtot,memmax;
  doublecomplex tempa, tempb;
  int i,i_or,j_or;

  extern int boxX,boxY,boxZ;
  extern int symX, symY;
  extern int jagged;
  extern double gridspace;
  extern int IntRelation;

  char bufstr[200];
  double *out;
  FILE *romb_file;

  /* initialize some variables */
  WaveNum     = TWO_PI / lambda;
  kd = TWO_PI / dpl;

  block_theta=16*nTheta;
  finish_avg=false;

  if (IntRelation == SOrd) ReadTables();
  
  D("init_Dmatrix started");
  init_Dmatrix();
  D("init_Dmatrix complete");
  
  /* allocate all the memory */
  if (prognose==false) if ((x = dCvector(nlocalRows)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc x");
  memory+=nlocalRows*sizeof(doublecomplex);

  if (prognose==false) if ((r = dCvector(nlocalRows))== NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc r");
  memory+=nlocalRows*sizeof(doublecomplex);

  if (prognose==false) if ((p = dCvector(nlocalRows)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc p");
  memory+=nlocalRows*sizeof(doublecomplex);

  if (prognose==false) if ((Einc = dCvector(nlocalRows)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Einc");
  memory+=nlocalRows*sizeof(doublecomplex);

  if (prognose==false) if ((Avecbuffer = dCvector(nlocalRows)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Avecbuffer");
  memory+=nlocalRows*sizeof(doublecomplex);

  if (IterMethod==IT_BICGSTAB || IterMethod==IT_QMR_CS) {
    /* additional vectors for iterative methods */
    if (prognose==false) if ((vec1 = dCvector(nlocalRows)) == NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc vec1");
    memory+=nlocalRows*sizeof(doublecomplex);
    if (prognose==false) if ((vec2 = dCvector(nlocalRows)) == NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc vec2");
    memory+=nlocalRows*sizeof(doublecomplex);
    if (prognose==false) if ((vec3 = dCvector(nlocalRows)) == NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc vec3");
    memory+=nlocalRows*sizeof(doublecomplex);
  }

  if (yzplane) {
    /* more exact definition is given at every call of CalculateE */
    if ((symY && symX) || orient_avg)
      ntheta_mod=nTheta;
    /* if not symmetric calculate for +- theta */
    else ntheta_mod=2*(nTheta-1);
    if (prognose == false) {
      /* Epar and Eper are allocated as parts of one vector,
         so they can be accumulated simultaneously */
      if ((Epar = dCvector(2*ntheta_mod)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc Epar");
      Eper=Epar+ntheta_mod; /* may be corrected in CalculateE */
#ifdef PARALLEL
      /* buffer for accumulate operation */
      if ((Eparper_buffer = dvector(0,4*ntheta_mod-1)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc Eparper_buffer");
#endif
    }
    memory+=4*(nTheta-1)*sizeof(doublecomplex);
#ifdef PARALLEL
    memory+=8*(nTheta-1)*sizeof(double);  /* extra memory for buffer */
#endif
  }

  if (all_dir) {
    ReadAlldirParms(alldir_parms);
    /* calculate size of vectors; 4 - because first it is used to store
       per and par components of the field, and only afterwards squares */
    temp_int=4*theta_int.N*phi_int.N;
    if (prognose == false) {
      if ((E2_alldir = dvector(0,temp_int-1)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc E2_alldir");
#ifdef PARALLEL
      /* buffer for accumulate operation */
      if ((E2_alldir_buffer = dvector(0,temp_int-1)) == NULL)
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
    if (prognose == false) {
      if ((EgridX = dCvector(temp_int)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc EgridX");
      if ((EgridY = dCvector(temp_int)) == NULL)
        LogError(EC_ERROR,ALL,POSIT,"Could not malloc EgridY");
      if (phi_integr && ringid==ROOT)
        if ((muel_phi = dvector(0,16*angles.phi.N+1)) == NULL)
          LogError(EC_ERROR,ONE,POSIT,"Could not malloc muel_phi");
        if ((muel_phi1 = dvector(0,16*angles.phi.N+1)) == NULL)
          LogError(EC_ERROR,ONE,POSIT,"Could not malloc muel_phi1");
#ifdef PARALLEL
      /* buffer for accumulate operation */
      if ((Egrid_buffer = dvector(0,2*temp_int-1)) == NULL)
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
    if (prognose==false && ringid==ROOT) {
      if ((ampl_alpha = dCvector(4*nTheta*alpha_int.N)) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc ampl_alpha");
      if ((muel_alpha = dvector(0,block_theta*alpha_int.N+1)) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc muel_alpha");
      muel_alpha+=2;
      if ((out = dvector(0,block_theta+1)) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc out");
    }
    memory+=4*nTheta*sizeof(doublecomplex)*alpha_int.N+sizeof(double)*(block_theta*(alpha_int.N+1)+4);
  }
  /* estimate of the memory (only the fastest scaling part):
       MatVec - 288*Ndip (more exactly gridX*gridY*gridZ*72)
       others - nvoid_Ndip*271(+144 for BiCGStab)
     PARALLEL: above is total; division over processors of MatVec is uniform,
               others - according to local_nvoid_Ndip  */
  memtot=memory/MBYTE;
  AccumulateMax(&memtot,&memmax);
  printz("Total memory usage:%.1f Mb\n",memtot);
  fprintz(logfile,"Total memory usage:%.1f Mb\n",memtot);
#ifdef PARALLEL
  fprintz(logfile,"Maximum memory usage of single processor:%.1f Mb\n",memmax);
#endif
  if (prognose==true) stop(0);

  if (!orient_avg) alpha_int.N=1;

  Timing_Init = clock() - tstart_main;

  if (orient_avg) {
    if (ringid==ROOT) {
      sprintf(bufstr,"%s/log_orient_avg",directory);
      Romberg2D(parms,orient_integrand,block_theta+2,out,bufstr);
      finish_avg=true;
      Bcast_orient(&finish_avg,&finish_avg,&finish_avg);
      save_mueller(&out[2]);
      save_CS(out[0],out[1]);
    }
    else while (!finish_avg) orient_integrand(0,0,NULL);
  }
  else calculate_one_orientation(NULL);

  /* save effective volume in_mie */
  if (ringid==ROOT) save_in_mie();

  /* tidy up */
  if (IntRelation == SOrd) FreeTables();
  free_FFT_Dmat();
  free_dCvector(x);
  free_dCvector(r);
  free_dCvector(p);
  free_dCvector(Einc);
  free_dCvector(Avecbuffer);
  if (IterMethod==IT_BICGSTAB || IterMethod==IT_QMR_CS) {
    free_dCvector(vec1);
    free_dCvector(vec2);
    free_dCvector(vec3);
  }
  if (yzplane) {
    free_dCvector(Epar);
#ifdef PARALLEL
    free_dvector(Eparper_buffer,0);
#endif
  }
  if (all_dir) {
    free(theta_int.val);
    free(phi_int.val);
    free_dvector(E2_alldir,0);
#ifdef PARALLEL
    free_dvector(E2_alldir_buffer,0);
#endif
  }
  if (scat_grid) {
    free(angles.theta.val);
    free(angles.phi.val);
    free_dCvector(EgridX);
    free_dCvector(EgridY);
    if (phi_integr && ringid==ROOT) {
      free(muel_phi);
      free(muel_phi1);
    }
#ifdef PARALLEL
    free_dvector(Egrid_buffer,0);
#endif
  }
  /* these 3 were allocated in make_particle */
  free_dvector(DipoleCoord,0);
  free(position);
  free(material);

  if (orient_avg) {
    if (ringid==ROOT) {
      free_dCvector(ampl_alpha);
      free_dvector(muel_alpha-2,0);
      free_dvector(out,0);
    }
    free(alpha_int.val);
    free(beta_int.val);
    free(gamma_int.val);
  }
}



