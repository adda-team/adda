/* FILE: calculator.c	
 * AUTH: Alfons Hoekstra
 * DESCR: All the initialization is done here before actually
 *        calculating internal fields
 *
 *        Currently is developed by Maxim Yurkin
 */
#include <stdio.h>
#include <time.h>
#include "cmplx.h"
#include "types.h"
#include "Romberg.h"
#include "const.h"
#include <math.h>
#include "comm.h"
#include "debug.h"
#include "read_cluster.h"

double dpl;
double LAMBDA;

extern doublecomplex *dCvector(int nl,int nh);
extern double **dmatrix(int nrl,int nrh,int ncl,int nch);
extern void free_dmatrix(double **m, int nrl, int nrh, int ncl);
extern int **imatrix (int nrl, int nrh, int ncl, int nch);

extern doublecomplex *x, *p, *r;
extern doublecomplex *buffer;	/* used to gather vectors in */
extern doublecomplex *Eper, *Epar, *Egrid, *Einc;
extern FILE *logfile;
extern int prognose;
extern doublecomplex ref_index[10];
extern Parms_1D parms[];
extern int yzplane;
extern int xzplane;
extern int all_dir;
extern char directory[200];
extern double prop[3];          /* propagation vector */
extern double incPolX[3],incPolY[3]; /* incident polarizations corresponding 
				        to X and Y for default incidence */
extern int orient_avg;
extern integr_parms alpha_int, beta_int, gamma_int;
extern double alph, bet, gam, alph_deg, bet_deg, gam_deg;

extern int shape;             /* shape of the particle */
extern int Nmat,mat_count[10];
extern int nvoidDip;
extern double gridspace;

extern int symR;
extern int PolRelation;

int nDip;             /* number of dipoles */
int nTheta;	      /* number of angles in scattering profile */
int nlocalRows;       /* number of local rows of decomposition */
doublecomplex *Avecbuffer; /* buffer to hold result of local Matvec */
double **DipoleCoord;   /* matrix to hold the coordinates of the dipoles */
short int *position;  /* position of the dipoles */
int *material;        /* material: index for cc */
double WaveNum;       /* wavenumber of incident light */
double eps;	      /* accuracy epsilon */
double lambda;                /* wavelenght of incident light */
int memory=0;         /* total memory usage in bytes */
double *incPol;		/* used polarization */
doublecomplex cc[10][3];   /* couple constants */
double kd;             /* k*d=2*PI/dpl */
double *tab1,*tab2,*tab3,*tab4,*tab5,*tab6,*tab7,*tab8,*tab9,*tab10; /* tables of integrals */
int **tab_index;      /* matrix for indexation of table arrays */

extern int avg_inc_pol;
doublecomplex CoupleConstant[3];

doublecomplex *ampl_alpha;  /* storing amplitude matrix for different values of alpha */
double *muel_alpha;        /* stroring mueller matrix for different values of alpha */
Parms_1D parms_alpha;
int block_theta;       /* size of one block of mueller matrix - 16*nTheta */

extern unsigned long TotalEval;
extern clock_t Timing_FileIO, Timing_Integration, Timing_Init;
extern clock_t tstart_main;

int finish_avg;
       
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
  r=(2.0*PI/dpl)*pow((0.75/PI)*nvoidDip,.333333333333);
  
  inner=0; outer=0; rin=r;
  if (shape==COATED) {
    inner=1;
    outer=0;
    rin=coat_ratio*r;
    rin=2.0*PI*pow((0.75/PI)*mat_count[1]/(dpl*dpl*dpl),.333333333333); /* !!!will not work properly for PARALLEL */
  }

  tstart=clock();

  for(i=0;i<Nmie;i++) {
    sprintz(buffer,"%s/in_mie%i",directory,i);
    
    if ((out=fopen (buffer, "wA"))==NULL) 
      LogError(EC_ERROR,ONE,POSIT,"File 'in_mie' write error");
    
    if (beamtype==PLANE || (beam_x0==0.0 && beam_y0==0.0 && beam_z0==0.0)) {
      fprintz(out,"%i\n",nTheta);
      fprintz(out,"%g\n%g\n",rin,r);
      fprintz(out,"%f\n%f\n%f\n%f\n",
	      ref_index[inner].r,ref_index[inner].i,
	      ref_index[outer].r,ref_index[outer].i);
      fprintz(out,"1.0\n%.5E\n",LAMBDA);
      if (beamtype==PLANE) {
	fprintz(out,"%.5E\n",1e6*LAMBDA);
	fprintz(out,"0\n");
      }
      else {
	fprintz(out,"%.5E\n",beam_w0);
	fprintz(out,"99\n");
      }
    } 
    else { /* focussed beam */
      fprintz(out,"%.7E\n",LAMBDA);
      fprintz(out,"%.7E\n",beam_w0);
      fprintz(out,"%.7E\n",r/PI*LAMBDA);
      fprintz(out,"%.7E\n",ref_index[0].r);
      fprintz(out,"%.7E\n1.0\n",ref_index[0].i);
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

doublecomplex *coupleconstant(doublecomplex mrel,char which)
{
  doublecomplex tempa,tempb,cm,m2,t1;
  double temp,V,b1,b2,b3;
  int i,j;
  double S,norm,prop2[3];
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
  cSquare(mrel,&m2);                    /* m2=m^2 */
  tempa.r = m2.r - 1.0;
  tempa.i = tempb.i = m2.i;                      
  tempb.r = m2.r + 2.0;
  cDiv(tempa,tempb,CoupleConstant);
  CoupleConstant->r *= temp;
  CoupleConstant->i *= temp;
  
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
    memcpy(&cm,CoupleConstant,sizeof(doublecomplex));
    for (i=0;i<j;i++) {
      t1.r=0.0;                     
      t1.i=2*kd*kd*kd/3;                       /* t1=2/3*i*kd^3         */
      if (PolRelation!=RADCOR) {
        if (PolRelation!=LDR) S=prop2[i];
        t1.r+=(b1+(b2+b3*S)*m2.r)*kd*kd;   /* t1+=(b1+(b2+b3*S)*m^2)*kd^2  */
        t1.i+=(b2+b3*S)*m2.i*kd*kd;
      }
      t1.r/=V;     
      t1.i/=V;                                
      cMult(t1,cm,&t1);
      t1.r=1-t1.r; 
      t1.i=-t1.i;
      cDiv(cm,t1,&(CoupleConstant[i]));        /* CC[i]=cm/(1-(cm/V)*t1) */
    }
  }
  if (asym) {
    if (!orient_avg) {
      printz("CoupleConstant:{%g%+gi, %g%+gi, %g%+gi}\n",CoupleConstant[0].r,CoupleConstant[0].i,
    	CoupleConstant[1].r,CoupleConstant[1].i,CoupleConstant[2].r,CoupleConstant[2].i);
      fprintz(logfile, "CoupleConstant:{%g%+gi, %g%+gi, %g%+gi}\n",CoupleConstant[0].r,CoupleConstant[0].i,
    	CoupleConstant[1].r,CoupleConstant[1].i,CoupleConstant[2].r,CoupleConstant[2].i);
    }
  }
  else {
    memcpy(&(CoupleConstant[1]),CoupleConstant,sizeof(doublecomplex));
    memcpy(&(CoupleConstant[2]),CoupleConstant,sizeof(doublecomplex));
    if (!orient_avg) {
      printz("CoupleConstant:%g%+gi\n",CoupleConstant[0].r,CoupleConstant[0].i);
      fprintz(logfile,"CoupleConstant:%g%+gi\n",CoupleConstant[0].r,CoupleConstant[0].i);
    }
  }
  return(CoupleConstant);
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
  /* allocate memory for tab_n */
  if ((tab_n = malloc(size*sizeof(double))) == NULL) 
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc tab_n");
  /* read file */
  strcpy(fname,TAB_PATH);
  strcat(fname,sh_fname);
  
  if ((ftab=fopen(fname,"r"))==NULL) 
    LogError(EC_ERROR,ALL,POSIT,"Could not open file '%s'",fname);
  
  for (i=0; i<size; i++) if (fscanf(ftab,"%lf\t",&(tab_n[i]))!=1) 
      LogError(EC_ERROR,ALL,POSIT,"Scan error in file '%s'. Probably file is too small",fname);

  if (!feof(ftab)) LogError(EC_WARN,ONE,POSIT,"File '%s' is longer than specified size (%d)",fname,size);
  fclose(ftab);
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
  /* printz("P[5,3]=%d (should be 41)\n",tab_index[5][3]); */
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
  
  if ((mueller = fopen (stringbuffer, "wA"))==NULL) 
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

  fprintf(logfile,"Total number of occupied dipoles %d\n",nvoidDip);
  a_eq = pow((0.75/PI)*nvoidDip,.333333333333)*gridspace;
  G = PI*a_eq*a_eq;
  fprintf(CCfile,"x=%lg\n\n",WaveNum*a_eq);

  printf("Cext\t= %12lg\nQext\t= %12lg\n",Cext,Cext/G);
  fprintf(CCfile,"Cext\t= %12lg\nQext\t= %12lg\n",Cext,Cext/G);
  printf("Cabs\t= %12lg\nQabs\t= %12lg\n",Cabs,Cabs/G);
  fprintf(CCfile,"Cabs\t= %12lg\nQabs\t= %12lg\n",Cabs,Cabs/G);
  
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
  
  for(i=0;i<Nmat;i++) memcpy(&cc[i][0],coupleconstant(ref_index[i],'Y'),3*sizeof(doublecomplex));

  if (symR==true&&all_dir==false) {
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
    if(PolRelation==LDR && !avg_inc_pol) for(i=0;i<Nmat;i++) 
	memcpy(&cc[i][0],coupleconstant(ref_index[i],'X'),3*sizeof(doublecomplex));
    
    CalculateE('X',NORMAL);
  }
  D("CalculateE complete");
  mueller_matrix();
  D("mueller_matrix complete");
  if (ringid==0 && orient_avg) {
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
  int tempdip;
  double dips_lambda;		/* number of dipoles per wavelenght */
  doublecomplex mrel;		/* relative refractive index of particle */
  double temp;
  doublecomplex tempa, tempb;
  int i,i_or,j_or;

  extern int boxX,boxY,boxZ;
  extern double eps0;
  extern int jagged;
  extern double gridspace;
  extern int IntRelation;

  char bufstr[200];
  double *out;
  FILE *romb_file;

  /* initialize some variables */
  lambda      = LAMBDA;
  dips_lambda = dpl;
  WaveNum     = 2.0 * PI / lambda;
  kd = 2.0 * PI / dpl;
  eps=eps0;
  
  block_theta=16*nTheta;
  finish_avg=false;
 
  /* determine the number dipoles in this processor */
  /* if (nDip%nprocs == 0) local_Ndip = nDip/nprocs;
  else if (ringid < nDip%nprocs) local_Ndip = nDip/nprocs + 1;
  else local_Ndip = nDip/nprocs; */  /* local_Ndip is used instead - defined in calculator */
  
  if (IntRelation == SOrd) ReadTables();
  
  D("init_Dmatrix started");
  init_Dmatrix(boxX,boxY,boxZ,lambda,dips_lambda,WaveNum);
  D("init_Dmatrix complete");
  
  nlocalRows = 3 * local_Ndip;
  
  if (prognose==false) if ((x = dCvector(0, nlocalRows -1)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc x");
  memory+=sizeof(doublecomplex)*nlocalRows;
  
  if (prognose==false) if ((r = dCvector(0, nlocalRows -1))== NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc r");
  memory+=sizeof(doublecomplex)*nlocalRows;

  if (prognose==false) if ((p = dCvector(0, nlocalRows -1)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc p");
  memory+=sizeof(doublecomplex)*nlocalRows;

  if (prognose==false) if ((Einc = dCvector(0, nlocalRows -1)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Einc");
  memory+=sizeof(doublecomplex)*nlocalRows;

  /* since we might need to hold nprocs*(nDip/nprocs+1) elements in */
  /* the buffer, we increase its size (was 3*nDip) */
  if (prognose==false) if ((buffer = dCvector(0, nlocalRows -1)) == NULL) 
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc buffer");
  memory+=nlocalRows*sizeof(doublecomplex);
  Avecbuffer=buffer;
  
  if (yzplane || xzplane)
    {
      if (prognose == false)
	{
	  if ((Eper = dCvector(0, 2*nTheta-3)) == NULL)
	    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Eper");
	  
	  if ((Epar = dCvector(0, 2*nTheta-3)) == NULL)
	    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Epar");
	}
      memory+=4*(nTheta-1)*sizeof(doublecomplex);
    }
  
  if (all_dir)
    {
      set_Parms();

      if (prognose == false) 
        if ((Egrid = dCvector(0,3*parms[THETA].Grid_size*parms[PHI].Grid_size*sizeof(doublecomplex))) == NULL) 
	  LogError(EC_ERROR,ALL,POSIT,"Could not malloc Egrid");
      memory+=3*parms[THETA].Grid_size*parms[PHI].Grid_size;
    }

  if (prognose==false) if ((DipoleCoord = dmatrix(0, local_Ndip, 0, 2)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc DipoleCoord");
  memory+=(sizeof(double)*3+sizeof(int))*local_Ndip;
  
  if (prognose==false) if ((position = (short int *) malloc(3*sizeof(short int)*local_Ndip)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc position");
  memory+=3*sizeof(short int)*local_Ndip;
  
  if (prognose==false) if ((material = (int *) malloc(sizeof(int)*local_Ndip)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc material");
  memory+=sizeof(char)*local_Ndip;

  if (orient_avg) {
    if (prognose==false && ringid==0) {
      if ((ampl_alpha = (doublecomplex *) malloc(sizeof(doublecomplex)*4*nTheta*alpha_int.N)) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc ampl_alpha");
      if ((muel_alpha = (double *) malloc(sizeof(double)*(block_theta*alpha_int.N+2))) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc muel_alpha");
      muel_alpha=&muel_alpha[2];
      if ((out = (double *) malloc(sizeof(double)*(block_theta+2))) == NULL)
        LogError(EC_ERROR,ONE,POSIT,"Could not malloc out");
    }
    memory+=4*nTheta*sizeof(doublecomplex)*alpha_int.N+sizeof(double)*(block_theta*(alpha_int.N+1)+4);
  }

  printz("Total memory usage for all matrices:%.1f Mb\n",memory/1048576.0);
  fprintz(logfile,"Total memory usage for all matrices:%.1f Mb\n",memory/1048576.0);
  if (prognose==true) {
    double flop,m;
    extern int gridX,gridY,gridZ,Ngrid;
    m=sqrt(ref_index[0].r*ref_index[0].r+ref_index[0].i*ref_index[0].i);
    printz("Prognose: total memory: %.1f Mb\n",memory/1048576.0);
    printz("Prognose: required ram: %.1f Mb\n",sizeof(doublecomplex)*Ngrid/1048576.0);
    flop=Ngrid*log((double) Ngrid)*(m-1)/dpl;
    printz("Prognose: required time: %.2g flop\n",flop);
    stop(1);
  }
  /* generate the matrix with dipole coordinates, do this in all processors */
  D("make_particle started");
  tempdip = make_particle (DipoleCoord,true,shape,jagged);
  D("make_particle complete");
  for(i=0;i<Nmat;i++) printz("material:%i has locally %i dipoles\n",i,mat_count[i]);
  
  if (nDip != tempdip) LogError(EC_ERROR,ALL,POSIT,"Wrong number of dipoles - %d",tempdip);

  if (!orient_avg) alpha_int.N=1;

  Timing_Init = clock() - tstart_main;
 
  if (orient_avg) {
    if (ringid==0) {
      sprintf(bufstr,"%s/romb_log",directory);
      Romberg2D(parms,&orient_integrand,block_theta+2,out,bufstr);
      finish_avg=true;
      Bcast_orient(&finish_avg,&finish_avg,&finish_avg);
      save_mueller(&out[2]);
      save_CS(out[0],out[1]);
    }
    else while (!finish_avg) orient_integrand(0,0,NULL);
  }
  else calculate_one_orientation(NULL);

  /* save effective volume in_mie */
  if (ringid==0) save_in_mie();

  /* tidy up */
  /* lets not do this for a while....
   * runing MPI with this sometimes gives strange results, 
   * no time to find out what goes wrong....
  free_FFT_Dmat();
  free(p);
  free(r);
  free(buffer);
  free(Avecbuffer);
  free(position);

  if (yzplane || xzplane)
    {
      free(Epar);
      free(Eper);
    }
  if (all_dir)
    free(Egrid);

  free(Einc);
  free(x);
  free_dmatrix(DipoleCoord,0, local_Ndip, 0);
  free(material);
  so, until this point we have skipped stuff..... */
}



