/* FILE: CD2main.c
 * AUTH: Alfons Hoekstra
 * DESCR: Main. Reads command line parameters, initialize variables, 
 *        output log.
 *
 *        sequential version, Michel Grimminck Jan 1995 
 *
 *        Currently is developed by Maxim Yurkin
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cmplx.h"
#include "types.h"
#include "const.h"
#include <math.h>
#include <limits.h>
#include "comm.h"
#include "debug.h"
#include "crosssec.h"

#define TO_SEC(p) (p / (double) CLOCKS_PER_SEC)

/* global variables */
extern int nDip;

int shape;                       /* particle shape definition */
int boxX,boxY,boxZ;              /* grid of calculation box */
int symX,symY,symZ;              /* symmetries in x y and z direction of
				    particle */
int symR;                        /* symmetry if 90 deg. rotation about z axes */
int NoSymmetry;                  /* do not use particle symmetries */
int symmetry_enforced;           /* enforce use of symmetry */
double prop_0[3],prop[3];         /* propagation vector of incident wave (initial and rotated) */
double incPolX_0[3],incPolX[3],incPolY_0[3],incPolY[3];  /* incident polariztions */
int PolRelation;                 /* type of polarization relation */
int ScatRelation;                /* type of formulae for scattering quantities */
int IntRelation;                 /* type of formula for interaction term */
int prognose;                    /* make a prognose about needed ram,virtmem */
int store_field;                 /* save internal fields to pbm file */
int store_int_field;		 /* save full internal fields to file */
double eps0;                     /* stop criterium */
char directory[200];             /* directory to save data in */
FILE *logfile;
double gridspace;                /* inter-dipole distance */
double sizeX;
int profile[10];
int maxiter;                     /* the maximal number of iterations */
int jagged;                      /* particle is jagged */
int beamtype;                    /* barton,davis,lminus */
double beam_w0,beam_x0,beam_y0,beam_z0=0.0;  /* beam properties (meters) */
doublecomplex ref_index[10];          /* a set of refraction indexes */
int Nmat;                        /* number of different refraction indexes */
int mat_count[10];               /* number of dipoles in a material */
double coat_ratio=1.0;
double coat_x=0.0,coat_y=0.0,coat_z=0.0;
double append_ratio, ratio_x, ratio_y, ratio_z, cratio_x, cratio_y, cratio_z; /* for limphocytes particle*/
double xc0,yc0,zc0,xc1,yc1,zc1,xc2,yc2,zc2,xc3,yc3,zc3,xc4,yc4,zc4,
       ratio0, ratio1, ratio2, ratio3, ratio4, ratio5, delta, inc_ratio; /* for leucocyte2 particle */ 
double diskratio;
double ellipsX, ellipsY, ellipsZ; /* for ellipsoidal particle */
double aspect_r, betaY, betaZ; /* for rotatable erithrocyte, for rotatable smoothed disk and rotatable oblate spheroid */
double aspect_ratio, alphaY, alphaX; /* for stick particle */

int yzplane = true;              /* Calculate the field in the yz-plane */ 
int xzplane = false;		 /* Calculate the field in the xz-plane */
int all_dir = false;             /* Calculate the field for all directions
				    on a theta-phi grid */
int store_all_dir = false;       /* Store the field for all directions */
int calc_Cext = false;           /* Calculate the extinction cross-section */ 
int calc_Cabs = false;           /* Calculate the absorption cross-section */
int calc_Csca = false;           /* Calculate the scattering cross-section
				    by integration */
int calc_vec = false;            /* Calculate the unnormalized 
				    asymmetry-parameter */
int calc_asym = false;           /* Calculate the asymmetry-parameter */
int calc_mat_force = false;      /* Calculate the scattering force 
				    by matrix-evaluation */
int store_force = false;         /* Write radiation pressure per dipole
				    to file */
int reduced_FFT;                 /* reduced number of storage for FFT, when matrix is symmetric */ 

char Romb_parms[200];
char avg_parms[200];
char aggregate_file[200];

double alph_deg, bet_deg, gam_deg;  /* Euler angles of particle orientation in degrees */
int orient_avg;               /* whether to use orientation averaging*/

char avg_string[200];          /* string for output of function that reads averaging parameters */

int avg_inc_pol;              /* whether to average polarization prescription over incident polarization */
                                 
/*
 * start adress of vectors, used by both the calculator and the router. 
 * The vectors are initialized by the Calculator
 */

doublecomplex *buffer, *x, *r, *p, *Epar, *Eper, *Einc, *Egrid;

/* Some timing variables */
clock_t Timing_TotalTime,Timing_FFT_Init,
  Timing_OneIter, Timing_OneIterCalc, Timing_OneIterComm, 
  Timing_EField, Timing_calc_EField, Timing_comm_EField,
  Timing_EField_ad, Timing_calc_EField_ad, Timing_comm_EField_ad, 
  Timing_FileIO, Timing_Integration, Timing_EFieldPlane,
  Timing_IntField, Timing_ScatQuan, Timing_Init, Timing_Particle,
  Timing_IntFieldOne,Timing_InitIter;

unsigned long TotalIter,TotalEval,TotalEFieldPlane;

clock_t  tstart_main;

int equivalent=true; /* whether to assume that min and max values of alpha,      */
                     /* gamma (phi) are completely equivalent */ 

/*============================================================*/

void main (int argc, char *argv[])
{
  char          *pname,sbuffer[500],logname[100],name[100],shapename[20];
  int		i,Nexp,Narg;
  FILE          *Nexpfile;
  time_t        start,end;
  char          procname[50];
  double 	temp,doubTmp; /*buffer variables*/
  
  extern double LAMBDA;	/* defined in calculator.c */
  extern double dpl;
  extern int    nDip,nTheta;
  
 
  time(&start);
  tstart_main = clock();

  for(i=0;i<10;i++) profile[i]=0;
  for(i=0;i<10;i++) {ref_index[i].r=1.0; ref_index[i].i=0.0;}
  /* defaults */
  prop_0[0]=0;        /*by default beam propagates along z-axis*/
  prop_0[1]=0;	
  prop_0[2]=1;
  directory[0]=0;
  LAMBDA=2*PI;
  ref_index[0].r=1.05;
  ref_index[0].i=0.0;
  boxX=boxY=boxZ=UNDEF;
  sizeX=UNDEF;
  dpl=10;
  strcpy(name,"run");
  nTheta=UNDEF; 
  eps0=/*1.0e-12*/1.0e-10;
  shape=SPHERE; strcpy(shapename,"sphere");
  store_field=false;
  store_int_field=false;
  PolRelation=LDR;
  ScatRelation=DRAINE;
  IntRelation=POINT;
  NoSymmetry=false;
  symmetry_enforced=false;
  prognose=false;
  maxiter=INT_MAX;
  jagged=1;
  procname[0]=0;
  beamtype=PLANE;
  strcpy(Romb_parms,"Romb5.dat");
  strcpy(avg_parms,"avg_params.dat");
  orient_avg=false;
  alph_deg=bet_deg=gam_deg=0.0;
  /*#ifdef PARALLEL*/
  Nmat=2;
  /*#else
  Nmat=1;
#endif*/
  test_interrupt();
  /* read command line */
  i=1; if (argc>1) do {
    /* get number of arguments */
    Narg=0;
    do {
      Narg++;
      if ((i+Narg)>=argc) break;
      if (argv[i+Narg][0]=='-') break;
    } while(true);
    Narg--;
    
    if (strcmp(argv[i],"-lambda")==0) {
      sscanf(argv[++i],"%lf",&LAMBDA);
    }
    else if (strcmp(argv[i],"-grid")==0) {
      if (Narg!=1 && Narg!=3) {
	printz("Illegal number of arguments (%d) to -grid option\n",Narg);
	stop(1);
      }
      sscanf(argv[++i],"%i",&boxX);
      if (Narg==3) {
        sscanf(argv[++i],"%i",&boxY);
        sscanf(argv[++i],"%i",&boxZ);
      }
    }
    else if (strcmp(argv[i],"-m")==0) {
      int j;
      
      if (Narg!=2 && Narg!=4) {
	printz("Illegal number of arguments (%d) to -m option\n",Narg);
	stop(1);
      }
      sscanf(argv[++i],"%lf",&ref_index[0].r);
      sscanf(argv[++i],"%lf",&ref_index[0].i);
      if (Narg==2) {
        ref_index[1].r=1.0;
        ref_index[1].i=0.0;
      }
      else {
        sscanf(argv[++i],"%lf",&ref_index[1].r);
        sscanf(argv[++i],"%lf",&ref_index[1].i);  
      }
      Narg=4;
      Nmat=(Narg+1)/2;
      /*#ifdef PARALLEL*/
      Nmat++;
      /*#endif*/
      /* if (Nmat>10) {printz("Too many materials\n"); stop(1);} */
    }
    else if (strcmp(argv[i],"-dpl")==0) {
      if (Narg!=1) {
        printz("Illegal number of arguments (%d) to -dpl option\n",Narg);
        stop(1);
      }
      sscanf(argv[++i],"%lf",&dpl);
    }
    else if (strcmp(argv[i],"-eps")==0) {
      if (Narg!=1) {
        printz("Illegal number of arguments (%d) to -eps option\n",Narg);
        stop(1);
      }
      sscanf(argv[++i],"%lf",&doubTmp);
      eps0=pow(10,-doubTmp);
    }
    else if (strcmp(argv[i],"-cylinder")==0) {
      if (Narg!=2) {
        printz("Illegal number of arguments (%d) to -particle option\n",Narg);
        stop(1);
      }
      shape=CYLINDER; strcpy(shapename,"cylinder");
      sscanf(argv[++i],"%lf",&aspect_ratio);
      sscanf(argv[++i],"%lf",&append_ratio);
    }
    else if (strcmp(argv[i],"-lymphocyte1")==0) {
      if (Narg!=5) {
        printz("Illegal number of arguments (%d) to -particle option\n",Narg);
        stop(1);
      }
      shape=LYMPHOCYTE1; strcpy(shapename,"lymphocyte1");
      coat_x=0.0; coat_y=0.0; coat_z=0.0;
      sscanf(argv[++i],"%lf",&coat_ratio);
      sscanf(argv[++i],"%lf",&coat_x);
      sscanf(argv[++i],"%lf",&coat_y);
      sscanf(argv[++i],"%lf",&coat_z);
      sscanf(argv[++i],"%lf",&append_ratio);
    }
    else if (strcmp(argv[i],"-lymphocyte2")==0) {
      if (Narg!=11) {
        printz("Illegal number of arguments (%d) to -particle option\n",Narg);
        stop(1);
      }
      shape=LYMPHOCYTE2; strcpy(shapename,"lymphocyte2");
      sscanf(argv[++i],"%lf",&coat_ratio);
      sscanf(argv[++i],"%lf",&coat_x);
      sscanf(argv[++i],"%lf",&coat_y);
      sscanf(argv[++i],"%lf",&coat_z);
      sscanf(argv[++i],"%lf",&append_ratio);
      sscanf(argv[++i],"%lf",&ratio_x);
      sscanf(argv[++i],"%lf",&ratio_y);
      sscanf(argv[++i],"%lf",&ratio_z);
      sscanf(argv[++i],"%lf",&cratio_x);
      sscanf(argv[++i],"%lf",&cratio_y);
      sscanf(argv[++i],"%lf",&cratio_z);
    }
    else if (strcmp(argv[i],"-leucocyte2")==0) {
      if (Narg!=22) {
        printz("Illegal number of arguments (%d) to -particle option\n",Narg);
        stop(1);
      }
      shape=LEUCOCYTE2; strcpy(shapename,"leucocyte2");
      sscanf(argv[++i],"%lf",&ratio0);
      sscanf(argv[++i],"%lf",&xc0);
      sscanf(argv[++i],"%lf",&yc0);
      sscanf(argv[++i],"%lf",&zc0);
      sscanf(argv[++i],"%lf",&ratio1);
      sscanf(argv[++i],"%lf",&xc1);
      sscanf(argv[++i],"%lf",&yc1);
      sscanf(argv[++i],"%lf",&zc1);
      sscanf(argv[++i],"%lf",&ratio2);
      sscanf(argv[++i],"%lf",&xc2);
      sscanf(argv[++i],"%lf",&yc2);
      sscanf(argv[++i],"%lf",&zc2);
      sscanf(argv[++i],"%lf",&ratio3);
      sscanf(argv[++i],"%lf",&xc3);
      sscanf(argv[++i],"%lf",&yc3);
      sscanf(argv[++i],"%lf",&zc3);
      sscanf(argv[++i],"%lf",&ratio4);
      sscanf(argv[++i],"%lf",&xc4);
      sscanf(argv[++i],"%lf",&yc4);
      sscanf(argv[++i],"%lf",&zc4);
      sscanf(argv[++i],"%lf",&inc_ratio);
      sscanf(argv[++i],"%lf",&delta);
    }
    else if (strcmp(argv[i],"-rbc_rot")==0) {
      if (Narg!=3) {
        printz("Illegal number of arguments (%d) to -particle option\n",Narg);
        stop(1);
      }
      shape=RBC_ROT; strcpy(shapename,"rbc_rot");
      sscanf(argv[++i],"%lf",&aspect_r);
      sscanf(argv[++i],"%lf",&betaY);
      sscanf(argv[++i],"%lf",&betaZ);
    }
    else if (strcmp(argv[i],"-sdisk_rot")==0) {
      if (Narg!=3) {
        printz("Illegal number of arguments (%d) to -particle option\n",Narg);
        stop(1);
      }
      shape=SDISK_ROT; strcpy(shapename,"sdisk_rot");
      sscanf(argv[++i],"%lf",&aspect_r);
      sscanf(argv[++i],"%lf",&betaY);
      sscanf(argv[++i],"%lf",&betaZ);
    }  
    else if (strcmp(argv[i],"-spheroid_rot")==0) {
      if (Narg!=3) {
        printz("Illegal number of arguments (%d) to -particle option\n",Narg);
        stop(1);
      }
      shape=SPHEROID_ROT; strcpy(shapename,"spheroid_rot");
      sscanf(argv[++i],"%lf",&aspect_r);
      sscanf(argv[++i],"%lf",&betaY);
      sscanf(argv[++i],"%lf",&betaZ);
    }   
    else if (strcmp(argv[i],"-stick")==0) {
      if (Narg!=4) {
        printz("Illegal number of arguments (%d) to -particle option\n",Narg);
	stop(1);
      }
      shape=STICK; strcpy(shapename,"stick");
      sscanf(argv[++i],"%lf",&aspect_ratio);
      sscanf(argv[++i],"%lf",&alphaY);
      sscanf(argv[++i],"%lf",&alphaX);
      sscanf(argv[++i],"%lf",&append_ratio);
    }
	
    else if (strcmp(argv[i],"-size")==0) {
      sscanf(argv[++i],"%lf",&sizeX);
      dpl=UNDEF;
    }
    else if (strcmp(argv[i],"-test")==0) {
      strcpy(name,"test");
    }
    else if (strcmp(argv[i],"-ntheta")==0) {
      sscanf(argv[++i],"%i",&nTheta);
      nTheta++;
    }
    else if (strcmp(argv[i],"-maxiter")==0) {
      sscanf(argv[++i],"%i",&maxiter);
    }
    else if (strcmp(argv[i],"-spherebox")==0) {
      sscanf(argv[++i],"%lf",&coat_ratio);
      shape=SPHEREBOX;
    }
    else if (strcmp(argv[i],"-coated")==0) {
      sscanf(argv[++i],"%lf",&coat_ratio);
      shape=COATED; strcpy(shapename,"coat");
      coat_x=0.0; coat_y=0.0; coat_z=0.0;
      {
	sscanf(argv[++i],"%lf",&coat_x);
	sscanf(argv[++i],"%lf",&coat_y);
	sscanf(argv[++i],"%lf",&coat_z);
      }
    }
    else if (strcmp(argv[i],"-ellipsoidal")==0) {
     sscanf(argv[++i],"%lf",&ellipsX);
     sscanf(argv[++i],"%lf",&ellipsY);
     sscanf(argv[++i],"%lf",&ellipsZ);
     ellipsY/=ellipsX;
     ellipsZ/=ellipsX;
     ellipsX=1;
     shape=ELLIPSOIDAL; strcpy(shapename,"ellipsoidal");
    }
    else if (strcmp(argv[i],"-sphere")==0) {
      shape=SPHERE; strcpy(shapename,"sphere");
    }
    else if (strcmp(argv[i],"-box")==0) {
      shape=BOX; strcpy(shapename,"box");
    }
    else if (strcmp(argv[i],"-disk")==0) {
      shape=DISK; strcpy(shapename, "disk");
      if (Narg!=1) {
        printz("Illegal number of arguments (%d) to -disk option\n",Narg);
        stop(1);
      }
      sscanf(argv[++i],"%lf",&diskratio);
    }
    else if (strcmp(argv[i],"-rbc")==0) {
      shape=RBC; strcpy(shapename, "rbc");
    }
    else if (strcmp(argv[i],"-aggregate")==0) {
      shape=AGGREGATE; strcpy(shapename,"aggregate");
      strcpy(aggregate_file,argv[++i]);
    }
    else if (strcmp(argv[i],"-read")==0) {
      shape=READ; strcpy(shapename,"read");
      if (Narg!=1) {
        printz("Illegal number of arguments (%d) to -read option\n",Narg);
        stop(1);
      }
      strcpy(aggregate_file,argv[++i]);
    }

    else if (strcmp(argv[i],"-dir")==0) {
      strcpy(directory,argv[++i]);
    }
    else if (strcmp(argv[i],"-xz")==0) {
      xzplane = true;
    }
    else if (strcmp(argv[i],"-Romb_inp")==0) {
      strcpy(Romb_parms,argv[++i]);
    }
    else if (strcmp(argv[i],"-all_dir")==0) {
      all_dir = true;
      /*yzplane = false;*/
      /*xzplane = false;*/
    }
    else if (strcmp(argv[i],"-store_all_dir")==0)
    {
      store_all_dir = true;
    }
    else if (strcmp(argv[i],"-Cext")==0) {
      calc_Cext = true;
    }    
    else if (strcmp(argv[i],"-Cabs")==0) {
      calc_Cabs = true;
    }
    else if (strcmp(argv[i],"-Csca")==0) {
      calc_Csca = true;
    }       
    else if (strcmp(argv[i],"-vec")==0) {
      calc_vec = true;
    }
    else if (strcmp(argv[i],"-asym")==0) {
      calc_asym = true;
      calc_vec = true;
      calc_Csca = true;
    }
    else if (strcmp(argv[i],"-Cpr_mat")==0) {
      calc_mat_force = true; 
    }
    else if (strcmp(argv[i],"-store_force")==0)
      {
	store_force = true;
      }
    else if (strcmp(argv[i],"-busy")==0) {
      char buffer[500];
      strcpy(procname,argv[++i]);
      sprintz(buffer,"touch /home/grimmink/runpar/busy.%s",procname);
      systemz(buffer);
    }
    else if (strcmp(argv[i],"-prisma")==0) {
      shape=PRISMA; strcpy(shapename,"prisma");
    }
    else if (strcmp(argv[i],"-line")==0) {
      shape=LINE; strcpy(shapename,"line");
    }
    else if (strcmp(argv[i],"-store")==0) {
      store_field=true;
    }
    else if (strcmp(argv[i],"-store_int_field")==0) {
      store_int_field=true;
    }
    else if (strcmp(argv[i],"-prognose")==0) {
      prognose=true;
    }
    else if (strcmp(argv[i],"-nosym")==0) {
      NoSymmetry=true;
    }
    else if (strcmp(argv[i],"-sym_enf")==0) {
      symmetry_enforced=true;
    }
    else if (strcmp(argv[i],"-prop")==0) {
      if (Narg!=3) {
        printz("Illegal number of arguments (%d) to -prop option\n",Narg);
        stop(1);
      }
      sscanf(argv[++i],"%lf",prop_0);
      sscanf(argv[++i],"%lf",&prop_0[1]);
      sscanf(argv[++i],"%lf",&prop_0[2]);     
      doubTmp=DotProd(prop_0,prop_0);
      if (doubTmp==0) {
        printz("Given propagation vector is null\n");
        stop(1);
      }
      doubTmp=1/sqrt(doubTmp);
      prop_0[0]*=doubTmp;
      prop_0[1]*=doubTmp;
      prop_0[2]*=doubTmp;
    }

    else if (strcmp(argv[i],"-beam")==0) {
      i++;
      if (strcmp(argv[i],"plane")==0) {beamtype=PLANE;}
      else {
	if (strcmp(argv[i],"lminus")==0) beamtype=LMINUS;
	else if (strcmp(argv[i],"buggy")==0) beamtype=BUGGY;
	else if (strcmp(argv[i],"davis1")==0) beamtype=DAVIS1;
	else if (strcmp(argv[i],"davis3")==0) beamtype=DAVIS3;
	else if (strcmp(argv[i],"barton1")==0) beamtype=BARTON1;
	else if (strcmp(argv[i],"barton3")==0) beamtype=BARTON3;
	else if (strcmp(argv[i],"barton5")==0) beamtype=BARTON5;
	else {
	  printz("Beam type '%s' is not supported\n",argv[i]);
          stop(1);
        }
        if (Narg!=5) {
          printz("Illegal number of arguments (%d) to '-beam %s' option\n",Narg-1,argv[i]);
          stop(1);
        }
	sscanf(argv[++i],"%lf",&beam_w0);
	sscanf(argv[++i],"%lf",&beam_x0);
	sscanf(argv[++i],"%lf",&beam_y0);
	sscanf(argv[++i],"%lf",&beam_z0);
      }
    }

    else if (strcmp(argv[i],"-pol")==0) {
      if (Narg>2) {
        printz("Illegal number of arguments (%d) to -pol option\n",Narg);
        stop(1);
      }
      i++;
      if (strcmp(argv[i],"cm")==0) PolRelation=CM;
      else if (strcmp(argv[i],"rrc")==0) PolRelation=RADCOR;
      else if (strcmp(argv[i],"ldr")==0) PolRelation=LDR;
      else if (strcmp(argv[i],"cldr")==0) PolRelation=CLDR;
      else if (strcmp(argv[i],"so")==0) PolRelation=SOrd;
      else {
        printz("Polarization Relation '%s' is not supported\n",argv[i]);
        stop(1);
      }
      if (Narg==2) {
        i++;
        if (strcmp(argv[i],"avgpol")==0) avg_inc_pol=true;
        else {
          printz("Unknown second argument '%s' to -pol option\n",argv[i]);
          stop(1);
        }
      }
    }

    else if (strcmp(argv[i],"-scat")==0) {
      if (Narg!=1) {
        printz("Illegal number of arguments (%d) to -scat option\n",Narg);
        stop(1);
      }
      i++;
      if (strcmp(argv[i],"dr")==0) ScatRelation=DRAINE;
      else if (strcmp(argv[i],"so")==0) ScatRelation=SOrd;
      else {
        printz("Scattering Quantities Relation '%s' is not supported\n",argv[i]);
        stop(1);
      }
    }

    else if (strcmp(argv[i],"-int")==0) {
      if (Narg!=1) {
        printz("Illegal number of arguments (%d) to -int option\n",Narg);
        stop(1);
      }
      i++;
      if (strcmp(argv[i],"poi")==0) IntRelation=POINT;
      else if (strcmp(argv[i],"so")==0) IntRelation=SOrd;
      else {
        printz("Interaction term prescription '%s' is not supported\n",argv[i]);
        stop(1);
      }
    }
    else if (strcmp(argv[i],"-orient")==0) {
      i++;
      if (strcmp(argv[i],"avg")==0) {
        if (Narg>2) {
          printz("Illegal number of arguments (%d) to '-orient avg' option\n",Narg-1);
          stop(1);
        }
        orient_avg=true;
        if (Narg==2) strcpy(avg_parms,argv[++i]);
      }
      else {
        if (Narg!=3) {
          printz("Illegal number of arguments (%d) to -orient option\n",Narg);
          stop(1);
        }
        sscanf(argv[i],"%lf",&alph_deg);
	sscanf(argv[++i],"%lf",&bet_deg);
	sscanf(argv[++i],"%lf",&gam_deg);
      }
    }
         
    else {
      printz("scatter: illegal option \"%s\"\n",argv[i]);
      printz("Usage: scatter [-lambda %%f] [-grid %%i [%%i %%i]] [-m %%f %%f %%f %%f]\n");
      printz("    [-dpl %%f] [-size %%f] [-test] [-ntheta %%i] [-maxiter %%i] [-dir %%s] \n");
      printz("    [-xz] [-Romb_inp %%s] [-sll_dir] [-store_all_dir] [-Cext] [-Cabs]\n");
      printz("    [-Csca] [-vec] [-asym] [-Cpr_mat] [-store_force] [-busy] [-store]\n");
      printz("    [-store_int_field] [-prognose] [-prop %%f %%f %%f] [-nosym] [-sym_enf] [-eps %%f]\n");
      printz("    [-cylinder|-lymphocyte1|-lymphocyte2|-leucocyte2|-rbc_rot|-sdisk_rot|\n");
      printz("     -spheroid_rot|-stick|-spherebox|-coated|-ellipsoidal|-sphere|-box|\n");
      printz("     -disk|-rbc|-aggregate|-prisma|-line[...]|-read filename] \n");
      printz("    [-beam plane|lminus|buggy|davis1|davis3|barton1|barton3|barton5 [...]]\n");
      printz("    [-pol cm|rrc|ldr|cldr|so] [-scat dr|so] [-int poi|so] [-orient avg|%%f %%f %%f]\n");
      printz("    \n");
      stop(1);
    }
    i++;
  } while(i<argc);
  
  if (prop_0[2]!=1 && orient_avg) {
    printz("-prop and '-orient avg' can not be used together\n");
    stop(1);
  }
  
  if (IntRelation==SOrd) reduced_FFT=false;
  else reduced_FFT=true;
    
  init_comm(&argc,&argv);    /* initialize communications */

  if (calc_Csca || calc_vec)
    {
      all_dir = true;
      /*yzplane = false;*/
    }

  if (store_all_dir)
    {
      all_dir = true;
      /*yzplane = false;*/
    }

  /*determine two incident polarizations. Equivalent to rotation of X,Y,Z basis by 
  angles Theta and Phi from (0,0,1) to given propagation vector*/
  if (abs(prop_0[2])==1) {
    incPolX_0[0]=prop_0[2];
    incPolY_0[1]=1;
    incPolX_0[1]=incPolX_0[2]=incPolY_0[0]=incPolY_0[2]=0.0;
  }
  else {
    temp=sqrt(1-prop_0[2]*prop_0[2]);
    incPolX_0[0]=prop_0[0]*prop_0[2]/temp;
    incPolX_0[1]=prop_0[1]*prop_0[2]/temp;
    incPolX_0[2]=-temp;
    incPolY_0[0]=-prop_0[1]/temp;
    incPolY_0[1]=prop_0[0]/temp;
    incPolY_0[2]=0.0;
  }
 
  if (orient_avg) {
    read_avg_parms(avg_parms);
    NoSymmetry=true;
    avg_inc_pol=true;
  }
  else {
    init_rotation();                  /* initialize rotation stuff */
    if (prop[2]!=1) NoSymmetry=true;
  }
  /* get number of dipoles */
  init_shape(shape);          /* initialize symmetries and box's */
  nDip=boxX*boxY*boxZ;

  if (dpl==UNDEF) {
    dpl=LAMBDA*boxX/sizeX;
 /*   if (shape==SPHERE) {*/ /* correct dpl to give equal volume sphere */
 /*     dpl=LAMBDA*pow((6/PI)*nDip,.333333333333)/sizeX;     is not that simple!!
    } */ 
  }
  gridspace=LAMBDA/dpl;
  
  if (nTheta==UNDEF) {
    if (boxX*boxY*boxZ<1000) nTheta=91;
    else if (boxX*boxY*boxZ<10000) nTheta=181;
    else if (boxX*boxY*boxZ<100000) nTheta=361;
    else nTheta=721;
  }
    
  if (prognose==false && directory[0]==0) { 
    if ((Nexpfile=fopen("ExpCount","r"))!=NULL) {
      fscanf(Nexpfile,"%i",&Nexp);
      fclose(Nexpfile);
    }
    else Nexp=0;
    
    synchronize();
    
    if (ringid==0) {
      if ((Nexpfile=fopen("ExpCount","w"))==NULL) 
        LogError(EC_ERROR,ONE,POSIT,"Failed to open file 'ExpCount'");
      fprintz(Nexpfile,"%i",Nexp+1);
      fclose(Nexpfile);
    }
    sprintf(directory,"%s%03i_%sg%id%.1f",name,Nexp,shapename,boxX,dpl);
  }
  else if (directory[0]==0 || prognose==true) sprintf(directory,"/tmp");
  strcpy(sbuffer,"mkdir ");                /* make new directory */
  strcat(sbuffer,directory);
  
  if (prognose==false) systemz(sbuffer);
  
  strcpy(logname,directory);              /* make logfile */
  strcat(logname,"/log");
 
  strcpy(sbuffer,"echo 'The program was run on:'$HOST >>");
  strcat(sbuffer,logname);
  systemz(sbuffer);
  synchronize();

  if (ringid==0) if ((logfile=fopen(logname,"aA"))==NULL) 
    LogError(EC_ERROR,ONE,POSIT,"Failed to open file 'logfile'"); 
  
  printz("all data is saved in '%s'\n",directory);
  printz("lambda:%g      m:%g%+gi     Dipoles/lambda:%f\n",LAMBDA,ref_index[0].r,ref_index[0].i,dpl);
  printz("epsilon square:%g\n",eps0);
  printz("Particle has %i dipoles\n",nDip);
  
  fprintz(logfile,"command: '");
  for(i=0;i<argc;i++) fprintz(logfile,"%s ",argv[i]);
  fprintz(logfile,"'\n");

  fprintz(logfile,"lambda:%g\n",LAMBDA);
  fprintz(logfile,"shape: '%s', size:%ix%ix%i\n",shapename,boxX,boxY,boxZ);
  fprintz(logfile,"m0:%g%+gi  m1:%g%+gi\n",ref_index[0].r,ref_index[0].i,ref_index[1].r,ref_index[1].i);
  
  fprintz(logfile,"\nIncident propagation vector:(%g,%g,%g)\n",prop_0[0],prop_0[1],prop_0[2]);
  fprintz(logfile,"Incident polarization Y(par):(%g,%g,%g)\n",incPolY_0[0],incPolY_0[1],incPolY_0[2]);
  fprintz(logfile,"Incident polarization X(per):(%g,%g,%g)\n\n",incPolX_0[0],incPolX_0[1],incPolX_0[2]);

  if (orient_avg) fprintz(logfile,"Particle orientation - averaged\n%s\n",avg_string);
  else fprintz(logfile,"Particle orientation (deg): alpha=%g, beta=%g, gamma=%g\n\n",alph_deg,bet_deg,gam_deg);

  if (alph_deg!=0 || bet_deg!=0 || gam_deg!=0) {
    fprintz(logfile,"After transformation to particle reference frame:\n");
    fprintz(logfile,"New incident propagation vector:(%g,%g,%g)\n",prop[0],prop[1],prop[2]);
    fprintz(logfile,"New incident polarization Y(par):(%g,%g,%g)\n",incPolY[0],incPolY[1],incPolY[2]);
    fprintz(logfile,"New incident polarization X(per):(%g,%g,%g)\n\n",incPolX[0],incPolX[1],incPolX[2]);
  }

  fprintz(logfile,"Dipoles/lambda:%f\n",dpl);
  fprintz(logfile,"epsilon square:%g\n",eps0);
  
  if (PolRelation==CM) fprintz(logfile,"Polarization relation: 'Clausius-Mossotti'\n");
  else if (PolRelation==RADCOR) fprintz(logfile,"Polarization relation: 'Radiative Reaction Correction'\n");
  else if (PolRelation==LDR) {
    fprintz(logfile,"Polarization relation: 'Lattice Dispersion Relation'");
    if (avg_inc_pol) fprintz(logfile," (averaged over incident polarization)");
    fprintz(logfile,"\n");
  }
  else if (PolRelation==CLDR) fprintz(logfile,"Polarization relation: 'Corrected Lattice Dispersion Relation'\n");
  else if (PolRelation==SOrd) fprintz(logfile,"Polarization relation: 'Second Order'\n");
  
  if (ScatRelation==DRAINE) fprintz(logfile,"Scattering quantities formulae: 'by Draine'\n");
  else if (ScatRelation==SOrd) fprintz(logfile,"Scattering quantities formulae: 'Second Order'\n");

  if (IntRelation==POINT) fprintz(logfile,"Interaction term prescription: 'as Point dipoles'\n");
  else if (IntRelation==SOrd) fprintz(logfile,"Interaction term prescription: 'Second Order'\n");
      
  if (symmetry_enforced) fprintz(logfile,"Symmetry is enforced by user (warning!)\n");
  else if (NoSymmetry) fprintz(logfile,"No symmetries are used\n",eps0);

  /* initialize times and counters */
  TotalIter=TotalEval=TotalEFieldPlane=0;
  Timing_EField=Timing_FileIO=Timing_IntField=Timing_ScatQuan=0;
  Timing_Integration=0;
  
  Calculator(); /* Main calculation part */
  
  Timing_TotalTime = clock() - tstart_main;

  time(&end);

  fprintz(logfile,"\n\n"
          "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
#if defined(MPI)
          "This was a problem with %d dipoles\n"\
          "ran on %d processors with MPI\n"
#else
          "This was a problem with %d dipoles\n"\
          "ran sequentially\n"
#endif
          "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
          nDip, nprocs);
  fprintz(logfile,"\n ---- TIMING RESULTS ----\n\n");
  
  if (orient_avg) fprintz(logfile,"Total number of single particle evaluations: %d\n",TotalEval);
  fprintz(logfile,
          "Total number of iterations: %d\n"\
          "Total planes of E field calculation (each %d points): %d\n\n"\
          "total time:            %4.4f\n"\
          "Wall time:             %.1f\n"\
	  "Initialization time:   %4.4f\n"\
	  "  FFT setup:             %4.4f\n"\
	  "  make particle:         %4.4f\n"\
          "Internal fields:       %4.4f\n"\
          "  one solution:          %4.4f\n"\
          "    init solver:           %4.4f\n"\
          "    one iteration:         %4.4f\n"\
          "      calculation:           %4.4f\n"\
          "      communication:         %4.4f\n"\
          "E field calculation:   %4.4f\n"\
          "  one plane:             %4.4f\n"\
          "    calculation:           %4.4f\n"\
          "    communication:         %4.4f\n"\
          "  one alldir:            %4.4f\n"\
          "    calculation:           %4.4f\n"\
          "    communication:         %4.4f\n"\
          "Other scat.quantities: %4.4f\n"\
          "file io:               %4.4f\n"\
          "Integration:           %4.4f\n",
          TotalIter,nTheta,TotalEFieldPlane,
          TO_SEC(Timing_TotalTime),difftime(end,start),
          TO_SEC(Timing_Init),TO_SEC(Timing_FFT_Init),TO_SEC(Timing_Particle),
          TO_SEC(Timing_IntField),TO_SEC(Timing_IntFieldOne),TO_SEC(Timing_InitIter),
          TO_SEC(Timing_OneIter),TO_SEC(Timing_OneIterCalc),TO_SEC(Timing_OneIterComm),
          TO_SEC(Timing_EField),
          TO_SEC(Timing_EFieldPlane),TO_SEC(Timing_calc_EField),TO_SEC(Timing_comm_EField),
          TO_SEC(Timing_EField_ad),TO_SEC(Timing_calc_EField_ad),TO_SEC(Timing_comm_EField_ad),
          TO_SEC(Timing_ScatQuan),TO_SEC(Timing_FileIO),TO_SEC(Timing_Integration));

  fclosez(logfile);
  stop(0);
}
