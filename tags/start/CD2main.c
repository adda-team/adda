/* FILE: CD2main.c
 * AUTH: Alfons Hoekstra (alfons@fwi.uva.nl)
 * VERS: express version 1.2
 * HIST: 1.1 creation
 */

/* sequential version, Michel Grimminck Jan 1995 */


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "cmplx.h"
#include "types.h"
#include "const.h"
#include <math.h>
#include <limits.h>
#include "comm.h"
#include "debug.h"

#ifndef PowerPVM
#define CLOCK_TICK 1000000   /* should be read from some header file */
#endif

#define TO_SEC(p) (p / (double) CLOCK_TICK)

/* global variables */
int        MyProcId=0,		/* ID of this processor */
           nProcs=1,		/* number of booted processors */
           TheRing;		/* ID of the Ring topology */
RingData_t *TheRing_data;	/* pointer to ringdata struct */

extern int __fake_nnodes;
extern int nDip;

int shape;                       /* particle shape definition */
int boxX,boxY,boxZ;              /* grid of calculation box */
int symX,symY,symZ;              /* symmetries in x y and z direction of
				    particle */
int symR;                        /* symmetry if 90 deg. rotation about z axes
				  */
int symD;                        /* symmetry in D-matrix */
int NoSymmetry;                  /* do not use particle symmetries */
int startclock;
int PolRelation;                 /* RRC,LDR or CM polerisation relation */
int fsc;                         /* use fase shift correction */
int prognose;                    /* make a prognose about needed ram,virtmem */
int store_field;                 /* save internal fields to pbm file */
int store_int_field;		/* save full internal fields to file */
double eps0;                     /* stop criterium */
char directory[200];             /* directory to save data in */
FILE *logfile;
double gridspaceX,gridspaceY,gridspaceZ;  /* inter-dipole distance */
double dplX,dplY,dplZ,sizeX,sizeY,sizeZ;
int profile[10];
int analyse=false;               /* automatic analyse */
char mark[500];                  /* mark the results with this string */
int maxiter;                     /* the maximal number of iterations */
int jagged;                      /* particle is jagged */
int beamtype;                    /* barton,davis,lminus */
double beam_w0,beam_x0,beam_y0,beam_z0=0.0; /* beam properties (meters) */
doublecomplex ref_index[10];          /* a set of refraction indexes */
doublecomplex cc[10];                 /* the corresponding coupleconstants */
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
int calc_Csca_diff = false;      /* Calculate the formal difference between
				    the direct and integrated scattering
				    cross section */
int calc_vec = false;            /* Calculate the unnormalized 
				    asymmetry-parameter */
int calc_asym = false;           /* Calculate the asymmetry-parameter */
int calc_mat_force = false;      /* Calculate the scattering force 
				    by matrix-evaluation */
int store_force = false;         /* Write radiation pressure per dipole
				    to file */
char Romb_parms[200];
char aggregate_file[200];

/*
 * start adress of vectors, used by both the calculator and the router. 
 * The vectors are initialized by the Calculator
 */

dcomplex *buffer, *x, *r, *p, *Epar, *Eper, *Einc, *Egrid;

/* Some timing variables */
unsigned long Timing_TotalTime,Timing_FFT_Init,
  Timing_OneIter, Timing_OneIterCalc, Timing_OneIterComm, 
  Timing_EField, Timing_calc_EField, Timing_comm_EField, 
  Timing_FileIO = 0, Timing_Integration;

void
main (int argc, char *argv[])
{
  char          *Fname = __FILE__;
  struct nodenv nodedata;
  unsigned int  tstart, tstop;	/* timers */
  char          *pname,sbuffer[500],logname[100],name[100],shapename[20];
  int		i,Nexp,Narg;
  FILE          *Nexpfile;
  time_t        start,end;
  char          procname[50];
  
  extern double LAMBDA;	/* defined in calculator.c */
  extern double DIPS_PER_LAMBDA;
  extern int    nDip,nTheta;
  
  printf("hello!!\n");
  init_comm(&argc,&argv);
  
  time(&start);
  
  for(i=0;i<10;i++) profile[i]=0;
  for(i=0;i<10;i++) {ref_index[i].r=1.0; ref_index[i].i=0.0;}
  /* defaults */
  directory[0]=0;
  LAMBDA=2*PI;
  ref_index[0].r=1.05;
  ref_index[0].i=0.0;
  boxX=boxY=boxZ=8;
  sizeX=sizeY=sizeZ=UNDEF;
  DIPS_PER_LAMBDA=10;
  dplX=dplY=dplZ=10;
  strcpy(name,"run");
  nTheta=-1; /* undecided */
  eps0=/*1.0e-12*/1.0e-10;
  shape=ELLIPS; strcpy(shapename,"ellips");
  store_field=false;
  store_int_field=false;
  PolRelation=LDR;
  NoSymmetry=false;
  fsc=false;
  prognose=false;
  analyse=false;
  maxiter=INT_MAX;
  jagged=1;
  procname[0]=0;
  beamtype=PLANE;
  strcpy(Romb_parms,"Romb5.dat");
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
      sscanf(argv[++i],"%i",&boxX);
      boxY=boxZ=boxX;
    }
    else if (strcmp(argv[i],"-m")==0) {
      int j;
      
      if (Narg<1) {
	printz("Illegal number of arguments to -m option\n");
	stop(1);
      }
      Narg=4;
      Nmat=(Narg+1)/2;
      /*#ifdef PARALLEL*/
      Nmat++;
      /*#endif*/
      if (Nmat>10) {printz("Too many materials\n"); stop(1);}
      sscanf(argv[++i],"%lf",&ref_index[0].r);
      sscanf(argv[++i],"%lf",&ref_index[0].i);
      
      sscanf(argv[++i],"%lf",&ref_index[1].r);
      sscanf(argv[++i],"%lf",&ref_index[1].i);
    }
    else if (strcmp(argv[i],"-dpl")==0) {
      if (Narg!=1) {
        printz("Illegal number of arguments to -dpl option\n");
        stop(1);
      }
      sscanf(argv[++i],"%lf",&dplX);
      DIPS_PER_LAMBDA=dplY=dplZ=dplX;
    }
    else if (strcmp(argv[i],"-cylinder")==0) {
      if (Narg!=2) {
        printz("Illegal number of arguments to -particle option\n");
        stop(1);
      }
      shape=CYLINDER; strcpy(shapename,"cylinder");
      sscanf(argv[++i],"%lf",&aspect_ratio);
      sscanf(argv[++i],"%lf",&append_ratio);
    }
    else if (strcmp(argv[i],"-lymphocyte1")==0) {
      if (Narg!=5) {
        printz("Illegal number of arguments to -particle option\n");
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
        printz("Illegal number of arguments to -particle option\n");
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
        printz("Illegal number of arguments to -particle option\n");
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
        printz("Illegal number of arguments to -particle option\n");
        stop(1);
      }
      shape=RBC_ROT; strcpy(shapename,"rbc_rot");
      sscanf(argv[++i],"%lf",&aspect_r);
      sscanf(argv[++i],"%lf",&betaY);
      sscanf(argv[++i],"%lf",&betaZ);
    }
    else if (strcmp(argv[i],"-sdisk_rot")==0) {
      if (Narg!=3) {
        printz("Illegal number of arguments to -particle option\n");
        stop(1);
      }
      shape=SDISK_ROT; strcpy(shapename,"sdisk_rot");
      sscanf(argv[++i],"%lf",&aspect_r);
      sscanf(argv[++i],"%lf",&betaY);
      sscanf(argv[++i],"%lf",&betaZ);
    }  
    else if (strcmp(argv[i],"-spheroid_rot")==0) {
      if (Narg!=3) {
        printz("Illegal number of arguments to -particle option\n");
        stop(1);
      }
      shape=SPHEROID_ROT; strcpy(shapename,"spheroid_rot");
      sscanf(argv[++i],"%lf",&aspect_r);
      sscanf(argv[++i],"%lf",&betaY);
      sscanf(argv[++i],"%lf",&betaZ);
    }   
    else if (strcmp(argv[i],"-stick")==0) {
      if (Narg!=4) {
        printz("Illegal number of arguments to -particle option\n");
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
      sizeY=sizeZ=sizeX;
      dplX=dplY=dplZ=DIPS_PER_LAMBDA=UNDEF;
    }
    else if (strcmp(argv[i],"-test")==0) {
      strcpy(name,"test");
    }
    else if (strcmp(argv[i],"-ntheta")==0) {
      sscanf(argv[++i],"%i",&nTheta);
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
        printz("Illegal number of arguments to -disk option\n");
        stop(1);
      }
      sscanf(argv[++i],"%lf",&diskratio);
    }
    else if (strcmp(argv[i],"-symdisk")==0) {
      shape=SDISK; strcpy(shapename, "symdisk");
      if (Narg!=1) {
        printz("Illegal number of arguments to -symdisk option\n");
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
    else if (strcmp(argv[i],"-mark")==0) {
      strcpy(mark,argv[++i]);
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
    else if (strcmp(argv[i],"-ellips")==0) {
      shape=ELLIPS; strcpy(shapename,"ellips");
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
    else if (strcmp(argv[i],"-fsc")==0) {
      fsc=true;
    }
    else if (strcmp(argv[i],"-prognose")==0) {
      prognose=true;
    }
    else if (strcmp(argv[i],"-cm")==0) {
      PolRelation=CM;
    }
    else if (strcmp(argv[i],"-ldr")==0) {
      PolRelation=LDR;
    }
    else if (strcmp(argv[i],"-rrc")==0) {
      PolRelation=RADCOR;
    }
    else if (strcmp(argv[i],"-nosym")==0) {
      NoSymmetry=true;
    }
    else if (strcmp(argv[i],"-beam")==0) {
      i++;
      if (strcmp(argv[i],"plane")==0) {beamtype=PLANE;}
      else {
	if (strcmp(argv[i],"lminus")==0) beamtype=LMINUS;
	if (strcmp(argv[i],"buggy")==0) beamtype=BUGGY;
	if (strcmp(argv[i],"davis1")==0) beamtype=DAVIS1;
	if (strcmp(argv[i],"davis3")==0) beamtype=DAVIS3;
	if (strcmp(argv[i],"barton1")==0) beamtype=BARTON1;
	if (strcmp(argv[i],"barton3")==0) beamtype=BARTON3;
	if (strcmp(argv[i],"barton5")==0) beamtype=BARTON5;
	sscanf(argv[++i],"%lf",&beam_w0);
	sscanf(argv[++i],"%lf",&beam_x0);
	sscanf(argv[++i],"%lf",&beam_y0);
	sscanf(argv[++i],"%lf",&beam_z0);
      }
    }
    else if (strcmp(argv[i],"-analyse")==0) {
      analyse=true;
    }
    else {
      printf("scatter: illegal option \"%s\"\n",argv[i++]);
      printz("Usage: CD2main.x [-lambda %%f] [-grid %%i [%%i %%i]] [-m %%f [%%f]]\n");
      printz("    [-dpl %%f[%%f %%f]] [-test] [-ntheta %%i] [-eps %%f] [-nosym]\n");
      printz("    [-sphere|-prisma|-box|-ellips|-coated %%f] [-store] [-cm|-ldr|-rrc]\n");
      printz("    [-fsc] [-prognose] [-analyse] [-size %%f [%%f %%f]]\n");
      printz("    [-mark %%s]\n");
      stop(1);
    }
    i++;
  } while(i<argc);
  
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

  /* get number of dipoles */
  gridspaceX=gridspaceY=gridspaceZ=1.0;
  nDip = make_particle (NULL,false,shape,jagged);
  printz("Particle has %i dipoles\n",nDip);
  
  if (DIPS_PER_LAMBDA==UNDEF) {
    dplX=boxX/sizeX;
    dplY=boxY/sizeY;
    dplZ=boxZ/sizeZ;
    DIPS_PER_LAMBDA=dplX;
    if (shape==SPHERE) { /* correct dpl to give equal volume sphere */
      double r;
      r=2.0*pow((0.75/PI)*nDip/dplX/dplY/dplZ,.333333333333);
      r/=sizeX;
      dplX*=r; dplY*=r; dplZ*=r;
    }
    DIPS_PER_LAMBDA=dplX;
  }
  gridspaceX=LAMBDA/dplX;
  gridspaceY=LAMBDA/dplY;
  gridspaceZ=LAMBDA/dplZ;
  if (nTheta==-1) {
    if (boxX*boxY*boxZ<1000) nTheta=90;
    else if (boxX*boxY*boxZ<10000) nTheta=180;
    else if (boxX*boxY*boxZ<100000) nTheta=360;
    else nTheta=720;
  }
  
  if (prognose==false && directory[0]==0) { 
    Nexpfile=fopen("ExpCount","r");
    if (Nexpfile!=NULL) {
      fscanf(Nexpfile,"%i",&Nexp);
      fclose(Nexpfile);
    }
    else Nexp=0;
    
    synchronize();
    
    if (me==0) {
      Nexpfile=fopen("ExpCount","w");
      if (Nexpfile==NULL) {printz("failed to open file 'ExpCount'\n"); stop(1);}
      fprintz(Nexpfile,"%i",Nexp+1);
      fclose(Nexpfile);
      
    }
    sprintf(directory,"%s%i_%sg%id%.1f",name,Nexp,shapename,boxX,DIPS_PER_LAMBDA);
  }
  else {
    if (directory[0]==0 || prognose==true) sprintf(directory,"/tmp");
  }
  strcpy(sbuffer,"mkdir ");                /* make new directory */
  strcat(sbuffer,directory);
  
  if (prognose==false) systemz(sbuffer);
  
  strcpy(logname,directory);              /* make logfile */
  strcat(logname,"/log");
 
  strcpy(sbuffer,"echo 'The program was run on:'$HOST >>");
  strcat(sbuffer,logname);
  systemz(sbuffer);
  synchronize();

  if (me==0) {
    logfile = fopen (logname, "aA");
    if (logfile == NULL) {printz("failed to open file 'logfile'\n"); stop(1);}
  }
  printz("all data is saved in '%s'\n",directory);
  printz("lambda:%g      m:%.3f+%.3fi     Dipoles/lambda:%f\n",LAMBDA,ref_index[0].r,ref_index[0].i,DIPS_PER_LAMBDA);
  printz("epsilon square:%g\n",eps0);
  printz("Particle has %i dipoles\n",nDip);
  
  tstart = extime ( );
  
  fprintz(logfile,"command: '");
  for(i=0;i<argc;i++) fprintz(logfile,"%s ",argv[i]);
  fprintz(logfile,"'\n");
  if (sizeof(REAL)==sizeof(float))
    fprintz(logfile,"Running in single precision\n");
  if (sizeof(REAL)==sizeof(double))
    fprintz(logfile,"Running in double precision\n");

  fprintz(logfile,"lambda:%g\n",LAMBDA);
  fprintz(logfile,"shape: '%s', size:%ix%ix%i\n",shapename,boxX,boxY,boxZ);
  fprintz(logfile,"m0:%f+%fi  m1:%f+%fi\n",ref_index[0].r,ref_index[0].i,ref_index[1].r,ref_index[1].i);
  fprintz(logfile,"Dipoles/lambda:%f\n",DIPS_PER_LAMBDA);
  fprintz(logfile,"epsilon square:%g\n",eps0);
  if (PolRelation==CM) fprintz(logfile,"Polarization relation: 'Clausius-Mossotti'\n");
  else if (PolRelation==LDR) fprintz(logfile,"Polarization relation: 'lattice dispersion'\n");
  else if (PolRelation==RADCOR) fprintz(logfile,"Polarization relation: 'radiative reaction correction'\n");
  if (fsc==true) fprintz(logfile,"Fase Shift Correction used\n");
  startclock=clock();
  
  nProcs = nprocs;
  MyProcId = ringid;  

  Calculator();
  tstop = extime ( );
  
  Timing_TotalTime = tstop - tstart;

  time(&end);

  fprintz(logfile,"\n\n"
          "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
#if defined(PVM)
          "This was a problem with %d dipoles\n"\
          "ran on %d processors with PVM\n"
#elif defined(MPI)
          "This was a problem with %d dipoles\n"\
          "ran on %d processors with MPI\n"
#else
          "This was a problem with %d dipoles\n"\
          "ran sequentially\n"
#endif
          "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
          nDip, nprocs);

  fprintz(logfile,"\n ---- TIMING RESULTS ----\n\n");
  fprintz(logfile,"total time:          %4.6f\n"\
          "Wall time:           %.1f\n"\
	  "FFT setup:           %4.4f\n"\
          "one iteration:       %4.4f\n"\
          "    calculation:     %4.4f\n"\
          "    communication:   %4.4f\n"\
          "E field calculation: %4.4f\n"\
          "    calculation:     %4.4f\n"\
          "    communication:   %4.4f\n"\
          "file io:             %4.4f\n"\
          "Integration:         %4.4f\n",
          TO_SEC(Timing_TotalTime), difftime(end,start),
	  TO_SEC(Timing_FFT_Init),
          TO_SEC(Timing_OneIter),
          TO_SEC(Timing_OneIterCalc),
          TO_SEC(Timing_OneIterComm),
          TO_SEC(Timing_EField),
          TO_SEC(Timing_calc_EField), /* is only used in crosssec.c */
          TO_SEC(Timing_comm_EField), /* idem */
          TO_SEC(Timing_FileIO),
          TO_SEC(Timing_Integration));/* is only used in crosssec.c */

  fclosez(logfile);
  stop(0);
}
