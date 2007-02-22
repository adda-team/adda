/* FILE: io.c
 * AUTH: Maxim Yurkin
 * DESCR: io routines
 *
 */
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#ifdef _WIN32  /* Windows */
# include <windows.h>
#endif

#include "io.h"
#include "const.h"
#include "comm.h"
#include "vars.h"
#include "crosssec.h"

/* definitions for file locking */
#ifdef USE_LOCK
# ifdef _WIN32  /* Windows */
#  include <windows.h>
#  define FILEHANDLE HANDLE
# else          /* UNIX - may depend on which exactly */
#  include <unistd.h>
#  include <fcntl.h>
#  define FILEHANDLE int
# endif
# define LOCK_WAIT 1              /* in seconds */
# define MAX_LOCK_WAIT_CYCLES 60
#else
# define FILEHANDLE int
#endif

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in crosssec.c */
extern char avg_string[];
/* defined and initialized in make_particle.c */
extern int volcor_used;
extern char sh_form_str[];

/* used in CalculateE.c */
int store_int_field; /* save full internal fields to text file */
int store_scat_grid; /* Store the scattered field for grid of angles */
int calc_Cext;       /* Calculate the extinction cross-section - allways do */
int calc_Cabs;       /* Calculate the absorption cross-section - allways do */
int calc_Csca;       /* Calculate the scattering cross-section by integration */
int calc_vec;        /* Calculate the unnormalized asymmetry-parameter */
int calc_asym;       /* Calculate the asymmetry-parameter */
int calc_mat_force;  /* Calculate the scattering force by matrix-evaluation */
int store_force;     /* Write radiation pressure per dipole to file */
int phi_int_type;    /* type of phi integration (each bit determines
                          whether to calculate with different multipliers) */
/* used in calculator.c */
int PolRelation;           /* type of polarization relation */
int avg_inc_pol;           /* whether to average CC over incident polarization */
char alldir_parms[200];    /* name of file with alldir parameters */
char scat_grid_parms[200]; /* name of file with parameters of scattering grid */
/* used in crosssec.c */
double prop_0[3];  /* initial incident direction (in laboratory reference frame) */
double incPolX_0[3],incPolY_0[3]; /* initial incident polarizations (in lab RF)*/
int ScatRelation;                 /* type of formulae for scattering quantities */
/* used in GenerateB.c */
int beamtype;                           /* type of incident beam */
double beam_w0,beam_x0,beam_y0,beam_z0; /* beam properties (microns) */
/* used in iterative.c */
double eps;                     /* relative error to reach */
/* used in make_particle.c */
int shape;                      /* particle shape definition */
int sh_Npars;                   /* number of shape parameters */
double sh_pars[MAX_N_SH_PARMS]; /* storage for shape parameters */
int NoSymmetry;            /* do not use particle symmetries */
int symmetry_enforced;     /* enforce use of all symmetries */
double sizeX;              /* size of particle along x-axis */
double dpl;                /* number of dipoles per lambda (wavelength) */
double lambda;             /* incident wavelength (in vacuum) */
int jagged;                /* size of big dipoles, used to construct a particle */
char aggregate_file[200];  /* name of aggregate file */
char save_geom_fname[200]; /* geometry file name to save dipole configuration */
char shapename[20];        /* name of the shape used */
int volcor;                /* whether to use volume correction */
int save_geom;             /* whether to save dipole configuration in .geom file */

/* LOCAL VARIABLES */

char run_name[10];      /* first part of the dir name ('run' or 'test') */
char logname[200];      /* name of logfile */
char avg_parms[200];    /* name of file with orientation averaging parameters */

/*============================================================*/

FILEHANDLE CreateLockFile(char *fname)
   /* create locks file
      works only if USE_LOCK is enabled */
{
#ifdef USE_LOCK
  FILEHANDLE fd;
  int i;

# ifdef _WIN32                /* Windows */
  i=0;
  while ((fd=CreateFile(fname,GENERIC_WRITE,FILE_SHARE_WRITE,NULL,CREATE_NEW,
                        FILE_ATTRIBUTE_NORMAL,NULL))==INVALID_HANDLE_VALUE) {
    Sleep(LOCK_WAIT*1000);
    if (i++ == MAX_LOCK_WAIT_CYCLES)
      LogError(EC_ERROR,ONE,POSIT,"Lock file %s permanently exists",fname);
  }
# else /* UNIX */
  struct flock lock;

  /* open file exclusively */
  i=0;
  while ((fd=open(fname,O_WRONLY | O_CREAT | O_EXCL,0666))==-1) {
    sleep(LOCK_WAIT);
    if (i++ == MAX_LOCK_WAIT_CYCLES)
      LogError(EC_ERROR,ONE,POSIT,"Lock file %s permanently exists",fname);
  }
  /* specify lock - file is additionally locked to work robustly over NFS */
  lock.l_type=F_WRLCK;
  lock.l_whence=SEEK_SET;
  lock.l_start=0;
  lock.l_len=0;
  /* obtain lock */
  if (fcntl(fd,F_SETLKW,&lock)==-1)
    LogError(EC_ERROR,ONE,POSIT,"Obtaining file lock failed");
# endif
  /* return file handle */
  return fd;
#else
  return 0;
#endif
}

/*============================================================*/

void RemoveLockFile(FILEHANDLE fd,char *fname)
   /* closes and remove lock file
      works only if USE_LOCK is enabled */
{
#ifdef USE_LOCK
# ifdef _WIN32                /* Windows */
  /* close file and remove it */
  CloseHandle(fd);
  DeleteFile(fname);
# else /* UNIX */
  char bufstr[100];
  /* close file; all locks are automatically released. Then remove it */
  close(fd);
  sprintf(bufstr,"rm %s",fname);
  system(bufstr);
# endif
#endif
}

/*============================================================*/

void LogError(int code, int who, char *fname, int lineN, char *fmt, ... )
   /* performs output of error specified by code at fname:lineN
    * fmt + args (...) specify error message
    * who specifies whether 1 (ringid=ROOT) or all processors should produce output
    * if code is EC_ERROR program aborts after output.
    * We use sprintf a couple of times, because we want each node to
    * generate an atomic message, not a couple of messages after
    * each other, since other nodes may then interfere with our output!
    */
{
  va_list args;
  char line[255];
  char id_str[30];
  extern char logname[];
  extern FILE *logfile;

  if (who==ALL || ringid==ROOT) {
    va_start(args, fmt);

    strcpy(id_str,"");
#ifdef PARALLEL
    if (who==ALL) sprintf(id_str," - ringID=%d",ringid);
#endif
    if (code==EC_ERROR) strcpy(line,"ERROR:");
    else if (code==EC_WARN) strcpy(line,"WARNING:");
    else sprintf(line,"Error code=%d:",code);
    sprintf(line+strlen(line)," (%s:%d%s) ", fname, lineN, id_str);
    vsprintf(line+strlen(line), fmt, args);
    strcat(line,"\n");
    fprintf(stderr,line);
    fflush(stderr);
    va_end(args);
  }
  if (code==EC_ERROR) {
    /* duplicate error message in logfile */
    if (logname[0]!=0) {  /* otherwise can't do anything */
      if (ringid==ROOT) {
        /* logfile==NULL with logname!=0 may only occur when error is
           in opening of the logfile itself */
        if (logfile!=NULL) {
          fprintf(logfile,line);
          fclose(logfile);
        }
      }
      else if (who==ALL) {
        logfile=fopen(logname,"a");
        fprintf(logfile,line);
        fclose(logfile);
      }
    }
    Stop(1);
  }
}

/*============================================================*/

void PrintBoth(FILE *file,char *fmt, ... )
   /* print anything both to file and to stdout */
{
  va_list args;
  char line[255];

  va_start(args,fmt);
  vsprintf(line,fmt,args);
  fprintf(file,line);
  printf(line);
  va_end(args);
}


/*============================================================*/

INLINE void IllNarg(char *string,int Nis,int Nreq)
      /* check if N==Nreq, otherwise produces error message */
{
  if (Nis!=Nreq) {
    printz("Illegal number of arguments (%d) to %s option (%d expected)\n",Nis,string,Nreq);
    Stop(1);
  }
}

/*============================================================*/

int TimeField(char c)
   /* analize one time multiplier */
{
  if (c=='d' || c=='D') return 86400;
  else if (c=='h' || c=='H') return 3600;
  else if (c=='m' || c=='M') return 60;
  else if (c=='s' || c=='S' || c==0) return 1;
  else {
    printz("Illegal time format specifier (%c)\n",c);
    Stop(1);
  }
  /* never reached */
  return 0;
}

/*============================================================*/

int ScanTime(char *str)
   /* scans time in seconds from a string "%d[d,D[%d]][h,H[%d]][m,M[%d]][s,S] */
{
#define TIME_N_TYPES 4   /* not that easy to change */
  int time,t[TIME_N_TYPES],n,i;
  char c[TIME_N_TYPES];

  for (i=0;i<TIME_N_TYPES;i++) c[i]=0;
  n=sscanf(str,"%d%c%d%c%d%c%d%c",t,c,t+1,c+1,t+2,c+2,t+3,c+3);
  if (n<1) {
    printz("Wrong time format\n",c);
    Stop(1);
  }
  time=0;
  i=0;
  while (n>0) {
    time+=t[i]*TimeField(c[i]);
    n-=2;
    i++;
  }
  return time;
#undef TIME_N_TYPES
}

/*============================================================*/

void PrintTime(char *s,time_t *time_ptr)
{
   struct tm *t;

   t=gmtime(time_ptr);
   s[0]=0; /* initialize string */
   if (t->tm_yday>0) sprintf(s,"%dd ",t->tm_yday);
   if (t->tm_hour>0) sprintf(s+strlen(s),"%dh ",t->tm_hour);
   if (t->tm_min>0) sprintf(s+strlen(s),"%dm ",t->tm_min);
   if (t->tm_sec>0) sprintf(s+strlen(s),"%ds ",t->tm_sec);
}

/*============================================================*/

void InitVariables(void)
{
  /* defaults */
  prop_0[0]=0;        /* by default beam propagates along z-axis */
  prop_0[1]=0;
  prop_0[2]=1;
  directory[0]=0;
  lambda=2*PI;
    /* initialize ref_index of scatterer */
  Nmat=1;
  ref_index[0][RE]=1.5;
  ref_index[0][IM]=0.0;
    /* initialize to null to determine further whether it is initialized */
  logfile=NULL;
  logname[0]=0;

  boxX=boxY=boxZ=UNDEF;
  sizeX=UNDEF;
  dpl=UNDEF;
  strcpy(run_name,"run");
  nTheta=UNDEF;
  eps=1.0e-5;
  shape=SH_SPHERE;
  strcpy(shapename,"sphere");
  store_int_field=FALSE;
  PolRelation=POL_LDR;
  ScatRelation=SQ_DRAINE;
  IntRelation=G_POINT_DIP;
  IterMethod=IT_CGNR;
  NoSymmetry=FALSE;
  symmetry_enforced=FALSE;
  prognose=FALSE;
  maxiter=UNDEF;
  jagged=1;
  beamtype=B_PLANE;
  strcpy(alldir_parms,"alldir_params.dat");
  strcpy(avg_parms,"avg_params.dat");
  strcpy(scat_grid_parms,"scat_params.dat");
  strcpy(chp_dir,"chpoint");
  chp_time=UNDEF;
  chp_type=CHP_NONE;
  orient_avg=FALSE;
  alph_deg=bet_deg=gam_deg=0.0;
  volcor=TRUE;
  reduced_FFT=TRUE;
  save_geom=FALSE;
  save_geom_fname[0]=0;
  yzplane=UNDEF;
  all_dir=FALSE;
  scat_grid=FALSE;
  phi_integr=FALSE;
  store_scat_grid=FALSE;
  calc_Cext=TRUE;
  calc_Cabs=TRUE;
  calc_Csca=FALSE;
  calc_vec=FALSE;
  calc_asym=FALSE;
  calc_mat_force=FALSE;
  store_force=FALSE;
  load_chpoint=FALSE;
  memory=0;
}

/*============================================================*/

void ParseParameters(int argc,char **argv)
  /* parses input parameters */
{
  int i,j,Narg;
  double doubTmp; /*buffer variable*/

  /* read command line */
  i=1;
  if (argc>1) do {
    /* get number of arguments
       parameter begins with '-' and then letter
       it enables use of negative numbers as subparameters */
    Narg=0;
    do {
      Narg++;
      if ((i+Narg)>=argc) break;
      if (argv[i+Narg][0]=='-') if ((argv[i+Narg][1]>='a' && argv[i+Narg][1]<='z') ||
         (argv[i+Narg][1]>='A' && argv[i+Narg][1]<='Z')) break;
    } while(TRUE);
    Narg--;

    if (strcmp(argv[i],"-alldir_inp")==0) {
      IllNarg("-alldir_inp",Narg,1);
      strcpy(alldir_parms,argv[++i]);
    }
    else if (strcmp(argv[i],"-asym")==0) {
      IllNarg("-asym",Narg,0);
      calc_asym = TRUE;
      calc_vec = TRUE;
      calc_Csca = TRUE;
    }
    else if (strcmp(argv[i],"-beam")==0) {
      if (Narg==0) {
        printz("No arguments are given to -beam option (at least 1 expected)\n");
        Stop(1);
      }
      i++;
      Narg--;
      if (strcmp(argv[i],"plane")==0) {
        beamtype=B_PLANE;
        IllNarg("'-beam plane'",Narg,0);
      }
      else {
	if (strcmp(argv[i],"buggy")==0) beamtype=B_BUGGY;
	else if (strcmp(argv[i],"barton1")==0) beamtype=B_BARTON1;
	else if (strcmp(argv[i],"barton3")==0) beamtype=B_BARTON3;
	else if (strcmp(argv[i],"barton5")==0) beamtype=B_BARTON5;
	else if (strcmp(argv[i],"davis1")==0) beamtype=B_DAVIS1;
	else if (strcmp(argv[i],"davis3")==0) beamtype=B_DAVIS3;
	else if (strcmp(argv[i],"lminus")==0) beamtype=B_LMINUS;
	else {
	  printz("Beam type '%s' is not supported\n",argv[i]);
          Stop(1);
        }
        if (Narg!=4) {
          printz("Illegal number of arguments (%d) to '-beam %s' option (4 expected)\n",
                 Narg,argv[i]);
          Stop(1);
        }
	sscanf(argv[++i],"%lf",&beam_w0);
	sscanf(argv[++i],"%lf",&beam_x0);
	sscanf(argv[++i],"%lf",&beam_y0);
	sscanf(argv[++i],"%lf",&beam_z0);
      }
    }
    else if (strcmp(argv[i],"-chp_dir")==0) {
      IllNarg("-chp_dir",Narg,1);
      strcpy(chp_dir,argv[++i]);
    }
    else if (strcmp(argv[i],"-chp_load")==0) {
      IllNarg("-chp_load",Narg,0);
      load_chpoint = TRUE;
    }
    else if (strcmp(argv[i],"-chp_type")==0) {
      IllNarg("-chp_type",Narg,1);
      i++;
      if (strcmp(argv[i],"normal")==0) chp_type=CHP_NORMAL;
      else if (strcmp(argv[i],"regular")==0) chp_type=CHP_REGULAR;
      else if (strcmp(argv[i],"always")==0) chp_type=CHP_ALWAYS;
      else {
        printz("Checkpoint type '%s' is not supported\n",argv[i]);
        Stop(1);
      }
    }
    else if (strcmp(argv[i],"-chpoint")==0) {
      IllNarg("-chpoint",Narg,1);
      chp_time=ScanTime(argv[++i]);
      if (chp_time<=0) {
        chp_time=UNDEF;
        if (chp_type==CHP_NONE) chp_type=CHP_ALWAYS;
      }
      else if (chp_type==CHP_NONE) chp_type=CHP_NORMAL;
    }
    else if (strcmp(argv[i],"-Cpr_mat")==0) {
      IllNarg("-Cpr_mat",Narg,0);
      calc_mat_force = TRUE;
    }
    else if (strcmp(argv[i],"-Csca")==0) {
      IllNarg("-Csca",Narg,0);
      calc_Csca = TRUE;
    }
    else if (strcmp(argv[i],"-dir")==0) {
      IllNarg("-dir",Narg,1);
      strcpy(directory,argv[++i]);
    }
    else if (strcmp(argv[i],"-dpl")==0) {
      IllNarg("-dpl",Narg,1);
      sscanf(argv[++i],"%lf",&dpl);
    }
    else if (strcmp(argv[i],"-eps")==0) {
      IllNarg("-eps",Narg,1);
      sscanf(argv[++i],"%lf",&doubTmp);
      eps=pow(10,-doubTmp);
    }
    else if (strcmp(argv[i],"-grid")==0) {
      if (Narg!=1 && Narg!=3) {
	printz("Illegal number of arguments (%d) to -grid option (1 or 3 expected)\n",Narg);
	Stop(1);
      }
      sscanf(argv[++i],"%i",&boxX);  /* boxes are further multiplied by jagged if needed */
      if (Narg==3) {
        sscanf(argv[++i],"%i",&boxY);
        sscanf(argv[++i],"%i",&boxZ);
      }
    }
    else if (strcmp(argv[i],"-int")==0) {
      IllNarg("-int",Narg,1);
      i++;
      if (strcmp(argv[i],"poi")==0) IntRelation=G_POINT_DIP;
      else if (strcmp(argv[i],"so")==0) IntRelation=G_SO;
      else {
        printz("Interaction term prescription '%s' is not supported\n",argv[i]);
        Stop(1);
      }
    }
    else if (strcmp(argv[i],"-iter")==0) {
      IllNarg("-iter",Narg,1);
      i++;
      if (strcmp(argv[i],"cgnr")==0) IterMethod=IT_CGNR;
      else if (strcmp(argv[i],"bicgstab")==0) IterMethod=IT_BICGSTAB;
      else if (strcmp(argv[i],"bicg")==0) IterMethod=IT_BICG_CS;
      else if (strcmp(argv[i],"qmr")==0) IterMethod=IT_QMR_CS;
      else {
        printz("Iterative method '%s' is not supported\n",argv[i]);
        Stop(1);
      }
    }
    else if (strcmp(argv[i],"-jagged")==0) {
      IllNarg("-jagged",Narg,1);
      sscanf(argv[++i],"%d",&jagged);
    }
    else if (strcmp(argv[i],"-lambda")==0) {
      IllNarg("-lambda",Narg,1);
      sscanf(argv[++i],"%lf",&lambda);
    }
    else if (strcmp(argv[i],"-m")==0) {
      if (Narg%2!=0 || Narg==0) {
	printz("Illegal number of arguments (%d) to -m option (even expected)\n",Narg);
	Stop(1);
      }
      Nmat=Narg/2;
      if (Nmat>MAX_NMAT) {
        printz("Too many materials (%d), maximum %d are supported.\n"\
               "Increase parameter MAX_NMAT and recompile\n",Nmat,MAX_NMAT);
        Stop(1);
      }
      for (j=0;j<Nmat;j++) {
        sscanf(argv[++i],"%lf",&ref_index[j][RE]);
        sscanf(argv[++i],"%lf",&ref_index[j][IM]);
      }
    }
    else if (strcmp(argv[i],"-maxiter")==0) {
      IllNarg("-maxiter",Narg,1);
      sscanf(argv[++i],"%i",&maxiter);
    }
    else if (strcmp(argv[i],"-no_reduced_fft")==0) {
      IllNarg("-no_reduced_fft",Narg,0);
      reduced_FFT=FALSE;
    }
    else if (strcmp(argv[i],"-no_vol_cor")==0) {
      IllNarg("-no_vol_cor",Narg,0);
      volcor=FALSE;
    }
    else if (strcmp(argv[i],"-nosym")==0) {
      IllNarg("-nosym",Narg,0);
      NoSymmetry=TRUE;
    }
    else if (strcmp(argv[i],"-ntheta")==0) {
      IllNarg("-ntheta",Narg,1);
      sscanf(argv[++i],"%i",&nTheta);
      nTheta++;
    }
    else if (strcmp(argv[i],"-orient")==0) {
      if (Narg==0) {
        printz("No arguments are given to -orient option (at least 1 expected)\n");
        Stop(1);
      }
      i++;
      if (strcmp(argv[i],"avg")==0) {
        if (Narg>2) {
          printz("Illegal number of arguments (%d) to '-orient avg' option (0 or 1 expected)\n",
                 Narg-1);
          Stop(1);
        }
        orient_avg=TRUE;
        if (Narg==2) strcpy(avg_parms,argv[++i]);
      }
      else {
        if (Narg!=3) {
          printz("Illegal number of numerical arguments (%d) to -orient option\n"\
                 "  or first argument (%s) is incorrect (can be 'avg')\n",Narg,argv[i]);
          Stop(1);
        }
        sscanf(argv[i],"%lf",&alph_deg);
	sscanf(argv[++i],"%lf",&bet_deg);
	sscanf(argv[++i],"%lf",&gam_deg);
      }
    }
    else if (strcmp(argv[i],"-phi_integr")==0) {
      IllNarg("-phi_integr",Narg,1);
      phi_integr = TRUE;
      sscanf(argv[++i],"%d",&phi_int_type);
    }
    else if (strcmp(argv[i],"-pol")==0) {
      if (Narg!=1 && Narg!=2) {
        printz("Illegal number of arguments (%d) to -pol option (1 or 2 expected)\n",Narg);
        Stop(1);
      }
      i++;
      if (strcmp(argv[i],"cm")==0) PolRelation=POL_CM;
      else if (strcmp(argv[i],"rrc")==0) PolRelation=POL_RR;
      else if (strcmp(argv[i],"ldr")==0) PolRelation=POL_LDR;
      else if (strcmp(argv[i],"cldr")==0) PolRelation=POL_CLDR;
      else if (strcmp(argv[i],"so")==0) PolRelation=POL_SO;
      else {
        printz("Polarization Relation '%s' is not supported\n",argv[i]);
        Stop(1);
      }
      if (Narg==2) {
        i++;
        if (strcmp(argv[i],"avgpol")==0) avg_inc_pol=TRUE;
        else {
          printz("Unknown argument '%s' to '-pol %s' option\n",argv[i],argv[i-1]);
          Stop(1);
        }
      }
    }
    else if (strcmp(argv[i],"-prognose")==0) {
      IllNarg("-prognose",Narg,0);
      prognose=TRUE;
      strcpy(run_name,"test");
    }
    else if (strcmp(argv[i],"-prop")==0) {
      IllNarg("-prop",Narg,3);
      sscanf(argv[++i],"%lf",prop_0);
      sscanf(argv[++i],"%lf",&prop_0[1]);
      sscanf(argv[++i],"%lf",&prop_0[2]);
      doubTmp=DotProd(prop_0,prop_0);
      if (doubTmp==0) {
        printz("Given propagation vector is null\n");
        Stop(1);
      }
      doubTmp=1/sqrt(doubTmp);
      prop_0[0]*=doubTmp;
      prop_0[1]*=doubTmp;
      prop_0[2]*=doubTmp;
    }
    else if (strcmp(argv[i],"-save_geom")==0) {
      if (Narg>1) {
        printz("Illegal number of arguments (%d) to -save_geom option (0 or 1 expected)\n",Narg);
        Stop(1);
      }
      save_geom=TRUE;
      if (Narg==1) sscanf(argv[++i],"%s",save_geom_fname);
    }
    else if (strcmp(argv[i],"-scat")==0) {
      IllNarg("-scat",Narg,1);
      i++;
      if (strcmp(argv[i],"dr")==0) ScatRelation=SQ_DRAINE;
      else if (strcmp(argv[i],"so")==0) ScatRelation=SQ_SO;
      else {
        printz("Scattering Quantities Relation '%s' is not supported\n",argv[i]);
        Stop(1);
      }
    }
    else if (strcmp(argv[i],"-scat_grid_inp")==0) {
      IllNarg("-scat_grid_inp",Narg,1);
      strcpy(scat_grid_parms,argv[++i]);
    }
    else if (strcmp(argv[i],"-shape")==0) {
      i++;
      Narg--;
      strcpy(shapename,argv[i]);
      if (strcmp(argv[i],"read")==0) {
        IllNarg("'-shape read'",Narg,1);
        shape=SH_READ;
        strcpy(aggregate_file,argv[++i]);
      }
      else if (strcmp(argv[i],"coated")==0) {
        if (Narg!=1 && Narg!=4) {
          printz("Illegal number of arguments (%d) to '-shape coated' option (1 or 4 expected)\n",
                 Narg);
          Stop(1);
        }
        shape=SH_COATED;
        sh_Npars=Narg;
        for (j=0;j<Narg;j++) sscanf(argv[++i],"%lf",sh_pars+j);
      }
      else {
/*    else if (strcmp(argv[i],"-sdisk_rot")==0) {
      if (Narg!=3) {
        printz("Illegal number of arguments (%d) to -sdisk_rot option\n",Narg);
        Stop(1);
      }
      shape=SH_SDISK_ROT;
      sscanf(argv[++i],"%lf",&aspect_r);
      sscanf(argv[++i],"%lf",&betaY);
      sscanf(argv[++i],"%lf",&betaZ);
    }  */
    /*    else if (strcmp(argv[i],"-prisma")==0) {
      shape=SH_PRISMA;
    }  */
        if (strcmp(argv[i],"box")==0) {
          sh_Npars=0;
          shape=SH_BOX;
        }
        else if (strcmp(argv[i],"cylinder")==0) {
          sh_Npars=1;
          shape=SH_CYLINDER;
        }
        else if (strcmp(argv[i],"ellipsoid")==0) {
          sh_Npars=2;
          shape=SH_ELLIPSOID;
        }
        else if (strcmp(argv[i],"line")==0) {
          sh_Npars=0;
          shape=SH_LINE;
        }
        else if (strcmp(argv[i],"rbc")==0) {
          sh_Npars=3;
          shape=SH_RBC;
        }
        else if (strcmp(argv[i],"sphere")==0) {
          sh_Npars=0;
          shape=SH_SPHERE;
        }
        else if (strcmp(argv[i],"spherebox")==0) {
          sh_Npars=1;
          shape=SH_SPHEREBOX;
        }
        else {
	  printz("Shape '%s' is not supported\n",argv[i]);
          Stop(1);
        }
        if (Narg!=sh_Npars) {
          printz("Illegal number of arguments (%d) to '-shape %s' option (%d expected)\n",
                 Narg,argv[i],sh_Npars);
          Stop(1);
        }
        for (j=0;j<Narg;j++) sscanf(argv[++i],"%lf",sh_pars+j);
      }
    }
    else if (strcmp(argv[i],"-size")==0) {
      IllNarg("-size",Narg,1);
      sscanf(argv[++i],"%lf",&sizeX);
    }
    else if (strcmp(argv[i],"-store_force")==0) {
      IllNarg("-store_force",Narg,0);
      store_force = TRUE;
    }
    else if (strcmp(argv[i],"-store_int_field")==0) {
      IllNarg("-store_int_field",Narg,0);
      store_int_field=TRUE;
    }
    else if (strcmp(argv[i],"-store_scat_grid")==0) {
      IllNarg("-store_scat_grid",Narg,0);
      store_scat_grid = TRUE;
    }
    else if (strcmp(argv[i],"-sym_enf")==0) {
      IllNarg("-sym_enf",Narg,0);
      symmetry_enforced=TRUE;
    }
    else if (strcmp(argv[i],"-test")==0) {
      IllNarg("-test",Narg,0);
      strcpy(run_name,"test");
    }
    else if (strcmp(argv[i],"-vec")==0) {
      IllNarg("-vec",Narg,0);
      calc_vec = TRUE;
    }
    else if (strcmp(argv[i],"-yz")==0) {
      IllNarg("-yz",Narg,0);
      yzplane = TRUE;
    }
    else {
      printz("adda: illegal option \"%s\"\n"\
             "Usage: adda [-lambda %%f] [-grid %%i [%%i %%i]] [-m %%f %%f [...]]\n"\
             "    [-dpl %%f] [-size %%f] [-test] [-ntheta %%i] [-maxiter %%i] [-dir %%s]\n"\
             "    [-yz] [-alldir_inp %%s] [-store_scat_grid] [-phi_integr %%d]\n"\
             "    [-scat_grid_inp %%s] [-Csca] [-vec] [-asym] [-Cpr_mat] [-store_force]\n"\
             "    [-jagged %%d] [-store_int_field] [-prognose] [-prop %%f %%f %%f]\n"\
             "    [-nosym] [-sym_enf] [-eps %%f]\n"\
             "    [-shape box|coated|cylinder|ellipsoid|line|rbc|sphere|spherebox\n"\
             "     [...]|read %%s] \n"\
             "    [-beam plane|lminus|buggy|davis1|davis3|barton1|barton3|barton5 [...]]\n"\
             "    [-pol cm|rrc|ldr|cldr|so] [-scat dr|so] [-int poi|so]\n"\
             "    [-orient avg[ %%s]|%%f %%f %%f] [-iter cgnr|bicgstab|bicg|qmr]\n"\
             "    [-no_vol_cor] [-no_reduced_fft] [-save_geom [%%s]]\n"\
             "    [-chpoint %%d[h,m,s]] [-chp_type normal|regular|always] [chp_load]\n"\
             "    [chp_dir %%s]\n"\
             "    \n",argv[i]);
      Stop(1);
    }
    i++;
  } while(i<argc);  /* end of reading command line arguments */
}

/*============================================================*/

void VariablesInterconnect(void)
  /* finish parameters initialization based on their interconnections */
{
  double temp;
  /* parameter interconnections */
  if (prop_0[2]!=1 && orient_avg) {
    printz("-prop and '-orient avg' can not be used together\n");
    Stop(1);
  }
  if (chp_time==UNDEF && chp_type!=CHP_NONE && chp_type!=CHP_ALWAYS) {
    printz("You must specify time for this checkpoint type\n");
    Stop(1);
  }
  if (IntRelation==G_SO) reduced_FFT=FALSE;
  /* scale boxes by jagged */
  if (jagged!=1) {
    if (boxX!=UNDEF) boxX*=jagged;
    if (boxY!=UNDEF) boxY*=jagged;
    if (boxZ!=UNDEF) boxZ*=jagged;
  }
  if (calc_Csca || calc_vec) all_dir = TRUE;
  if (store_scat_grid || phi_integr) {
    scat_grid = TRUE;
    if (yzplane==UNDEF) yzplane = FALSE;
  }
  else if (yzplane==UNDEF) yzplane = TRUE;

  /*determine two incident polarizations. Equivalent to rotation of X,Y,Z basis by
  angles Theta and Phi from (0,0,1) to given propagation vector */
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
  /* initialize averaging over orientation */
  if (orient_avg) {
    ReadAvgParms(avg_parms);
    NoSymmetry=TRUE;
    avg_inc_pol=TRUE;
  }
  else {
    /* else - initialize rotation stuff */
    InitRotation();
    if (prop[2]!=1) NoSymmetry=TRUE;
  }
}

/*============================================================*/

void DirectoryLog(int argc,char **argv)
   /* create input directory and start logfile */
{
  int  i,Nexp;
  FILE *Nexpfile;
  char *ptmp,*ptmp2;
  FILEHANDLE lockid;
  char sbuffer[500];
#ifdef _WIN32  /* Windows; for obtaining computer name */
  TCHAR cname[MAX_COMPUTERNAME_LENGTH+1];
  DWORD cname_size=MAX_COMPUTERNAME_LENGTH+1;
#endif

  /* devise directory name (for output files) */
  if (directory[0]==0) {
    /* create lock file */
    if (ringid==ROOT) lockid=CreateLockFile("ExpCount.lck");
    Synchronize();
    /* read ExpCount */
    if ((Nexpfile=fopen("ExpCount","r"))!=NULL) {
      fscanf(Nexpfile,"%i",&Nexp);
      fclose(Nexpfile);
    }
    else Nexp=0;
    /* wait for all processors to finish reading Nexpfile */
    Synchronize();
    /* put new number in Nexpfile */
    if (ringid==ROOT) {
      if ((Nexpfile=fopen("ExpCount","w"))==NULL)
        LogError(EC_ERROR,ONE,POSIT,"Failed to open file 'ExpCount'");
      fprintz(Nexpfile,"%i",Nexp+1);
      fclose(Nexpfile);
      /* unlock */
      RemoveLockFile(lockid,"ExpCount.lck");
    }
    /* create directory name */
    sprintf(sbuffer,"m%.4g",ref_index[0][RE]);
    ptmp=strchr(sbuffer,'.');
    if (ptmp!=NULL) *ptmp='_';
    sprintf(directory,"%s%03i_%s_g%i%s",run_name,Nexp,shapename,boxX,sbuffer);
#ifdef PARALLEL
    /* add PBS job id to the directory name if available */
    ptmp=getenv("PBS_JOBID");
    if (ptmp!=NULL && *ptmp!=0) {
      /* jobid is truncated at first "." */
      if ((ptmp2=strchr(ptmp,'.'))!=NULL) *ptmp2=0;
      sprintf(directory+strlen(directory),"id%s",ptmp);
    }
#endif
  }
  /* make new directory and print info */
  if (ringid==ROOT) {
    strcpy(sbuffer,"mkdir ");
    strcat(sbuffer,directory);
    system(sbuffer);
    printf("all data is saved in '%s'\n",directory);
  }
  /* make logname; do it for all processors to enable additional logging in LogError */
  strcpy(logname,directory);
  strcat(logname,"/log");
  /* start logfile */
  if (ringid==ROOT) {
    /* open logfille */
    if ((logfile=fopen(logname,"w"))==NULL)
      LogError(EC_ERROR,ONE,POSIT,"Failed to open file 'logfile'");
    /* log version number */
    fprintf(logfile,"Generated by ADDA v.%s\n",version);
#ifdef PARALLEL
    /* write number of processors */
    fprintf(logfile,"The program was run on: %d processors",nprocs);
    /* add PBS host name if present */
    if ((ptmp=getenv("PBS_O_HOST"))!=NULL) fprintf(logfile," from %s\n",ptmp);
    else fprintf(logfile,"\n");
#else /* sequential */
    /* write computer name */
# ifdef _WIN32 /* Windows */
    GetComputerName(cname,&cname_size);
    fprintf(logfile,"The program was run on: %s\n",cname);
# else /* UNIX */
    if ((ptmp=getenv("HOST"))!=NULL)
      fprintf(logfile,"The program was run on: %s\n",ptmp);
# endif
#endif
    /* log command line */
    fprintf(logfile,"command: '");
    for(i=0;i<argc;i++) fprintf(logfile,"%s ",argv[i]);
    fprintf(logfile,"'\n");
  }
}

/*============================================================*/

void PrintInfo(void)
  /* print info to stdout and logfile */
{
   int  i;
   char sbuffer[200];

   if (ringid==ROOT) {
    /* print basic parameters */
    printf("lambda: %g   m0: %g%+gi   Dipoles/lambda: %g\n",
           lambda,ref_index[0][RE],ref_index[0][IM],dpl);
    printf("Required relative error: %g\n",eps);
    printf("Total number of occupied dipoles: %d\n",nvoid_Ndip);
    /* log basic parameters */
    fprintf(logfile,"lambda: %g\n",lambda);
    fprintf(logfile,"shape: ");
    fprintf(logfile,sh_form_str,sizeX);
    fprintf(logfile,"\nbox dimensions: %ix%ix%i\n",boxX,boxY,boxZ);
    fprintf(logfile,"refractive index: ");
    if (Nmat==1) fprintf(logfile,"%g%+gi\n",ref_index[0][RE],ref_index[0][IM]);
    else {
      fprintf(logfile,"1. %g%+gi\n",ref_index[0][RE],ref_index[0][IM]);
      for (i=1;i<Nmat;i++) fprintf(logfile,
         "                  %d. %g%+gi\n",i+1,ref_index[i][RE],ref_index[i][IM]);
    }
    fprintf(logfile,"Dipoles/lambda: %g\n",dpl);
    if (volcor_used) fprintf(logfile,"      (Volume correction used)\n");
    fprintf(logfile,"Required relative error: %g\n",eps);
    fprintf(logfile,"Total number of occupied dipoles: %d\n",nvoid_Ndip);
    fprintf(logfile,"Volume-equivalent size parameter: %g\n",ka_eq);
    /* log incident polarization */
    fprintf(logfile,"\nIncident propagation vector: (%g,%g,%g)\n",
            prop_0[0],prop_0[1],prop_0[2]);
    fprintf(logfile,"Incident polarization Y(par): (%g,%g,%g)\n",
            incPolY_0[0],incPolY_0[1],incPolY_0[2]);
    fprintf(logfile,"Incident polarization X(per): (%g,%g,%g)\n\n",
            incPolX_0[0],incPolX_0[1],incPolX_0[2]);
    /* log particle orientation */
    if (orient_avg) fprintf(logfile,"Particle orientation - averaged\n%s\n",avg_string);
    else fprintf(logfile,"Particle orientation (deg): alpha=%g, beta=%g, gamma=%g\n\n",
                 alph_deg,bet_deg,gam_deg);
    /* log incident polarization after transformation */
    if (alph_deg!=0 || bet_deg!=0 || gam_deg!=0) {
      fprintf(logfile,"After transformation to particle reference frame:\n");
      fprintf(logfile,"New incident propagation vector: (%g,%g,%g)\n",
              prop[0],prop[1],prop[2]);
      fprintf(logfile,"New incident polarization Y(par): (%g,%g,%g)\n",
              incPolY[0],incPolY[1],incPolY[2]);
      fprintf(logfile,"New incident polarization X(per): (%g,%g,%g)\n\n",
              incPolX[0],incPolX[1],incPolX[2]);
    }
    /* log Polarization relation */
    if (PolRelation==POL_CM)
      fprintf(logfile,"Polarization relation: 'Clausius-Mossotti'\n");
    else if (PolRelation==POL_RR)
      fprintf(logfile,"Polarization relation: 'Radiative Reaction Correction'\n");
    else if (PolRelation==POL_LDR) {
      fprintf(logfile,"Polarization relation: 'Lattice Dispersion Relation'");
      if (avg_inc_pol) fprintf(logfile," (averaged over incident polarization)");
      fprintf(logfile,"\n");
    }
    else if (PolRelation==POL_CLDR)
      fprintf(logfile,"Polarization relation: 'Corrected Lattice Dispersion Relation'\n");
    else if (PolRelation==POL_SO)
      fprintf(logfile,"Polarization relation: 'Second Order'\n");
    /* log Scattering Quantities formulae */
    if (ScatRelation==SQ_DRAINE)
      fprintf(logfile,"Scattering quantities formulae: 'by Draine'\n");
    else if (ScatRelation==SQ_SO)
      fprintf(logfile,"Scattering quantities formulae: 'Second Order'\n");
    /* log Interaction term prescription */
    if (IntRelation==G_POINT_DIP)
      fprintf(logfile,"Interaction term prescription: 'as Point dipoles'\n");
    else if (IntRelation==G_SO)
      fprintf(logfile,"Interaction term prescription: 'Second Order'\n");
    /* log FFT method */
#ifdef FFTW3
    fprintf(logfile,"FFT algorithm: FFTW3\n");
#else
# ifdef FFT_TEMPERTON
    fprintf(logfile,"FFT algorithm: by C.Temperton\n");
# endif
#endif
    /* log Iterative Method */
    if (IterMethod==IT_CGNR)
      fprintf(logfile,"Iterative Method: CGNR\n");
    else if (IterMethod==IT_BICGSTAB)
      fprintf(logfile,"Iterative Method: Bi-CG Stabilized\n");
    else if (IterMethod==IT_BICG_CS)
      fprintf(logfile,"Iterative Method: Bi-CG (complex symmetric)\n");
    else if (IterMethod==IT_QMR_CS)
      fprintf(logfile,"Iterative Method: QMR (complex symmetric)\n");
    /* log Symmetry options */
    if (symmetry_enforced) fprintf(logfile,"Symmetry is enforced by user (warning!)\n");
    else if (NoSymmetry) fprintf(logfile,"No symmetries are used\n");
    /* log Checkpoint options */
    if (load_chpoint) fprintf(logfile,"Simulation is continued from a checkpoint\n");
    if (chp_type!=CHP_NONE) {
      fprintf(logfile,"Checkpoint is turned on:\n");
      if (chp_type==CHP_NORMAL) fprintf(logfile,"    type = normal\n");
      else if (chp_type==CHP_REGULAR) fprintf(logfile,"    type = regular\n");
      else if (chp_type==CHP_ALWAYS) fprintf(logfile,"    type = always\n");
      if (chp_time==UNDEF) fprintf(logfile,"    time = no limit\n");
      else {
        PrintTime(sbuffer,&chp_time);
        fprintf(logfile,"    time = %s(%d sec)\n",sbuffer,chp_time);
      }
    }
    if (load_chpoint || chp_type!=CHP_NONE)
      fprintf(logfile,"    directory = '%s'\n",chp_dir);
  }
}

