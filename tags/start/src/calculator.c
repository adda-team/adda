/* FILE: calculator.c
 * AUTH: Alfons Hoekstra (alfons@fwi.uva.nl)
 * HIST:
 * DATE:
 */
#include <stdio.h>
#include "cmplx.h"
#include "types.h"
#include "Romberg.h"
#include "const.h"
#include <math.h>
#include "comm.h"
#include "debug.h"

/* LEVELS, LAMBDA, DIPS_PER_LAMBDA are read from the
 * command-line by main
 */
double DIPS_PER_LAMBDA;
double LAMBDA;

extern dcomplex *dCvector(int nl,int nh);
extern REAL **dmatrix(int nrl,int nrh,int ncl,int nch);
extern void free_dmatrix(REAL **m, int nrl, int nrh, int ncl);

extern int MyProcId;
extern dcomplex *x, *p, *r;
extern dcomplex *buffer;	/* used to gather vectors in */
extern dcomplex *Eper, *Epar, *Egrid, *Einc;
extern FILE *logfile;
extern int prognose;
extern doublecomplex ref_index[10];
extern doublecomplex cc[10];
extern Parms_1D parms[];
extern int yzplane;
extern int xzplane;
extern int all_dir;
extern char directory[200];

dcomplex *eper_buffer;
int nDip;             /* number of dipoles */
int nTheta;	      /* number of angles in scattering profile */
int nRing;            /* number of processors in Ring */
int RingId;           /* ID of this processor in the Ring */
int nlocalDip;        /* number of local dipoles */
int nlocalRows;       /* number of local rows of decomposition */
dcomplex *Avecbuffer; /* buffer to hold result of local Matvec */
REAL **DipoleCoord;   /* matrix to hold the coordinates of the dipoles */
short int *position;  /* position of the dipoles */
int *material;        /* material: index for cc */
double WaveNum;       /* wavenumber of incident light */
double eps;	      /* accuracy epsilon */
doublecomplex CoupleConstant; /* coupling between Polarization and E-field */
double lambda;                /* wavelenght of incident light */
int memory=0;         /* total memory usage in bytes */

doublecomplex coupleconstant(doublecomplex mrel)
{
  doublecomplex CoupleConstant,tempa,tempb;
  double temp;
  extern double dplX,dplY,dplZ;
  extern int PolRelation;
  extern double gridspaceX,gridspaceY,gridspaceZ;
  double dips_lambda = DIPS_PER_LAMBDA;

  /* calculate the couple constant */
  temp = (3 * lambda*lambda*lambda/dplX/dplY/dplZ) / (4 * PI);

  CoupleConstant.r = temp;
  CoupleConstant.i = temp;
  tempa.r = mrel.r * mrel.r - mrel.i * mrel.i -1.0;
  tempa.i = 2.0 * mrel.r * mrel.i;
  tempb.r = mrel.r * mrel.r - mrel.i * mrel.i + 2.0;
  tempb.i = 2.0 * mrel.r * mrel.i;
  temp = tempb.r * tempb.r + tempb.i * tempb.i;
  CoupleConstant.r *= (tempa.r * tempb.r + tempa.i * tempb.i) / temp;
  CoupleConstant.i *= (-tempa.r * tempb.i + tempa.i * tempb.r) / temp;

  fprintz(logfile,"CoupleConstant:%g+%gi\n",CoupleConstant.r,CoupleConstant.i);
  if (PolRelation==LDR || PolRelation==RADCOR) {
    double gridspace;
    dcomplex c;                  /* correction to CM */
    double S,kd,norm;
    dcomplex t1,t2;
    dcomplex m2,cc;
    
    if (gridspaceX!=gridspaceY || gridspaceX!=gridspaceZ) {
       printz("unable to use ldr with unequal gridspacing\n");
       stop(1);
     }
    cc.r=CoupleConstant.r; cc.i=CoupleConstant.i;

    gridspace=lambda / dips_lambda;
    kd=2*PI/dips_lambda;
    S=0;
    if (PolRelation==RADCOR) S=0;

    m2.r=mrel.r*mrel.r-mrel.i*mrel.i;      m2.i=2*mrel.r*mrel.i;
    t1.r=LDR_B1;                   t1.i=0.0;                   /* t1=b1                      */
    t1.r+=LDR_B2*m2.r;             t1.i+=LDR_B2*m2.i;          /* t1=b1+b2*m^2               */
    t1.r+=LDR_B3*m2.r*S;           t1.i+=LDR_B3*m2.i*S;        /* t1=b1+b2*m^2+b3*m^2*S      */
    t1.r*=kd*kd;                   t1.i*=kd*kd;                /* t1*=kd^2                   */

    if (PolRelation==RADCOR) t1.r = t1.i = 0.0; /* ignore LDR corrections */

    t1.r-=0.0;                     t1.i-=2*kd*kd*kd/3;         /* t1=t1-2/3*i*kd^3           */
    t1.r/=gridspace*gridspace*gridspace;
    t1.i/=gridspace*gridspace*gridspace;                       /* t1/=d*d*d                  */
    t2.r=t1.r*cc.r - t1.i*cc.i +1; t2.i=t1.i*cc.r + t1.r*cc.i; /* t2=cc*t1+1                 */
    printz("t2:%g %g\n",t2.r,t2.i);

    t1.r=cc.r*t2.r + cc.i* t2.i;   t1.i=cc.i*t2.r - cc.r*t2.i; /* t1=cc*t2*                  */
    norm=t2.r*t2.r - t2.i*t2.i;
    CoupleConstant.r=t1.r/norm;    CoupleConstant.i=t1.i/norm; /* cc=cc/t2                   */
    
  }
  printz("CoupleConstant:%g+%gi\n",CoupleConstant.r,CoupleConstant.i);
  fprintz(logfile,"CoupleConstant:%g+%gi\n",CoupleConstant.r,CoupleConstant.i);
  return(CoupleConstant);
}

void
Calculator (void)
{
  int tempdip;
  double dips_lambda;		/* number of dipoles per wavelenght */
  dcomplex mrel;		/* relative refractive index of particle */
  char *Fname = __FILE__;
  double temp;
  dcomplex tempa, tempb;
  int i;

  extern int shape;             /* shape of the particle */
  extern int boxX,boxY,boxZ;
  extern double eps0;
  extern int PolRelation;
  extern int symR,jagged;
  extern double dplX,dplY,dplZ,sizeX,sizeY,sizeZ;
  extern double gridspaceX,gridspaceY,gridspaceZ;
  extern int Nmat,mat_count[10];

  /* initialize some variables */
  lambda      = LAMBDA;
  dips_lambda = DIPS_PER_LAMBDA;
  WaveNum     = 2.0 * PI / lambda;
  eps=eps0;
  for(i=0;i<Nmat;i++) cc[i]=coupleconstant(ref_index[i]);
  /* find position in the ring */
  nRing  = nprocs;
  RingId = ringid;

  /* determine the number dipoles in this processor */
  if (nDip%nRing == 0)
    nlocalDip = nDip/nRing;
  else if (RingId < nDip%nRing)
    nlocalDip = nDip/nRing + 1;
  else
    nlocalDip = nDip/nRing;
  init_Dmatrix(boxX,boxY,boxZ,lambda,dips_lambda,WaveNum);
  
  nlocalRows = 3 * local_Ndip;
  
  if (prognose==false) if ((x = dCvector(0, nlocalRows -1)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc x",
              MyProcId, RingId);
    AbortServer (1);
  }
  memory+=sizeof(dcomplex)*nlocalRows;
  
  if (prognose==false) if ((r = dCvector(0, nlocalRows -1))== NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc r",
              MyProcId, RingId);
    AbortServer (1);
  }
  memory+=sizeof(dcomplex)*nlocalRows;

  if (prognose==false) if ((p = dCvector(0, nlocalRows -1)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc p",
              MyProcId, RingId);
    AbortServer (1);
  }
  memory+=sizeof(dcomplex)*nlocalRows;

  if (prognose==false) if ((Einc = dCvector(0, nlocalRows -1)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc Einc",
              MyProcId, RingId);
    AbortServer (1);
  }
  memory+=sizeof(dcomplex)*nlocalRows;

  /* since we might need to hold nRing*(nDip/nRing+1) elements in */
  /* the buffer, we increase its size (was 3*nDip) */
  if (prognose==false) if ((buffer = dCvector(0, nlocalRows -1)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc buffer",
              MyProcId, RingId);
    AbortServer (1);
  }
  memory+=nlocalRows*sizeof(dcomplex);
  Avecbuffer=buffer;
  if ((eper_buffer = dCvector(0, nTheta*2)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc Eperb",
              MyProcId, RingId);
    AbortServer (1);
  }
  
  if (yzplane || xzplane)
    {
      if (prognose == false)
	{
	  if ((Eper = dCvector(0, nTheta*2)) == NULL) {
	    LogError (EC_ERROR, Fname,
		      "processor %d, ringID %d, could not malloc Eper",
		      MyProcId, RingId);
	    AbortServer (1);
	  }
	  
	  if ((Epar = dCvector(0, nTheta*2)) == NULL) {
	    LogError (EC_ERROR, Fname,
		      "processor %d, ringID %d, could not malloc Epar",
		      MyProcId, RingId);
	    AbortServer (1);
	  }
	}
      memory+=2*nTheta*sizeof(dcomplex);
    }
  
  if (all_dir)
    {
      set_Parms();

      if (prognose == false) 
	{
	  if ((Egrid = dCvector(0,3*parms[THETA].Grid_size*parms[PHI].Grid_size*sizeof(dcomplex))) == NULL) {
	    LogError (EC_ERROR, Fname,
		      "processor %d, ringID %d, could not malloc Egrid",
		  MyProcId, RingId);
	    AbortServer (1);
	  }
	}
      memory+=3*parms[THETA].Grid_size*parms[PHI].Grid_size;
    }

  if (prognose==false) if ((DipoleCoord = dmatrix(0, local_Ndip, 0, 2)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc DipoleCoord",
              MyProcId, RingId);
    AbortServer (1);
  }
  memory+=(sizeof(REAL)*3+sizeof(int))*local_Ndip;
  
  if (prognose==false) if ((position = (short int *) malloc(3*sizeof(short int)*local_Ndip)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc position",
              MyProcId, RingId);
    AbortServer (1);
  }
  memory+=3*sizeof(short int)*local_Ndip;
  
  if (prognose==false) if ((material = (int *) malloc(sizeof(int)*local_Ndip)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc position",
              MyProcId, RingId);
    AbortServer (1);
  }
  memory+=sizeof(char)*local_Ndip;

  printz("Total memory usage for all matrices:%.1f Mb\n",memory/1048576.0);
  if (prognose==true) {
    double flop,m;
    extern int gridX,gridY,gridZ,Ngrid;
    m=sqrt(ref_index[0].r*ref_index[0].r+ref_index[0].i*ref_index[0].i);
    printz("Prognose: total memory: %.1f Mb\n",memory/1048576.0);
    printz("Prognose: required ram: %.1f Mb\n",sizeof(dcomplex)*Ngrid/1048576.0);
    flop=Ngrid*log((double) Ngrid)*(m-1)/dplZ;
    printz("Prognose: required time: %.2g flop\n",flop);
    stop(1);
  }

  /* generate the matrix with dipole coordinates, do this in all processors */
  tempdip = make_particle (DipoleCoord,true,shape,jagged);
  for(i=0;i<Nmat;i++) mat_count[i]=0;
  for(i=0;i<local_Ndip;i++) mat_count[material[i]]++;
  for(i=0;i<Nmat;i++) printz("material:%i has locally %i dipoles\n",i,mat_count[i]);
  
  
  if (nDip != tempdip) { 
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, wrong number of dipoles: %d",
              MyProcId, RingId,tempdip);
    AbortServer (1);
  }
  
  /* calculate scattered field for y - polarized incident light */
  if(RingId == 0 ) fprintz(logfile,"here we go, calc Y\n\n");
  
  if (symR==true&&all_dir==false) {
    CalculateE('Y',PAR_AND_PER);
  }
  else { /* no rotational symmetry */
  /* in case of all_dir we run twice to get the full electric field */
  /* with incoming light in X and Y direction. In case of rotational */
  /* symmetry this is not needed but requires lots more programming */
  /* so we leave this optimization to a later time. */

    CalculateE('Y',NORMAL);
    if(RingId == 0 ) fprintz(logfile,"\nhere we go, calc X\n\n");
    CalculateE('X',NORMAL);
  }
  mueller_matrix();

  /* tidy up */
  /* lets not do this for a while....
   * runing MPI with this sometimes gives strange results, 
   * no time to find out what goes wrong....
  free_FFT_Dmat();
  free(p);
  free(r);
  free(eper_buffer);
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



