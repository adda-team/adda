/* FILE: fft.c
 * AUTH: Michel Grimminck (using old code by Alfons Hoekstra)
 * DESCR: Initialization of all FFT for matrix-vector products
 *        and FFT itself     
 *        A lot of indirect indexing used - way to optimize
 *
 *        Currently is developed by Maxim Yurkin
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "cmplx.h"
#include "const.h"
#include "comm.h"
#include "memory.h"
#include "debug.h"
#include "fft.h"
#include "timing.h"

extern clock_t Timing_FFT_Init,Timing_Dm_Init;

extern double gridspace;
extern int prognose;
extern int boxX,boxY,boxZ;
extern double WaveNum;
extern FILE *logfile;
extern int memory;
extern double dpl;
extern int reduced_FFT;

doublecomplex *Dmatrix;      /* holds FFT of the interaction matrix */
/* used in matvec */
doublecomplex *Xmatrix;      /* holds input vector (on expanded grid) to matvec */
doublecomplex *slices;       /* used in inner cycle of matvec - holds 3 components (for fixed x) */
doublecomplex *slices_tr;    /* additional storage space for slices to accelerate transpose */

/* D2 matrix and its two slices; used only temporary for init_Dmatrix */
doublecomplex *slice,*slice_tr,*D2matrix;

double *BT_buffer, *BT_rbuffer; /* buffers for block_transpose */

int DsizeX,DsizeY,DsizeZ;    /* size of the 'matrix' D */
int gridX,gridY,gridZ,Ngrid; /* size of the 'matrix' X */
int gridXY,gridYZ,gridXZ;    /* derived */
int D2sizeX,D2sizeY,D2sizeZ; /* size of the 'matrix' D2 */
int DsizeXY,DsizeYZ,DsizeXZ; /* derived */
int smallY,smallZ;           /* the size of the reduced matrix X */
int local_Nsmall;            /* number of  points of expanded grid per one processor */
int NDcomp;                  /* number of components of D; equals to 6 */

int blockTr=BLOCK;           /* block size for transposeYZ; see fft.h */

#ifdef FFTW3
  /* FFTW3 plans: f - FORWARD; b - BACKWARD */
  fftw_plan planXf,planXb,planYf,planYb,planZf,planZb,planXf_Dm,planYf_Dm,planZf_Dm;
#elif defined(TEMPERTON)
  double *trigsX,*trigsY,*trigsZ,*work;     /* arrays for Temperton FFT */
  int ifaxX[20],ifaxY[20],ifaxZ[20];
#endif

/*============================================================*/

INLINE int index_tDmatrix(int x,int y,int z)
   /* index D matrix to store final result */
{
  if (y>=DsizeY) y=gridY-y;
  if (z>=DsizeZ) z=gridZ-z;

  return(NDcomp*(x*DsizeYZ+z*DsizeY+y));
}

/*============================================================*/

INLINE int index_garbledD(int x,int y,int z, int lengthN)
   /* index D2 matrix after block_transpose */
{
  if (y<0) y+=D2sizeY;
  if (z<0) z+=D2sizeZ;
#ifdef PARALLEL  
  return(((z%lengthN)*D2sizeY+y)*gridX+(z/lengthN)*local_Nx+x%local_Nx);
#else
  return((z*D2sizeY+y)*gridX+x);
#endif
}

/*============================================================*/

INLINE int index_D2matrix(int x,int y,int z, int nnn)
     /* index D2 matrix to store calculate elements */
{
  if (x<0) x+=gridX;
  if (y<0) y+=D2sizeY;
  if (z<0) z+=D2sizeZ;
  return(((z-nnn*local_z0)*D2sizeY+y)*gridX+x);
}

/*============================================================*/

INLINE int index_sliceD2matrix(int y,int z)
    /* index slice of D2 matrix */
{
  if (y<0) y+=gridY;
  if (z<0) z+=gridZ;
  
  return(y*gridZ+z);
}

/*============================================================*/

INLINE int index_slice_zyD2matrix(int y,int z)
   /* index transposed slice of D2 matrix */
{
  return (z*gridY+y);
}

/*============================================================*/

void transposeYZ(int direction)
     /* optimised routine to transpose y and z
        forward: slices -> slices_tr
        backward: slices_tr -> slices  */
{
  int y,z,Y,Z,y1,y2,z1,z2,i,j,y0,z0,Xcomp;
  doublecomplex *t0,*t1,*t2,*t3,*t4,*w0,*w1,*w2,*w3;

  if (direction==FORWARD) { 
    Y=gridY;
    Z=gridZ;
    w0=slices;
    t0=slices_tr-Y;
  }
  else {    /* direction==BACKWARD */
    Y=gridZ;
    Z=gridY;
    w0=slices_tr;
    t0=slices-Y;
  } 

  y1=Y/blockTr;
  y2=Y%blockTr;
  z1=Z/blockTr;
  z2=Z%blockTr;

  for(Xcomp=0;Xcomp<3;Xcomp++) {
    w1=w0+Xcomp*gridYZ;
    t1=t0+Xcomp*gridYZ;
    for(i=0;i<=y1;i++) {
      if (i==y1) y0=y2;
      else y0=blockTr;
      w2=w1;
      t2=t1;
      for(j=0;j<=z1;j++) {
        if (j==z1) z0=z2;
        else z0=blockTr;
        w3=w2;
        t3=t2;
        for (y=0;y<y0;y++) {
          t4=t3+y;
          for (z=0;z<z0;z++) {
            memcpy(t4+=Y,w3+z,sizeof(doublecomplex));
          }
          w3+=Z;
        } 
        w2+=blockTr;
        t2+=blockTr*Y;
      }
      w1+=blockTr*Z;
      t1+=blockTr;
    }
  }
}
/*============================================================*/

void transposeYZ_Dm(doublecomplex *data,doublecomplex *trans)
     /* optimised routine to transpose y and z for Dmatrix: data -> trans */
{
  int y,z,Y,Z,y1,y2,z1,z2,i,j,y0,z0;
  doublecomplex *t1,*t2,*t3,*t4,*w1,*w2,*w3;

  Y=gridY;
  Z=gridZ;

  y1=Y/blockTr;
  y2=Y%blockTr;
  z1=Z/blockTr;
  z2=Z%blockTr;

  w1=data;
  t1=trans-Y;

  for(i=0;i<=y1;i++) {
    if (i==y1) y0=y2;
    else y0=blockTr;
    w2=w1;
    t2=t1;
    for(j=0;j<=z1;j++) {
      if (j==z1) z0=z2;
      else z0=blockTr;
      w3=w2;
      t3=t2;
      for (y=0;y<y0;y++) {
        t4=t3+y;
        for (z=0;z<z0;z++) {
          memcpy(t4+=Y,w3+z,sizeof(doublecomplex));
        }
        w3+=Z;
      } 
      w2+=blockTr;
      t2+=blockTr*Y;
    }
    w1+=blockTr*Z;
    t1+=blockTr;
  }
}

/*============================================================*/

void fftX(int isign)
  /* FFT three components of Xmatrix(x) for all y,z; called from matvec */
{

#ifdef FFTW3
  if (isign==FORWARD) fftw_execute(planXf);
  else fftw_execute(planXb);
#elif defined(TEMPERTON)
  int nn=gridX,z,inc=1,jump=nn,lot=boxY;

  for (z=0;z<3*local_Nz_f;z++)  /* -f */
    cfft99_ ((double *)(Xmatrix+z*gridX*smallY),work,trigsX,ifaxX,&inc,&jump,&nn,&lot,&isign);
#endif
}

/*============================================================*/

void fftY(int isign)
  /* FFT three components of slices_tr(y) for all z; called from matvec */
{
#ifdef FFTW3
  if (isign==FORWARD) fftw_execute(planYf);
  else fftw_execute(planYb);
#elif defined(TEMPERTON)
  int nn=gridY,inc=1,jump=nn,lot=smallZ,j;
  for(j=0;j<6;j++)
    cfft99_ ((double *)(slices_tr+j*gridY*smallZ),work,trigsY,ifaxY,&inc,&jump,&nn,&lot,&isign);
     /* cfft99_ slows down rapidly when lot is big, hence a small loop */
#endif
}

/*============================================================*/

void fftZ(int isign)
  /* FFT three components of slices(z) for all y; called from matvec */
{
#ifdef FFTW3
  if (isign==FORWARD) fftw_execute(planZf);
  else fftw_execute(planZb);
#elif defined(TEMPERTON)
  int nn=gridZ,y,inc=1,jump=nn,lot=boxY,Xcomp;

  for (Xcomp=0;Xcomp<3;Xcomp++)
    cfft99_ ((double *)(slices+gridYZ*Xcomp),work,trigsZ,ifaxZ,&inc,&jump,&nn,&lot,&isign);
#endif
}

/*============================================================*/

void fftX_Dm(int lengthZ)
  /* FFT(forward) D2matrix(x) for all y,z; used for Dmatrix calculation */
{
#ifdef FFTW3
  fftw_execute(planXf_Dm);
#elif defined(TEMPERTON)
  int nn=gridX,z,inc=1,jump=nn,lot=D2sizeY,isign=FORWARD;

  for (z=0;z<lengthZ;z++)
    cfft99_ ((double *)(D2matrix+z*gridX*D2sizeY),work,trigsX,ifaxX,&inc,&jump,&nn,&lot,&isign);
#endif
}

/*============================================================*/

void fftY_Dm(void)
   /* FFT(forward) slice_tr(y) for all z; used for Dmatrix calculation */
{
#ifdef FFTW3
  fftw_execute(planYf_Dm);
#elif defined(TEMPERTON)
  int nn=gridY,inc=1,jump=nn,lot=gridZ,isign=FORWARD;

  cfft99_ ((double *)slice_tr,work,trigsY,ifaxY,&inc,&jump,&nn,&lot,&isign);
#endif
}

/*============================================================*/

void fftZ_Dm(void)
   /* FFT(forward) slice(z) for all y; used for Dmatrix calculation */
{
#ifdef FFTW3
  fftw_execute(planZf_Dm);
#elif defined(TEMPERTON)
  int nn=gridZ,inc=1,jump=nn,lot=gridY,isign=FORWARD;

  cfft99_ ((double *)slice,work,trigsZ,ifaxZ,&inc,&jump,&nn,&lot,&isign);
#endif
}

/*============================================================*/

int fft_fit(int x, int div)
   /* find the first number >=x divisible by 2,3,5 and 7 (if FFTW3)
      only, and divisible by 2 and div */
{
  int y;

  do {
    y=x;
    while (y%2==0) y/=2;
    while (y%3==0) y/=3;
    while (y%5==0) y/=5;
#ifdef FFTW3
    while (y%7==0) y/=7;
#endif
    if (y==1 && x%2==0 && x%div==0) return(x);
    x++;
  } while(true);
}

/*=============================================================*/

void fft_init_beforeD(int lengthZ)
    /* initialize fft before initialization of Dmatrix */
{
#ifdef FFTW3
  planYf_Dm=fftw_plan_many_dft(1,&gridY,gridZ,slice_tr,NULL,1,gridY,slice_tr,NULL,1,gridY,FORWARD,PLAN_FFTW_DM);
  planZf_Dm=fftw_plan_many_dft(1,&gridZ,gridY,slice,NULL,1,gridZ,slice,NULL,1,gridZ,FORWARD,PLAN_FFTW_DM);
  planXf_Dm=fftw_plan_many_dft(1,&gridX,lengthZ*D2sizeY,D2matrix,NULL,1,gridX,D2matrix,NULL,1,gridX,FORWARD,PLAN_FFTW_DM);
#elif defined(TEMPERTON)
  int size,nn;

  /* allocate memory */
  if ((trigsX = (double *) malloc(2*gridX*sizeof(double))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"could not malloc trigsX");
  if ((trigsY = (double *) malloc(2*gridY*sizeof(double))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"could not malloc trigsY");
  if ((trigsZ = (double *) malloc(2*gridZ*sizeof(double))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"could not malloc trigsZ");
  size=MAX(gridX*D2sizeY,3*gridYZ);
  if ((work = (double *) malloc(2*size*sizeof(double))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"could not malloc work");
  /* initialize ifax and trigs */
  nn=gridX;
  cftfax_ (&nn,ifaxX,trigsX);
  nn=gridY;
  cftfax_ (&nn,ifaxY,trigsY);
  nn=gridZ;
  cftfax_ (&nn,ifaxZ,trigsZ);
#endif
}

/*============================================================*/
void fft_init_afterD(void)
    /* second part of fft initialization */
{
#ifdef FFTW3
  int lot;
  fftw_iodim dims,howmany_dims[2];
# ifdef PRECISE_TIMING
  char stat[500];
  SYSTEM_TIME tvp[13];
# endif
  printz("Initializing FFTW3\n");
  fflushz(stdout);
# ifdef PRECISE_TIMING
  GetTime(tvp);
# endif
  lot=3*gridZ;
  planYf=fftw_plan_many_dft(1,&gridY,lot,slices_tr,NULL,1,gridY,slices_tr,NULL,1,gridY,FORWARD,PLAN_FFTW);
# ifdef PRECISE_TIMING
  GetTime(tvp+1);
# endif
  planYb=fftw_plan_many_dft(1,&gridY,lot,slices_tr,NULL,1,gridY,slices_tr,NULL,1,gridY,BACKWARD,PLAN_FFTW);
# ifdef PRECISE_TIMING
  GetTime(tvp+2);
# endif
  dims.n=gridZ;
  dims.is=dims.os=1;
  howmany_dims[0].n=3;
  howmany_dims[0].is=howmany_dims[0].os=gridZ*gridY;
  howmany_dims[1].n=boxY;
  howmany_dims[1].is=howmany_dims[1].os=gridZ;
  planZf=fftw_plan_guru_dft(1,&dims,2,howmany_dims,slices,slices,FORWARD,PLAN_FFTW);
# ifdef PRECISE_TIMING  
  GetTime(tvp+3);
# endif
  planZb=fftw_plan_guru_dft(1,&dims,2,howmany_dims,slices,slices,BACKWARD,PLAN_FFTW);
# ifdef PRECISE_TIMING  
  GetTime(tvp+4);
# endif
  dims.n=gridX;
  dims.is=dims.os=1;
  howmany_dims[0].n=3*local_Nz_f;
  howmany_dims[0].is=howmany_dims[0].os=smallY*gridX; 
  howmany_dims[1].n=boxY;
  howmany_dims[1].is=howmany_dims[1].os=gridX;
  planXf=fftw_plan_guru_dft(1,&dims,2,howmany_dims,Xmatrix,Xmatrix,FORWARD,PLAN_FFTW);
# ifdef PRECISE_TIMING  
  GetTime(tvp+5);
# endif
  planXb=fftw_plan_guru_dft(1,&dims,2,howmany_dims,Xmatrix,Xmatrix,BACKWARD,PLAN_FFTW);
# ifdef PRECISE_TIMING
  GetTime(tvp+6);
  /* print precise timing of FFT planning */
  SetTimerFreq();
  sprintz(stat,
         "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"\
         "         FFTW3 planning       \n"\
         "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"\
         "Yf = %4.4f  Total = %4.4f\n"\
         "Yb = %4.4f\n"\
         "Zf = %4.4f\n"\
         "Zb = %4.4f\n"\
         "Xf = %4.4f\n"\
         "Xb = %4.4f\n\n",
         DiffSec(tvp,tvp+1),DiffSec(tvp,tvp+6),DiffSec(tvp+1,tvp+2),DiffSec(tvp+2,tvp+3),
         DiffSec(tvp+3,tvp+4),DiffSec(tvp+4,tvp+5),DiffSec(tvp+5,tvp+6));
  printz(stat);
  fprintz(logfile,stat);
# endif
  /* destroy old plans */
  fftw_destroy_plan(planXf_Dm);
  fftw_destroy_plan(planYf_Dm);
  fftw_destroy_plan(planZf_Dm);
#endif
}

/**************************************************************************/

void CalcInterTerm(int i, int j, int k, int mu, int nu, doublecomplex result)
    /* calculates interaction term between two dipoles; given integer distance vector {i,j,k}
        (in units of d), and component indices mu,nu */
{
  extern double *tab1,*tab2,*tab3,*tab4,*tab5,*tab6,*tab7,*tab8,*tab9,*tab10; /* tables of integrals */
  extern double prop[3];
  extern int IntRelation;
  extern double kd;
  extern doublecomplex ref_index[MAXNMAT];
  extern int **tab_index;

  double rr, rtemp[3], qvec[3], q2[3], invr, invr3, qavec[3], av[3];
  double rr2, kr, kr2, kr3, kd2, q4, rn;
  double temp, qmunu, qa, qamunu, invrn, invrn2, invrn3, invrn4, dmunu;
  doublecomplex expval, br, br1, m, m2, Gf1, Gm0, Gm1, Gc1, Gc2;
  int ind0, ind1, ind2, ind2m, ind3, ind4, indmunu;
  int sigV[3], ic, sig, ivec[3], ord[3], invord[3];
  double t3q, t3a, t4q, t4a, t5tr, t5aa, t6tr, t6aa;
  int pr;

  if (i==0 && j==0 && k==0) { /* self interaction */
    result[re]=result[im]=0.0;
    return;
  }

/*  pr=(i==1 && j==1 && k==1);
  if (pr) printz("%d,%d: ",mu,nu);    /* for debugging */

  rtemp[0]=i*gridspace;
  rtemp[1]=j*gridspace;
  rtemp[2]=k*gridspace;

  rr2 = DotProd(rtemp,rtemp);
  rr = sqrt(rr2);
  invr = 1/rr;
  invr3 = invr*invr*invr;
  MultScal(invr,rtemp,qvec);
  kr = WaveNum * rr;
  kr2 = kr*kr;
  qmunu=qvec[mu]*qvec[nu];

  br[re]=(3-kr2)*qmunu;
  br[im]=-3*kr*qmunu;
  if(mu==nu) {                       /* br=delta[mu,nu]*(-1+ikr+kr^2)  */
    br[re]+=kr2-1;                     /*    -qmunu*(-3+3ikr+kr^2)       */
    br[im]+=kr;
  }

  cExp(kr,expval);
  cMultReal(invr3,expval,expval);      /* expval=Exp(ikr)/rr^3 */
  cMult(br,expval,result);             /* result=Gp            */

  if (IntRelation==SOrd) {
    kd2=kd*kd;
    kr3=kr2*kr;
    rn=rr/gridspace;              /* normalized r */
    memcpy(m,ref_index[0],sizeof(doublecomplex));       /* only one refractive index can be used for FFT-compatible algorithm !!! */
    cSquare(m,m2);
    qa=DotProd(qvec,prop);

    qamunu=qvec[mu]*prop[nu];          /* qamunu=qvec[mu]*prop[nu] + qvec[nu]*prop[mu] */
    if (mu==nu) qamunu*=2;
    else qamunu+=qvec[nu]*prop[mu];

    if (kr*rn < 1) {       /* G close */
      memcpy(av,prop,3*sizeof(double));   /* av is copy of propagation vector */
      ivec[0]=i;
      ivec[1]=j;
      ivec[2]=k;
      /* transformation of negative coordinates */
      for (ic=0;ic<3;ic++) {
        if (ivec[ic]<0) {
          sigV[ic]=-1;
          av[ic]*=-1;
          qvec[ic]*=-1;
          ivec[ic]*=-1;
        }
        else sigV[ic]=1;
      }
      i=ivec[0];
      j=ivec[1];
      k=ivec[2];
      sig=sigV[mu]*sigV[nu];           /* sign of some terms below */
      /* transformation to case i>=j>=k>=0 */
        /* building of ord */
      if (i>=j) {
        if (i>=k) {
          ord[0]=0;
          if (j>=k) {
            ord[1]=1;
            ord[2]=2;
          }
          else {
            ord[1]=2;
            ord[2]=1;
          }
        }
        else {
          ord[0]=2;
          ord[1]=0;
          ord[2]=1;          /* ord[x] is x-th largest coordinate (0-th - the largest) */
        }
      }
      else {
        if (i>=k) {
          ord[0]=1;
          ord[1]=0;
          ord[2]=2;
        }
        else {
          ord[2]=0;
          if (j>=k) {
            ord[0]=1;
            ord[1]=2;
          }
          else {
            ord[0]=2;
            ord[1]=1;
          }
        }
      }
      /* change parameters according to coordinate transforms */
      permutate(qvec,ord);
      permutate(av,ord);
      permutate_i(ivec,ord);
      i=ivec[0];
      j=ivec[1];
      k=ivec[2];
      /* compute inverse permutation */
      memcpy(invord,ord,3*sizeof(int));
      permutate_i(invord,ord);
      if (invord[0]==0 && invord[1]==1 && invord[2]==2) memcpy(invord,ord,3*sizeof(int));
      /* compute transformed indices mu and nu */
      mu=invord[mu];
      nu=invord[nu];
      /* indexes for tables of different dimensions */
      indmunu=mu+nu;
      if (mu==2 || nu==2) indmunu++;  /* indmunu is a number of component[mu,nu] in symmetric matrix */

      ind0=tab_index[i][j]+k;
      ind1=3*ind0;
      ind2m=6*ind0;
      ind2=ind2m+indmunu;
      ind3=3*ind2;
      ind4=6*ind2;
      /* computing several quantities with table integrals */
      t3q=DotProd(qvec,tab3+ind1);
      t3a=DotProd(av,tab3+ind1);
      t4q=DotProd(qvec,tab4+ind3);
      t4a=DotProd(av,tab4+ind3);
      t5tr=TrSym(tab5+ind2m);
      t5aa=QuadForm(tab5+ind2m,av);
      t6tr=TrSym(tab6+ind4);
      t6aa=QuadForm(tab6+ind4,av);
      /* computing Gc0 */
      temp=kr/24;                      /* temp = kr/12 */

      br[re]=sig*(3*(tab10[ind2]/2+tab8[ind2])-2*t4q-t6tr)+temp*qmunu*kr;
      br[im]=3*temp*qmunu;                                   /* br=delta[mu,nu]*(-I7-I9/2-kr*(i+kr)/24+2*t3q+t5tr) */
      if (mu==nu) {                                        /*   -(-3I8[mu,nu]-3I10[mu,nu]/2-qmunu*kr*(i+kr)/24   */
        br[re]+=2*t3q+t5tr-temp*kr-tab9[ind0]/2-tab7[ind0];  /*     +2*t4q+t6tr)                                   */
        br[im]-=temp;
      }

      cMultReal(kd2,br,br);          /* br*=kd^2 */

      br[re]+=sig*tab2[ind2]*(3-kr2);   /* br+=I1*delta[mu,nu]*(-1+ikr+kr^2) */
      br[im]-=sig*tab2[ind2]*3*kr;      /*     -sig*I2[mu,nu]*(-3+3ikr+kr^2) */
      if (mu==nu) {
        br[re]+=tab1[ind0]*(kr2-1);
        br[im]+=tab1[ind0]*kr;
      }

      cMult(expval,br,result);        /* Gc0=expval*br */
      /* computing Gc1 */
      br[re]=6*qmunu;
      br[im]=-kr*qmunu;      /* br=(kd*kr/24)*(qa*(delta[mu,nu]*(-2+ikr)-qmunu*(-6+ikr)) */
      if (mu==nu) {        /*                -qamunu)                                  */
        br[re]-=2;
        br[im]+=kr;
      }
      cMultReal(qa,br,br);
      br[re]-=qamunu;
      cMultReal(2*temp*kd,br,br);

      br1[re]=3*sig*t4a;
      br1[im]=-kr*br1[re];
      if (mu==nu) {        /*  br1=(d/r)*(delta[mu,nu]*t3h*(-1+ikr)-sig*t4h*(-3+3ikr)) */
        br1[re]-=t3a;
        br1[im]+=t3a*kr;
      }
      cMultReal(1/rn,br1,br1);

      cAdd(br,br1,Gc1);
      cMultSelf(Gc1,m);
      cMultReal(kd,Gc1,Gc1);
      cMultSelf(Gc1,expval);
      cMult_i(Gc1);                   /* Gc1=expval*i*m*kd*(br1+br) */
      /* computing Gc2 */
      br[re]=-kr*qmunu;
      br[im]=-3*qmunu;      /* br=delta[mu,nu]*t5aa-3*sig*t6aa                 */
      if (mu==nu) {       /*    -(kr/12)*(delta[mu,nu]*(i+kr)-qmunu*(3i+kr)) */
        br[re]+=kr;
        br[im]+=1;
      }
      cMultReal(-2*temp,br,br);
      br[re]-=3*sig*t6aa;
      if (mu==nu) br[re]+=t5aa;

      cMult(m2,br,Gc2);
      cMultReal(kd2/2,Gc2,Gc2);     /* Gc2=expval*(kd^2/2)*m^2*br */
      cMultSelf(Gc2,expval);

      cAdd(Gc1,Gc2,Gc1);
      cAdd(Gc1,result,result);        /* result = Gc0 + Gc1 + Gc2 */
    }
    else {              /* G median and G far */
      temp=kd2/24;             /* temp=kd^2/24;        */
      br[re]=1-(1+m2[re])*temp;    /* br=1-(1+m^2)*kd^2/24 */
      br[im]=-m2[im]*temp;
      cMultSelf(result,br);   /* result = Gp*br */

      br[re]=(6*kr2-15)*qmunu;
      br[im]=(15*kr-kr3)*qmunu;
      if(mu==nu) {                       /* br={delta[mu,nu]*(3-3ikr-2kr^2+ikr^3)  */
        br[re]+=3-2*kr2;                   /*    -qmunu*(15-15ikr-6kr^2+ikr^3)}*qa   */
        br[im]+=kr3-3*kr;                  /*    +qamunu*(3-3ikr-kr^2)               */
      }
      cMultReal(qa,br,br);
      br[re]+=(3-kr2)*qamunu;
      br[im]-=3*kr*qamunu;

      temp*=2/kr;                /* temp = kd^2/(12*kr) */
      cMult(m,br,Gf1);
      cMultReal(temp,Gf1,Gf1);
      cMultSelf(Gf1,expval);
      cMult_i(Gf1);                     /* Gf1=expval*i*m*temp*br */
      cAdd(Gf1,result,result);          /* result = Gf  */

      if (kr < 1) {                      /* G median */
        vMult(qvec,qvec,q2);
        q4=DotProd(q2,q2);
        invrn=1/rn;
        invrn2=invrn*invrn;
        invrn3=invrn2*invrn;
        invrn4=invrn2*invrn2;

        temp=qmunu*(33*q4-7-12*(q2[mu]+q2[nu]));
        if (mu == nu) temp+=(1-3*q4+4*q2[mu]);
        temp*=7*invrn4/64;
        br[re]=-1;
        br[im]=kr;
        cMultReal(temp,br,Gm0);
        cMultSelf(Gm0,expval);                     /* Gm0=expval*br*temp */

        vMult(qvec,prop,qavec);
        if (mu == nu) dmunu=1;
        else dmunu=0;
        temp=3*qa*(dmunu-7*qmunu)+6*dmunu*qvec[mu]*prop[mu]-7*(dmunu-9*qmunu)*DotProd(qavec,q2)+
          3*(prop[mu]*qvec[nu]*(1-7*q2[mu])+prop[nu]*qvec[mu]*(1-7*q2[nu]));
        temp*=kd*invrn3/48;
        cMultReal(temp,m,Gm1);
        cMult_i(Gm1);
        cMultSelf(Gm1,expval);    /* Gm1=expval*i*m*temp */

        cAdd(Gm0,Gm1,Gm0);
        cAdd(Gm0,result,result);   /* result=Gf+Gm0+Gm1 */
      }

    }
  }
  /* if (pr) printz("%d,%d: %f+%fi\n",mu,nu,result[re],result[im]); */
}

/*============================================================*/

void init_Dmatrix(void)
     /* initialises the matrix D. D[i][j][k]=A[i1-i2][j1-j2][k1-k2]
	The routine is called only once, so needs not to be very fast */
{
  int y,z,i,j,k,ind,kcor,index,mem,Dcomp,indexfrom,indexto,D2sizeTot;
  double invNgrid;

  int nnn; /* multiplier used for reduced_FFT or not reduced; 1 or 2 */
  int jstart, kstart;
  int lengthN, bufsize;

  int mu, nu;     /* indices for interaction term */
  clock_t start,time1;
#ifdef PRECISE_TIMING
  /* precise timing of the Dmatrix computation */
  char stat[500];
  SYSTEM_TIME tvp[13];
  SYSTEM_TIME Timing_fftX,Timing_fftY,Timing_fftZ,Timing_ar1,Timing_ar2,Timing_ar3,
              Timing_BT,Timing_TYZ,Timing_beg;
  double t_fftX,t_fftY,t_fftZ,t_ar1,t_ar2,t_ar3,
         t_TYZ,t_beg,t_Arithm,t_FFT,t_BT;

  InitTime(&Timing_fftX);
  InitTime(&Timing_fftY);
  InitTime(&Timing_fftZ);
  InitTime(&Timing_ar1);
  InitTime(&Timing_ar2);
  InitTime(&Timing_ar3);
  InitTime(&Timing_BT);
  InitTime(&Timing_TYZ);
  GetTime(tvp);
#endif
  start=clock();
  /* initialize sizes of D and D2 matrices */
  DsizeX=gridX;
  D2sizeX=gridX;
  if (reduced_FFT) {
    D2sizeY=gridY/2;
    D2sizeZ=gridZ/2;
    DsizeY=gridY/2+1;
    DsizeZ=gridZ/2+1;
    nnn=1;
    jstart=0;
    kstart=0;
  }
  else {
    D2sizeY=DsizeY=gridY;
    D2sizeZ=DsizeZ=gridZ;
    nnn=2;
    jstart=1-boxY;
    kstart=1-boxZ;
  }
  /* auxiliary parameters */
  lengthN=nnn*local_Nz_f;
  DsizeXY=DsizeX*DsizeY;
  DsizeYZ=DsizeY*DsizeZ;
  DsizeXZ=DsizeX*DsizeZ;
  Ngrid=gridX*gridY*gridZ;
  invNgrid=1.0/Ngrid;
  local_Nsmall=gridX*gridY*gridZ/(4*nprocs); /* size of X vector (for 1 component) */
  NDcomp=6;
  /* calculate size of matvec matrices (X,D,slices,slices_tr) and BT buffers (if parallel) */
  mem=sizeof(doublecomplex)*(3*(local_Nsmall+2*gridYZ)+NDcomp*local_Nx*DsizeY*DsizeZ);
#ifdef PARALLEL
  mem+=12*smallY*local_Nz_f*local_Nx*sizeof(double);
#endif
  /* printout some information */
  fprintz(logfile,"The 3d grid is:%ix%ix%i\n",gridX,gridY,gridZ);
#ifdef PARALLEL
  printz("Memory usage for MatVec matrices (per processor):%.1f Mb\n",mem/MBYTE);
  fprintz(logfile,"Memory usage for MatVec matrices (per processor):%.1f Mb\n",mem/MBYTE);
#else
  printz("Memory usage for MatVec matrices:%.1f Mb\n",mem/MBYTE);
  fprintz(logfile,"Memory usage for MatVec matrices:%.1f Mb\n",mem/MBYTE);
#endif
  fflushz(logfile);
  memory+=mem;
  if (prognose==true) return;
  /* allocate memory for Dmatrix */
  if ((Dmatrix = dCvector(NDcomp*local_Nx*DsizeY*DsizeZ)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Dmatrix");
  /* allocate memory for D2matrix components */
  D2sizeTot=nnn*local_Nz_f*D2sizeY*D2sizeX;
  if ((D2matrix = dCvector(D2sizeTot)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc D2matrix");
  if ((slice = dCvector(gridYZ)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc slice");
  if ((slice_tr = dCvector(gridYZ)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc slice_tr");
  /* actually allocation of Xmatrix, slices, slices_tr is below;
     after freeing of Dmatrix and its slice */
#ifdef PARALLEL
  /* allocate buffer for block_transpose_Dm */
  bufsize = 2*lengthN*D2sizeY*local_Nx;
  if ((BT_buffer = dvector(0,bufsize-1)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc BT_buffer");
  if ((BT_rbuffer = dvector(0,bufsize-1)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc BT_rbuffer");
#endif
  D("init FFT (1st part)");
  fft_init_beforeD(lengthN);

#ifdef PRECISE_TIMING
  GetTime(tvp+1);
  elapsed(tvp,tvp+1,&Timing_beg);
#endif
  printz("Calculating Dmatrix");
  fflushz(stdout);

  for(Dcomp=0;Dcomp<NDcomp;Dcomp++) {   /* main cycle over components of Dmatrix */
#ifdef PRECISE_TIMING
    GetTime(tvp+2);
#endif
    switch((char) Dcomp) {    /* determine mu,nu */
      case 0: {
        mu=0;
        nu=0;
	break;
      }
      case 1: {
        mu=0;
        nu=1;
	break;
      }
      case 2: {
        mu=0;
        nu=2;
	break;
      }
      case 3: {
        mu=1;
        nu=1;
	break;
      }
      case 4: {
        mu=1;
        nu=2;
	break;
      }
      case 5: {
        mu=2;
        nu=2;
	break;
      }
    } /* end of switch */

    /* fill D2matrix with 0.0 */
    for (i=0;i<D2sizeTot;i++) D2matrix[i][re]=D2matrix[i][im]=0.0;

    /* fill D (F'i-j) */
    for(k=nnn*local_z0;k<nnn*local_z1_f;k++) {
      if (k>smallZ) kcor=k-gridZ;
      else kcor=k;
      for (j=jstart;j<boxY;j++) for (i=1-boxX;i<boxX;i++) {
        index=index_D2matrix(i,j,k,nnn);
        CalcInterTerm(i,j,kcor,mu,nu,D2matrix[index]);   /* calculate F[mu][nu] */
      }
    } /* end of i,j,k loop */
#ifdef PRECISE_TIMING
    GetTime(tvp+3);
    ElapsedInc(tvp+2,tvp+3,&Timing_ar1);
#endif
    fftX_Dm(lengthN); /* fftX D2matrix */
#ifdef PRECISE_TIMING
    GetTime(tvp+4);
    ElapsedInc(tvp+3,tvp+4,&Timing_fftX);
#endif
    block_transpose_Dm(D2matrix,D2sizeY,lengthN);
#ifdef PRECISE_TIMING
    GetTime(tvp+5);
    ElapsedInc(tvp+4,tvp+5,&Timing_BT);
#endif

    for(i=local_x0;i<local_x1;i++) {
#ifdef PRECISE_TIMING
      GetTime(tvp+6);
#endif
      for (ind=0;ind<gridYZ;ind++) slice[ind][re]=slice[ind][im]=0.0;  /* fill slice with 0.0 */

      for(j=jstart;j<boxY;j++) for(k=kstart;k<boxZ;k++) {
	indexfrom=index_garbledD(i,j,k,lengthN);
	indexto=index_sliceD2matrix(j,k);
	memcpy(slice[indexto],D2matrix[indexfrom],sizeof(doublecomplex));
      }

      if (reduced_FFT) {
        for(j=1;j<boxY;j++) for(k=0;k<boxZ;k++) {
          /* mirror along y */
          indexfrom=index_sliceD2matrix(j,k);
          indexto=index_sliceD2matrix(-j,k);
          if (Dcomp==1 || Dcomp==4) cInvSign2(slice[indexfrom],slice[indexto]);
          else memcpy(slice[indexto],slice[indexfrom],sizeof(doublecomplex));
        }
        for(j=1-boxY;j<boxY;j++) for(k=1;k<boxZ;k++) {
          /* mirror along z */
          indexfrom=index_sliceD2matrix(j,k);
          indexto=index_sliceD2matrix(j,-k);
          if (Dcomp==2 || Dcomp==4) cInvSign2(slice[indexfrom],slice[indexto]);
          else memcpy(slice[indexto],slice[indexfrom],sizeof(doublecomplex));
        }
      }
#ifdef PRECISE_TIMING
      GetTime(tvp+7);
      ElapsedInc(tvp+6,tvp+7,&Timing_ar2);
#endif
      fftZ_Dm();  /* fftZ slice */
#ifdef PRECISE_TIMING
      GetTime(tvp+8);
      ElapsedInc(tvp+7,tvp+8,&Timing_fftZ);
#endif
      transposeYZ_Dm(slice,slice_tr);
#ifdef PRECISE_TIMING
      GetTime(tvp+9);
      ElapsedInc(tvp+8,tvp+9,&Timing_TYZ);
#endif
      fftY_Dm();  /* fftY slice_tr */
#ifdef PRECISE_TIMING
      GetTime(tvp+10);
      ElapsedInc(tvp+9,tvp+10,&Timing_fftY);
#endif
      for(z=0;z<DsizeZ;z++) for(y=0;y<DsizeY;y++) {
	indexto=index_tDmatrix(i-local_x0,y,z)+Dcomp;
	indexfrom=index_slice_zyD2matrix(y,z);
	Dmatrix[indexto][re]=-slice_tr[indexfrom][re]*invNgrid;
	Dmatrix[indexto][im]=-slice_tr[indexfrom][im]*invNgrid;
      }
#ifdef PRECISE_TIMING
      GetTime(tvp+11);
      ElapsedInc(tvp+10,tvp+11,&Timing_ar3);
#endif
    } /* end slice X */
    printz(".");
    fflushz(stdout);
  } /* end of Dcomp */
  /* free vectors used for computation of Dmatrix */
  free_dCvector(D2matrix);
  free_dCvector(slice);
  free_dCvector(slice_tr);
#ifdef PARALLEL
  /* deallocate buffers for block_transpose_Dm */
  free_dvector(BT_buffer,0);
  free_dvector(BT_rbuffer,0);
  /* allocate buffers for block_transpose */
  bufsize = 6*smallY*local_Nz_f*local_Nx; /* in doubles */
  if ((BT_buffer = dvector(0,bufsize-1)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc BT_buffer");
  if ((BT_rbuffer = dvector(0,bufsize-1)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc BT_rbuffer");
#endif
  /* allocate memory for Xmatrix, slices and slices_tr - used in matvec */
  if ((Xmatrix = dCvector(3*local_Nsmall)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Xmatrix");
  if ((slices = dCvector(3*gridYZ)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc slices");
  if ((slices_tr = dCvector(3*gridYZ)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc slices_tr");

  printz("\n");
  time1=clock();
  Timing_Dm_Init=time1-start;
#ifdef PRECISE_TIMING
  GetTime(tvp+12);
  /* analyze and print precise timing information */
  SetTimerFreq();
  t_beg=TimerToSec(&Timing_beg);
  t_ar1=TimerToSec(&Timing_ar1);
  t_ar2=TimerToSec(&Timing_ar2);
  t_ar3=TimerToSec(&Timing_ar3);
  t_fftX=TimerToSec(&Timing_fftX);
  t_fftY=TimerToSec(&Timing_fftY);
  t_fftZ=TimerToSec(&Timing_fftZ);
  t_TYZ=TimerToSec(&Timing_TYZ);
  t_BT=TimerToSec(&Timing_BT);
  t_Arithm=t_beg+t_ar1+t_ar2+t_ar3+t_TYZ;
  t_FFT=t_fftX+t_fftY+t_fftZ;

  sprintz(stat,
         "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"\
         "            Init Dmatrix timing            \n"\
         "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"\
         "Begin  = %4.4f    Arithmetics = %4.4f\n"\
         "Arith1 = %4.4f    FFT         = %4.4f\n"\
         "FFTX   = %4.4f    Comm        = %4.4f\n"\
         "BT     = %4.4f\n"\
         "Arith2 = %4.4f          Total = %4.4f\n"\
         "FFTZ   = %4.4f\n"\
         "TYZ    = %4.4f\n"\
         "FFTY   = %4.4f\n"\
         "Arith3 = %4.4f\n\n",
         t_beg,t_Arithm,t_ar1,t_FFT,t_fftX,t_BT,t_BT,
         t_ar2,DiffSec(tvp,tvp+12),t_fftZ,t_TYZ,t_fftY,t_ar3);
  printz(stat);
  fprintz(logfile,stat);
#endif

  fft_init_afterD();

  Timing_FFT_Init = clock()-time1;
}

/*============================================================*/

void free_FFT_Dmat(void)
   /* free all vectors that were allocated in fft.c
       (all used for FFT and MatVec) */
{
  free_dCvector(Dmatrix);
  free_dCvector(Xmatrix);
  free_dCvector(slices);
  free_dCvector(slices_tr);
#ifdef PARALLEL
  free_dvector(BT_buffer,0);
  free_dvector(BT_rbuffer,0);
#endif
#ifdef FFTW3
  fftw_destroy_plan(planXf);
  fftw_destroy_plan(planXb);
  fftw_destroy_plan(planYf);
  fftw_destroy_plan(planYb);
  fftw_destroy_plan(planZf);
  fftw_destroy_plan(planZb);
#elif defined(TEMPERTON)
  free(work);
  free(trigsX);
  free(trigsY);
  free(trigsZ);
#endif
}
