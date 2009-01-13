/* FILE: matvec.c
 * AUTH: Michel Grimminck
 * DESCR: calculate local matrix vector product of decomposed interaction
 *        matrix with rk or pk, using a FFT based convolution algorithm
 *
 *        Currently is developed by Maxim Yurkin
 */

#include <stdio.h>
#include <time.h>
#include <math.h>
#include "cmplx.h"
#include "const.h"
#include "comm.h"
#include "debug.h"
#include "memory.h"
#include "fft.h"
#include "timing.h"
#include "linalg.h"

extern doublecomplex *Dmatrix;
extern int smallY,smallZ;
extern int gridX,gridY,gridZ,Ngrid;
extern int gridXY,gridXZ,gridYZ;
extern int boxX,boxY,boxZ;
extern int DsizeY,DsizeZ,DsizeYZ,NDcomp;

/*============================================================*/

INLINE int index_slice_zy(int y,int z)
{
  return (z*gridY+y);
}

/*============================================================*/

INLINE int index_slice_yz(int y,int z)
{
  return(y*gridZ+z);
}

/*============================================================*/

INLINE int index_garbledX(int x,int y,int z)
{
#ifdef PARALLEL
  return(((z%local_Nz_f)*smallY+y)*gridX+(z/local_Nz_f)*local_Nx+x%local_Nx);
#else
  return((z*smallY+y)*gridX+x);
#endif
}

/*============================================================*/

INLINE int index_Xmatrix(int x,int y,int z)
{
  return((z*smallY+y)*gridX+x);
}

/*============================================================*/

INLINE int index_Dmatrix_mv(int x,int y,int z, int transposed)
{
  if (transposed) {   /* used only for SOrd */
    if (x>0) x=gridX-x;
    if (y>0) y=gridY-y;
    if (z>0) z=gridZ-z;
  }
  else {
    if (y>=DsizeY) y=gridY-y;
    if (z>=DsizeZ) z=gridZ-z;
  }
  
  return(NDcomp*(x*DsizeYZ+z*DsizeY+y));
}

/*============================================================*/

/* This function implements both MatVec_nim and MatVecAndInp_nim.
 * the difference is that when we want to calculate the inproduct
 * as well, we pass 'inprod' as a non-null pointer. if 'inprod' is
 * a NULL, we don't calculate it.
 */
void MatVec  (doublecomplex *argvec,     /* the argument vector */
              doublecomplex *resultvec,  /* the result vector */
              double *inprod,            /* the result inner product */
              int her)                   /* 0 for non-hermitic, 1 for hermetic */
{
  int i, j;

  doublecomplex fmat[6],xvec[3],yvec[3];
  doublecomplex temp;
  int index,x,y,z,Xcomp,ipr;
  char mat;

  extern doublecomplex *slices,*slices_tr,*Xmatrix;
  extern int reduced_FFT;
  extern int local_Nsmall;
  extern double gridspace;
  extern short int *position;
  extern doublecomplex cc_sqrt[MAXNMAT][3];
  extern char *material;
  extern clock_t Timing_OneIterComm;
  extern FILE *logfile;

  clock_t tstart;
#ifdef PRECISE_TIMING  
  char stat[500];
  SYSTEM_TIME tvp[18];
  SYSTEM_TIME Timing_FFTXf,Timing_FFTYf,Timing_FFTZf,Timing_FFTXb,Timing_FFTYb,Timing_FFTZb,
              Timing_Mult1,Timing_Mult2,Timing_Mult3,Timing_Mult4,Timing_Mult5,
              Timing_BTf,Timing_BTb,Timing_TYZf,Timing_TYZb,Timing_ipr;
  double t_FFTXf,t_FFTYf,t_FFTZf,t_FFTXb,t_FFTYb,t_FFTZb,
         t_Mult1,t_Mult2,t_Mult3,t_Mult4,t_Mult5,t_ipr,
         t_BTf,t_BTb,t_TYZf,t_TYZb,t_Arithm,t_FFT,t_Comm;

#endif
  static int mvcounter=0;
  int transposed;

  /* A = I + S.D.S
   * S = sqrt(C)
   * A*x = x + S.D.(C.x)
   * A(H)*x = x + (S(T).D(T).S(T).x(*))(*)
   * C,S - diagonal => symmetric
   * D - symmetric (except for SOrd)
   *
   * D.x=F(-1)(F(D).F(X))
   * F(D) is just a vector
   *
   * SOrd: F(D(T)) (k) =  F(D) (-k)
   *       k - vector index
   *
   *   For (her) three additional operations of nConj are used. Should not be a problem,
   *     but can be avoided by a more complex code.
   */

  transposed=(!reduced_FFT) && her;
  if (inprod) ipr=true;
  else ipr=false;

#ifdef PRECISE_TIMING
  InitTime(&Timing_FFTYf);
  InitTime(&Timing_FFTZf);
  InitTime(&Timing_FFTYb);
  InitTime(&Timing_FFTZb);
  InitTime(&Timing_Mult2);
  InitTime(&Timing_Mult3);
  InitTime(&Timing_Mult4);
  InitTime(&Timing_TYZf);
  InitTime(&Timing_TYZb);
  GetTime(tvp);
#endif
  /* FFT_matvec code */
  if (ipr) *inprod = 0.0;

  /* fill Xmatrix with 0.0 */
  for (i=0;i<3*local_Nsmall;i++) Xmatrix[i][re]=Xmatrix[i][im]=0.0;

  /* transform from coordinates to grid and multiply with coupling constant */
  if (her) nConj(argvec); /* conjugated back afterwards */

  for (i=0;i<local_nvoid_Ndip;i++) {
    /* fill grid with E*sqrt_cc */
    j=3*i;
    mat=material[i];
    index=index_Xmatrix(position[j],position[j+1],position[j+2]-local_z0);
    for (Xcomp=0;Xcomp<3;Xcomp++)
      cMult(cc_sqrt[mat][Xcomp],argvec[j+Xcomp],Xmatrix[index+Xcomp*local_Nsmall]);  /* Xmat=cc_sqrt*argvec */
  }
#ifdef PRECISE_TIMING
  GetTime(tvp+1);
  elapsed(tvp,tvp+1,&Timing_Mult1);
#endif
  /* FFT X */
  fftX(FORWARD);    /* fftX Xmatrix */
#ifdef PRECISE_TIMING
  GetTime(tvp+2);
  elapsed(tvp+1,tvp+2,&Timing_FFTXf);
#endif
  block_transpose(Xmatrix);
#ifdef PRECISE_TIMING
  GetTime(tvp+3);
  elapsed(tvp+2,tvp+3,&Timing_BTf);
#endif
  /* following is done by slices */
  for(x=local_x0;x<local_x1;x++) {
#ifdef PRECISE_TIMING
    GetTime(tvp+4);
#endif
    /* clear slice */
    for(i=0;i<3*gridYZ;i++) slices[i][re]=slices[i][im]=0.0;

    /* fill slices with values from Xmatrix */
    for(y=0;y<boxY;y++) for(z=0;z<boxZ;z++) {
      i=index_slice_yz(y,z);
      j=index_garbledX(x,y,z);
      for (Xcomp=0;Xcomp<3;Xcomp++)
        memcpy(slices[i+Xcomp*gridYZ],Xmatrix[j+Xcomp*local_Nsmall],sizeof(doublecomplex));
    }
#ifdef PRECISE_TIMING
    GetTime(tvp+5);
    ElapsedInc(tvp+4,tvp+5,&Timing_Mult2);
#endif
    /* fftZ&Y */
    fftZ(FORWARD);   /* fftZ slices */
#ifdef PRECISE_TIMING
    GetTime(tvp+6);
    ElapsedInc(tvp+5,tvp+6,&Timing_FFTZf);
#endif
    transposeYZ(FORWARD);
#ifdef PRECISE_TIMING
    GetTime(tvp+7);
    ElapsedInc(tvp+6,tvp+7,&Timing_TYZf);
#endif
    fftY(FORWARD);   /* fftY slices_tr */
#ifdef PRECISE_TIMING
    GetTime(tvp+8);
    ElapsedInc(tvp+7,tvp+8,&Timing_FFTYf);
#endif
    /* do the product D~*X~ */
    for(z=0;z<gridZ;z++) for(y=0;y<gridY;y++) {
      i=index_slice_zy(y,z);
      for (Xcomp=0;Xcomp<3;Xcomp++)
        memcpy(xvec[Xcomp],slices_tr[i+Xcomp*gridYZ],sizeof(doublecomplex));

      j=index_Dmatrix_mv(x-local_x0,y,z,transposed);
      memcpy(fmat,Dmatrix[j],6*sizeof(doublecomplex));
      if (reduced_FFT) {
        if (y>smallY) {
          cInvSign(fmat[1]);                  /* fmat[1]*=-1 */
          if (z>smallZ) cInvSign(fmat[2]);    /* fmat[2]*=-1 */
          else cInvSign(fmat[4]);             /* fmat[4]*=-1 */
        }
        else if (z>smallZ) {
          cInvSign(fmat[2]);     /* fmat[2]*=-1 */
          cInvSign(fmat[4]);     /* fmat[4]*=-1 */
        }
      }
      csymMatrVec(fmat,xvec,yvec);    /* yvec=fmat*xvec */
      for (Xcomp=0;Xcomp<3;Xcomp++)
        memcpy(slices_tr[i+Xcomp*gridYZ],yvec[Xcomp],sizeof(doublecomplex));
    }
#ifdef PRECISE_TIMING
    GetTime(tvp+9);
    ElapsedInc(tvp+8,tvp+9,&Timing_Mult3);
#endif
    /* fft_invY&Z */
    fftY(BACKWARD);   /* fftY slices_tr */
#ifdef PRECISE_TIMING
    GetTime(tvp+10);
    ElapsedInc(tvp+9,tvp+10,&Timing_FFTYb);
#endif
    transposeYZ(BACKWARD);
#ifdef PRECISE_TIMING
    GetTime(tvp+11);
    ElapsedInc(tvp+10,tvp+11,&Timing_TYZb);
#endif
    fftZ(BACKWARD);   /* fftZ slices */
#ifdef PRECISE_TIMING
    GetTime(tvp+12);
    ElapsedInc(tvp+11,tvp+12,&Timing_FFTZb);
#endif
    /* copy slice back to Xmatrix */
    for(y=0;y<boxY;y++) for(z=0;z<boxZ;z++) {
      i=index_slice_yz(y,z);
      j=index_garbledX(x,y,z);
      for (Xcomp=0;Xcomp<3;Xcomp++)
        memcpy(Xmatrix[j+Xcomp*local_Nsmall],slices[i+Xcomp*gridYZ],sizeof(doublecomplex));
    }
#ifdef PRECISE_TIMING
    GetTime(tvp+13);
    ElapsedInc(tvp+12,tvp+13,&Timing_Mult4);
#endif
  } /* end of loop over slices */
  /* FFT-X back the result */
  block_transpose(Xmatrix);
#ifdef PRECISE_TIMING
  GetTime(tvp+14);
  elapsed(tvp+13,tvp+14,&Timing_BTb);
#endif
  fftX(BACKWARD);  /* fftX Xmatrix */
#ifdef PRECISE_TIMING
  GetTime(tvp+15);
  elapsed(tvp+14,tvp+15,&Timing_FFTXb);
#endif
  /* fill resultvec */
  for (i=0;i<local_nvoid_Ndip;i++) {
    j=3*i;
    mat=material[i];
    index=index_Xmatrix(position[j],position[j+1],position[j+2]-local_z0);
    for (Xcomp=0;Xcomp<3;Xcomp++) {
      cMult(cc_sqrt[mat][Xcomp],Xmatrix[index+Xcomp*local_Nsmall],temp);
      cAdd(argvec[j+Xcomp],temp,resultvec[j+Xcomp]);        /* result=argvec+cc_sqrt*Xmat */
    }
    if (ipr) *inprod+=cvNorm2(resultvec+j); /* norm is unaffected by conjugation, hence can be computed here */
  }
  if (her) {
    nConj(resultvec);
    nConj(argvec);  /* conjugate back argvec, so it remains unchanged after MatVec */
  }
#ifdef PRECISE_TIMING
  GetTime(tvp+16);
  elapsed(tvp+15,tvp+16,&Timing_Mult5);
#endif
  if (ipr) {
    tstart=clock();
    my_inner_product(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
  mvcounter++;
#ifdef PRECISE_TIMING
  GetTime(tvp+17);
  elapsed(tvp+16,tvp+17,&Timing_ipr);

  SetTimerFreq();
  t_Mult1=TimerToSec(&Timing_Mult1);
  t_Mult2=TimerToSec(&Timing_Mult2);
  t_Mult3=TimerToSec(&Timing_Mult3);
  t_Mult4=TimerToSec(&Timing_Mult4);
  t_Mult5=TimerToSec(&Timing_Mult5);
  t_TYZf=TimerToSec(&Timing_TYZf);
  t_TYZb=TimerToSec(&Timing_TYZb);
  t_BTf=TimerToSec(&Timing_BTf);
  t_BTb=TimerToSec(&Timing_BTb);
  t_FFTXf=TimerToSec(&Timing_FFTXf);
  t_FFTXb=TimerToSec(&Timing_FFTXb);
  t_FFTYf=TimerToSec(&Timing_FFTYf);
  t_FFTYb=TimerToSec(&Timing_FFTYb);
  t_FFTZf=TimerToSec(&Timing_FFTZf);
  t_FFTZb=TimerToSec(&Timing_FFTZb);
  t_ipr=TimerToSec(&Timing_ipr);

  t_Arithm=t_Mult1+t_Mult2+t_Mult3+t_Mult4+t_Mult5+t_TYZf+t_TYZb;
  t_FFT=t_FFTXf+t_FFTYf+t_FFTZf+t_FFTXb+t_FFTYb+t_FFTZb;
  t_Comm=t_BTf+t_BTb+t_ipr;

  sprintz(stat,
         "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"\
         "                MatVec timing              \n"\
         "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"\
         "Arith1 = %4.4f    Arithmetics = %4.4f\n"\
         "FFTXf  = %4.4f    FFT         = %4.4f\n"\
         "BTf    = %4.4f    Comm        = %4.4f\n"\
         "Arith2 = %4.4f\n"\
         "FFTZf  = %4.4f          Total = %4.4f\n"\
         "TYZf   = %4.4f\n"\
         "FFTYf  = %4.4f\n"\
         "Arith3 = %4.4f\n"\
         "FFTYb  = %4.4f\n"\
         "TYZb   = %4.4f\n"\
         "FFTZb  = %4.4f\n"\
         "Arith4 = %4.4f\n"\
         "BTb    = %4.4f\n"\
         "FFTXb  = %4.4f\n"\
         "Arith5 = %4.4f\n"\
         "InProd = %4.4f\n\n",
         t_Mult1,t_Arithm,t_FFTXf,t_FFT,t_BTf,t_Comm,t_Mult2,
         t_FFTZf,DiffSec(tvp,tvp+16),t_TYZf,t_FFTYf,t_Mult3,t_FFTYb,t_TYZb,t_FFTZb,
         t_Mult4,t_BTb,t_FFTXb,t_Mult5,t_ipr);
  printz(stat);
  fprintz(logfile,stat);
  stop(1);
#endif
}
