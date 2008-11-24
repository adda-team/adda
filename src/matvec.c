/* FILE: matvec.c
 * AUTH: Maxim Yurkin
 * DESCR: calculate local matrix vector product of decomposed interaction
 *        matrix with rk or pk, using a FFT based convolution algorithm
 *
 *        Previous version by Michel Grimminck
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdio.h>
#include <string.h>
#include "vars.h"
#include "cmplx.h"
#include "const.h"
#include "comm.h"
#include "fft.h"
#include "prec_time.h"
#include "linalg.h"
#include "function.h"
#include "io.h"

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in fft.c */
extern const doublecomplex *Dmatrix;
extern doublecomplex *Xmatrix,*slices,*slices_tr;
extern const int DsizeY,DsizeZ,DsizeYZ,NDcomp;

/*============================================================*/

INLINE int IndexSliceZY(const int y,const int z)
{
  return (z*gridY+y);
}

/*============================================================*/

INLINE int IndexSliceYZ(const int y,const int z)
{
  return(y*gridZ+z);
}

/*============================================================*/

INLINE int IndexGarbledX(const int x,const int y,const int z)
{
#ifdef PARALLEL
  return(((z%local_Nz)*smallY+y)*gridX+(z/local_Nz)*local_Nx+x%local_Nx);
#else
  return((z*smallY+y)*gridX+x);
#endif
}

/*============================================================*/

INLINE int IndexXmatrix(const int x,const int y,const int z)
{
  return((z*smallY+y)*gridX+x);
}

/*============================================================*/

INLINE int IndexDmatrix_mv(int x,int y,int z,const int transposed)
{
  if (transposed) {   /* used only for G_SO */
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

void MatVec  (doublecomplex *argvec,    /* the argument vector */
              doublecomplex *resultvec, /* the result vector */
              double *inprod,           /* the result inner product */
              const int her)            /* 0 for non-hermitic, 1 for hermetic */
/* This function implements both MatVec_nim and MatVecAndInp_nim.
   the difference is that when we want to calculate the inproduct
   as well, we pass 'inprod' as a non-null pointer. if 'inprod' is
   a NULL, we don't calculate it.
   argvec allways remains unchanged afterwards, however it is not
   strictly const - some manipulations may occur during the execution */
{
  int i, j;

  doublecomplex fmat[6],xvec[3],yvec[3];
  doublecomplex temp;
  int index,x,y,z,Xcomp,ipr;
  unsigned char mat;
  int transposed;
#ifdef PRECISE_TIMING
  SYSTEM_TIME tvp[18];
  SYSTEM_TIME Timing_FFTXf,Timing_FFTYf,Timing_FFTZf,Timing_FFTXb,Timing_FFTYb,Timing_FFTZb,
              Timing_Mult1,Timing_Mult2,Timing_Mult3,Timing_Mult4,Timing_Mult5,
              Timing_BTf,Timing_BTb,Timing_TYZf,Timing_TYZb,Timing_ipr;
  double t_FFTXf,t_FFTYf,t_FFTZf,t_FFTXb,t_FFTYb,t_FFTZb,
         t_Mult1,t_Mult2,t_Mult3,t_Mult4,t_Mult5,t_ipr,
         t_BTf,t_BTb,t_TYZf,t_TYZb,t_Arithm,t_FFT,t_Comm;

#endif

  /* A = I + S.D.S
   * S = sqrt(C)
   * A.x = x + S.D.(S.x)
   * A(H).x = x + (S(T).D(T).S(T).x(*))(*)
   * C,S - diagonal => symmetric
   * (!! will change if tensor (non-diagonal) polarizability is used !!) 
   * D - symmetric (except for G_SO)
   *
   * D.x=F(-1)(F(D).F(X))
   * F(D) is just a vector
   *
   * G_SO: F(D(T)) (k) =  F(D) (-k)
   *       k - vector index
   *
   *   For (her) three additional operations of nConj are used. Should not be a problem,
   *     but can be avoided by a more complex code.
   */

  transposed=(!reduced_FFT) && her;
  if (inprod) ipr=TRUE;
  else ipr=FALSE;

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
  for (i=0;i<3*local_Nsmall;i++) Xmatrix[i][RE]=Xmatrix[i][IM]=0.0;

  /* transform from coordinates to grid and multiply with coupling constant */
  if (her) nConj(argvec); /* conjugated back afterwards */

  for (i=0;i<local_nvoid_Ndip;i++) {
    /* fill grid with argvec*sqrt_cc */
    j=3*i;
    mat=material[i];
    index=IndexXmatrix(position[j],position[j+1],position[j+2]-local_z0);
    for (Xcomp=0;Xcomp<3;Xcomp++)     /* Xmat=cc_sqrt*argvec */
      cMult(cc_sqrt[mat][Xcomp],argvec[j+Xcomp],Xmatrix[index+Xcomp*local_Nsmall]);
  }
#ifdef PRECISE_TIMING
  GetTime(tvp+1);
  elapsed(tvp,tvp+1,&Timing_Mult1);
#endif
  /* FFT X */
  fftX(FFT_FORWARD);    /* fftX Xmatrix */
#ifdef PRECISE_TIMING
  GetTime(tvp+2);
  elapsed(tvp+1,tvp+2,&Timing_FFTXf);
#endif
  BlockTranspose(Xmatrix);
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
    for(i=0;i<3*gridYZ;i++) slices[i][RE]=slices[i][IM]=0.0;

    /* fill slices with values from Xmatrix */
    for(y=0;y<boxY;y++) for(z=0;z<boxZ;z++) {
      i=IndexSliceYZ(y,z);
      j=IndexGarbledX(x,y,z);
      for (Xcomp=0;Xcomp<3;Xcomp++)
        memcpy(slices[i+Xcomp*gridYZ],Xmatrix[j+Xcomp*local_Nsmall],sizeof(doublecomplex));
    }
#ifdef PRECISE_TIMING
    GetTime(tvp+5);
    ElapsedInc(tvp+4,tvp+5,&Timing_Mult2);
#endif
    /* fftZ&Y */
    fftZ(FFT_FORWARD);   /* fftZ slices */
#ifdef PRECISE_TIMING
    GetTime(tvp+6);
    ElapsedInc(tvp+5,tvp+6,&Timing_FFTZf);
#endif
    TransposeYZ(FFT_FORWARD);
#ifdef PRECISE_TIMING
    GetTime(tvp+7);
    ElapsedInc(tvp+6,tvp+7,&Timing_TYZf);
#endif
    fftY(FFT_FORWARD);   /* fftY slices_tr */
#ifdef PRECISE_TIMING
    GetTime(tvp+8);
    ElapsedInc(tvp+7,tvp+8,&Timing_FFTYf);
#endif
    /* do the product D~*X~ */
    for(z=0;z<gridZ;z++) for(y=0;y<gridY;y++) {
      i=IndexSliceZY(y,z);
      for (Xcomp=0;Xcomp<3;Xcomp++)
        memcpy(xvec[Xcomp],slices_tr[i+Xcomp*gridYZ],sizeof(doublecomplex));

      j=IndexDmatrix_mv(x-local_x0,y,z,transposed);
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
      cSymMatrVec(fmat,xvec,yvec);    /* yvec=fmat*xvec */
      for (Xcomp=0;Xcomp<3;Xcomp++)
        memcpy(slices_tr[i+Xcomp*gridYZ],yvec[Xcomp],sizeof(doublecomplex));
    }
#ifdef PRECISE_TIMING
    GetTime(tvp+9);
    ElapsedInc(tvp+8,tvp+9,&Timing_Mult3);
#endif
    /* fft_invY&Z */
    fftY(FFT_BACKWARD);   /* fftY slices_tr */
#ifdef PRECISE_TIMING
    GetTime(tvp+10);
    ElapsedInc(tvp+9,tvp+10,&Timing_FFTYb);
#endif
    TransposeYZ(FFT_BACKWARD);
#ifdef PRECISE_TIMING
    GetTime(tvp+11);
    ElapsedInc(tvp+10,tvp+11,&Timing_TYZb);
#endif
    fftZ(FFT_BACKWARD);   /* fftZ slices */
#ifdef PRECISE_TIMING
    GetTime(tvp+12);
    ElapsedInc(tvp+11,tvp+12,&Timing_FFTZb);
#endif
    /* copy slice back to Xmatrix */
    for(y=0;y<boxY;y++) for(z=0;z<boxZ;z++) {
      i=IndexSliceYZ(y,z);
      j=IndexGarbledX(x,y,z);
      for (Xcomp=0;Xcomp<3;Xcomp++)
        memcpy(Xmatrix[j+Xcomp*local_Nsmall],slices[i+Xcomp*gridYZ],sizeof(doublecomplex));
    }
#ifdef PRECISE_TIMING
    GetTime(tvp+13);
    ElapsedInc(tvp+12,tvp+13,&Timing_Mult4);
#endif
  } /* end of loop over slices */
  /* FFT-X back the result */
  BlockTranspose(Xmatrix);
#ifdef PRECISE_TIMING
  GetTime(tvp+14);
  elapsed(tvp+13,tvp+14,&Timing_BTb);
#endif
  fftX(FFT_BACKWARD);  /* fftX Xmatrix */
#ifdef PRECISE_TIMING
  GetTime(tvp+15);
  elapsed(tvp+14,tvp+15,&Timing_FFTXb);
#endif
  /* fill resultvec */
  for (i=0;i<local_nvoid_Ndip;i++) {
    j=3*i;
    mat=material[i];
    index=IndexXmatrix(position[j],position[j+1],position[j+2]-local_z0);
    for (Xcomp=0;Xcomp<3;Xcomp++) {
      cMult(cc_sqrt[mat][Xcomp],Xmatrix[index+Xcomp*local_Nsmall],temp);
      cAdd(argvec[j+Xcomp],temp,resultvec[j+Xcomp]);      /* result=argvec+cc_sqrt*Xmat */
    }
    /* norm is unaffected by conjugation, hence can be computed here */
    if (ipr) *inprod+=cvNorm2(resultvec+j);
  }
  if (her) {
    nConj(resultvec);
    nConj(argvec);  /* conjugate back argvec, so it remains unchanged after MatVec */
  }
#ifdef PRECISE_TIMING
  GetTime(tvp+16);
  elapsed(tvp+15,tvp+16,&Timing_Mult5);
#endif
  if (ipr) MyInnerProduct(inprod,double_type,1,&Timing_OneIterComm);
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

  PRINTBOTHZ(logfile,
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
  PRINTZ("\nPrecise timing is complete. Finishing execution.\n");
  Stop(0);
#endif
}
