/* FILE: cmplx.h
 * AUTH: Maxim Yurkin
 * DESCR: inline complex functions
 *        plus few auxiliary functions
 */
#ifndef __cmplx_h
#define __cmplx_h

#include <string.h>
#include <math.h>
#include "const.h"

/*============================================================*/
/* operations on complex numbers */

INLINE double cAbs2(doublecomplex a)
     /* square of absolute value of complex number; |a|^2 */
{
  return (a[re]*a[re] + a[im]*a[im]);
}

/*============================================================*/

INLINE void cConj(doublecomplex a,doublecomplex b)
     /* complex conjugate; b=a* */
{
  b[re] = a[re];
  b[im] = - a[im];
}

/*============================================================*/

INLINE void cAdd(doublecomplex a,doublecomplex b,doublecomplex c)
     /* add two complex numbers; c=a+b */
{
  c[re] = a[re] + b[re];
  c[im] = a[im] + b[im];
}

/*============================================================*/

INLINE void cSubtr(doublecomplex a,doublecomplex b,doublecomplex c)
     /* subtract two complex numbers; c=a-b */
{
  c[re] = a[re] - b[re];
  c[im] = a[im] - b[im];
}

/*============================================================*/

INLINE void cSquare(doublecomplex a,doublecomplex b)
     /* square of complex number; b=a^2 */
{
  b[re]=a[re]*a[re] - a[im]*a[im];
  b[im]=2*a[im]*a[re];
}

/*============================================================*/

INLINE void cMultReal(double a,doublecomplex b,doublecomplex c)
     /* complex multiplication by real; c=ab */
{
  c[re]=a*b[re];
  c[im]=a*b[im];
}

/*============================================================*/

INLINE void cMult_i(doublecomplex c)
     /* complex multiplication by i; c=i*c */
{
  double tmp;
  tmp=c[re];
  c[re]=-c[im];
  c[im]=tmp;
}

/*============================================================*/

INLINE void cMult(doublecomplex a,doublecomplex b,doublecomplex c)
     /* complex multiplication; c=ab */
     /* !!! c should be different from a and b !!! */
{
  c[re]=a[re]*b[re] - a[im]*b[im];
  c[im]=a[im]*b[re] + a[re]*b[im];
}

/*============================================================*/

INLINE void cMultSelf(doublecomplex a,doublecomplex b)
     /* complex multiplication; a*=b */
{
  double tmp;
  tmp=a[re];
  a[re]=a[re]*b[re] - a[im]*b[im];
  a[im]=a[im]*b[re] + tmp*b[im];
}

/*============================================================*/

INLINE double cMultConRe(doublecomplex a,doublecomplex b)
     /* complex multiplication; returns real(a*b_conjugated) */
{
  return (a[re]*b[re] + a[im]*b[im]);
}

/*============================================================*/

INLINE double cMultConIm(doublecomplex a,doublecomplex b)
     /* complex multiplication; returns imag(a*b_conjugated) */
{
  return (a[im]*b[re] - a[re]*b[im]);
}

/*============================================================*/

INLINE void cLinComb(doublecomplex a,doublecomplex b, double  c1, double c2, doublecomplex c)
    /* linear combination of two complex numbers; c=c1*a+c2*b */
{
    c[re]=c1*a[re]+c2*b[re];
    c[im]=c1*a[im]+c2*b[im];
}

/*============================================================*/

INLINE void cInvSign(doublecomplex a)
     /* change sign of complex number; a*=-1; */
{
  a[re] = - a[re];
  a[im] = - a[im];
}

/*============================================================*/

INLINE void cInvSign2(doublecomplex a,doublecomplex b)
     /* change sign of complex number and store to different address; b=-a; */
{
  b[re] = - a[re];
  b[im] = - a[im];
}

/*============================================================*/

INLINE void cInv(doublecomplex a,doublecomplex b)
     /* complex inversion; b=1/a */
{
  double tmp;
  tmp=1/(a[re]*a[re] + a[im]*a[im]);
  b[re]= a[re] * tmp;
  b[im]= - a[im] * tmp;
}

/*============================================================*/

INLINE void cDiv(doublecomplex a,doublecomplex b,doublecomplex c)
     /* complex division; c=a/b */
     /* !!! c should be different from a and b !!! */
{
  double tmp;
  tmp=1/(b[re]*b[re] + b[im]*b[im]);
  c[re]=(a[re]*b[re] + a[im]*b[im])*tmp;
  c[im]=(a[im]*b[re] - a[re]*b[im])*tmp;
}

/*============================================================*/

INLINE void cDivSelf(doublecomplex a,doublecomplex b)
     /* complex division; a/=b */
{
  double tmp, tmp2;
  tmp=1/(b[re]*b[re] + b[im]*b[im]);
  tmp2=a[re];
  a[re]=(a[re]*b[re] + a[im]*b[im])*tmp;
  a[im]=(a[im]*b[re] - tmp2*b[im])*tmp;
}

/*============================================================*/

INLINE void cSqrt(doublecomplex a, doublecomplex b)
     /* complex square root; b=a
        branch cut discontinuity is (-inf,0) - b[re]>=0 */
{
  double tmp;

  tmp=sqrt(a[re]*a[re]+a[im]*a[im]);
  b[re]=sqrt((tmp+a[re])/2.0);
  b[im]=sqrt((tmp-a[re])/2.0);
  if (a[im]<0) b[re]=-b[re];
}

/*============================================================*/

INLINE void cExp(double arg, doublecomplex c)
     /* complex exponent c=Exp(i*arg)*/
{
  while (arg>=TWO_PI) arg-=TWO_PI;
  while (arg<0) arg+=TWO_PI;

  c[re] = cos(arg);
  c[im] = sqrt(1-c[re]*c[re]);
  if (arg>PI) c[im]=-c[im];
}

/*============================================================*/
/* operations on complex vectors */

INLINE void cvMultScal(double a,doublecomplex *b,doublecomplex *c)
     /* multiplication of vector by real scalar; c=ab */
{
  c[0][re] = a*b[0][re];
  c[0][im] = a*b[0][im];
  c[1][re] = a*b[1][re];
  c[1][im] = a*b[1][im];
  c[2][re] = a*b[2][re];
  c[2][im] = a*b[2][im];
}

/*============================================================*/

INLINE void cScalMultRVec(double *a, doublecomplex b, doublecomplex *c)
    /* complex scalar- real vector[3] multiplication; c=b*a */
{
  c[0][re] = b[re]*a[0];
  c[0][im] = b[im]*a[0];
  c[1][re] = b[re]*a[1];
  c[1][im] = b[im]*a[1];
  c[2][re] = b[re]*a[2];
  c[2][im] = b[im]*a[2];
}

/*============================================================*/

INLINE double cvNorm2(doublecomplex *a)
    /* square of the norm of a complex vector[3] */
{
  return (a[0][re]*a[0][re] + a[0][im]*a[0][im] +
          a[1][re]*a[1][re] + a[1][im]*a[1][im] +
          a[2][re]*a[2][re] + a[2][im]*a[2][im]);
}


/*============================================================*/

INLINE void cDotProd(doublecomplex *a,doublecomplex *b,doublecomplex c)
    /* conjugate dot product of two complex vector[3]; c=a.b */
{
  c[re] = a[0][re]*b[0][re] + a[0][im]*b[0][im] +
          a[1][re]*b[1][re] + a[1][im]*b[1][im] +
          a[2][re]*b[2][re] + a[2][im]*b[2][im];
  c[im] = a[0][im]*b[0][re] - a[0][re]*b[0][im] +
          a[1][im]*b[1][re] - a[1][re]*b[1][im] +
          a[2][im]*b[2][re] - a[2][re]*b[2][im];
}

/*============================================================*/

INLINE void cDotProd_conj(doublecomplex *a,doublecomplex *b,doublecomplex c)
    /* dot product of two complex vector[3]; c=a.b* */
{
  c[re] = a[0][re]*b[0][re] - a[0][im]*b[0][im] +
          a[1][re]*b[1][re] - a[1][im]*b[1][im] +
          a[2][re]*b[2][re] - a[2][im]*b[2][im];
  c[im] = a[0][im]*b[0][re] + a[0][re]*b[0][im] +
          a[1][im]*b[1][re] + a[1][re]*b[1][im] +
          a[2][im]*b[2][re] + a[2][re]*b[2][im];
}

/*============================================================*/

INLINE void cvSubtr(doublecomplex *a,doublecomplex *b,doublecomplex *c)
    /* subtraction of two complex vector[3]; c=a-b */
{
  c[0][re] = a[0][re] - b[0][re];
  c[0][im] = a[0][im] - b[0][im];
  c[1][re] = a[1][re] - b[1][re];
  c[1][im] = a[1][im] - b[1][im];
  c[2][re] = a[2][re] - b[2][re];
  c[2][im] = a[2][im] - b[2][im];
}

/*============================================================*/

INLINE void crDotProd(doublecomplex *a, double *b, doublecomplex c)
    /* dot product of complex and real vectors[3]; c=a.b */
{
  c[re] = a[0][re]*b[0] + a[1][re]*b[1] + a[2][re]*b[2];
  c[im] = a[0][im]*b[0] + a[1][im]*b[1] + a[2][im]*b[2];
}

/*============================================================*/

INLINE void cvIncremScaled_cmplx(doublecomplex *a,doublecomplex b, doublecomplex *c)
    /* increment of complex vectors[3] by complex-scaled other vector; c+=b*a */
{
  c[0][re] += b[re]*a[0][re] - b[im]*a[0][im];
  c[0][im] += b[re]*a[0][im] + b[im]*a[0][re];
  c[1][re] += b[re]*a[1][re] - b[im]*a[1][im];
  c[1][im] += b[re]*a[1][im] + b[im]*a[1][re];
  c[2][re] += b[re]*a[2][re] - b[im]*a[2][im];
  c[2][im] += b[re]*a[2][im] + b[im]*a[2][re];
}
/*============================================================*/

INLINE void cvMultAdd(doublecomplex *a,doublecomplex b, doublecomplex *c)
    /* multiply complex vectors[3] with complex constant and add another vector;
       second coef is unity; c=c1*a+b   !!! c=b*c+a */
{
  double tmp;
  tmp=c[0][re];
  c[0][re] = c[0][re]*b[re] - c[0][im]*b[im] + a[0][re];
  c[0][im] = tmp*b[im] + c[0][im]*b[re] + a[0][im];
  tmp=c[1][re];
  c[1][re] = c[1][re]*b[re] - c[1][im]*b[im] + a[1][re];
  c[1][im] = tmp*b[im] + c[1][im]*b[re] + a[1][im];
  tmp=c[2][re];
  c[2][re] = c[2][re]*b[re] - c[2][im]*b[im] + a[2][re];
  c[2][im] = tmp*b[im] + c[2][im]*b[re] + a[2][im];
}

/*============================================================*/

INLINE void cvLinComb1(doublecomplex *a,doublecomplex *b, double c1, doublecomplex *c)
    /* linear combination of complex vectors[3]; second coef is unity; c=c1*a+b */
{
  c[0][re] = c1*a[0][re] + b[0][re];
  c[0][im] = c1*a[0][im] + b[0][im];
  c[1][re] = c1*a[1][re] + b[1][re];
  c[1][im] = c1*a[1][im] + b[1][im];
  c[2][re] = c1*a[2][re] + b[2][re];
  c[2][im] = c1*a[2][im] + b[2][im];
}

/*============================================================*/

INLINE void cvLinComb1_cmplx(doublecomplex *a,doublecomplex *b, doublecomplex c1, doublecomplex *c)
    /* linear combination of complex vectors[3] with complex coefficients;
       second coef is unity; c=c1*a+b   !!! c!=a */
{
  c[0][re] = a[0][re]*c1[re] - a[0][im]*c1[im] + b[0][re];
  c[0][im] = a[0][re]*c1[im] + a[0][im]*c1[re] + b[0][im];
  c[1][re] = a[1][re]*c1[re] - a[1][im]*c1[im] + b[1][re];
  c[1][im] = a[1][re]*c1[im] + a[1][im]*c1[re] + b[1][im];
  c[2][re] = a[2][re]*c1[re] - a[2][im]*c1[im] + b[2][re];
  c[2][im] = a[2][re]*c1[im] + a[2][im]*c1[re] + b[2][im];
}

/*============================================================*/

INLINE void csymMatrVec(doublecomplex *matr, doublecomplex *vec, doublecomplex *res)
     /* multiplication of complex symmetric matrix[6] by complex vec[3]
        res=matr.vec */
{
  res[0][re] = matr[0][re] * vec[0][re] - matr[0][im] * vec[0][im] +
    matr[1][re] * vec[1][re] - matr[1][im] * vec[1][im] +
    matr[2][re] * vec[2][re] - matr[2][im] * vec[2][im];
  res[0][im] = matr[0][re] * vec[0][im] + matr[0][im] * vec[0][re] +
    matr[1][re] * vec[1][im] + matr[1][im] * vec[1][re] +
    matr[2][re] * vec[2][im] + matr[2][im] * vec[2][re];

  res[1][re] = matr[1][re] * vec[0][re] - matr[1][im] * vec[0][im] +
    matr[3][re] * vec[1][re] - matr[3][im] * vec[1][im] +
    matr[4][re] * vec[2][re] - matr[4][im] * vec[2][im];
  res[1][im] = matr[1][re] * vec[0][im] + matr[1][im] * vec[0][re] +
    matr[3][re] * vec[1][im] + matr[3][im] * vec[1][re] +
    matr[4][re] * vec[2][im] + matr[4][im] * vec[2][re];

  res[2][re] = matr[2][re] * vec[0][re] - matr[2][im] * vec[0][im] +
    matr[4][re] * vec[1][re] - matr[4][im] * vec[1][im] +
    matr[5][re] * vec[2][re] - matr[5][im] * vec[2][im];
  res[2][im] = matr[2][re] * vec[0][im] + matr[2][im] * vec[0][re] +
    matr[4][re] * vec[1][im] + matr[4][im] * vec[1][re] +
    matr[5][re] * vec[2][im] + matr[5][im] * vec[2][re];
}

/*============================================================*/
/* operations on real vectors */

INLINE void MultScal(double a,double *b,double *c)
     /* multiplication of vector by scalar; c=a*b */
{
  c[0]=a*b[0];
  c[1]=a*b[1];
  c[2]=a*b[2];
}

/*============================================================*/

INLINE void vMult(double *a,double *b,double *c)
     /* multiplication of two vectors (by elements); c[i]=a[i]*b[i] */
{
  c[0]=a[0]*b[0];
  c[1]=a[1]*b[1];
  c[2]=a[2]*b[2];
}

/*============================================================*/

INLINE double DotProd(double *a, double *b)
    /* dot product of two real vectors[3] */
{
  return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

/*============================================================*/

INLINE void LinComb(double *a,double *b, double  c1, double c2, double *c)
    /* linear combination of real vectors[3]; c=c1*a+c2*b */
{
  c[0]=c1*a[0]+c2*b[0];
  c[1]=c1*a[1]+c2*b[1];
  c[2]=c1*a[2]+c2*b[2];
}

/*============================================================*/

INLINE double TrSym(double *a)
    /* trace of a symmetric matrix stored as a vector of size 6 */
{
  return (a[0]+a[2]+a[5]);
}

/*============================================================*/

INLINE double QuadForm(double *matr, double *vec)
    /* value of a quadratic form matr (symmetric matrix stored as
       a vector of size 6) over a vector vec */
{
  return (vec[0]*vec[0]*matr[0]+vec[1]*vec[1]*matr[2]+vec[2]*vec[2]*matr[5]+
    2*(vec[0]*vec[1]*matr[1]+vec[0]*vec[2]*matr[3]+vec[1]*vec[2]*matr[4]));
}

/*============================================================*/
INLINE void MatrVec(double matr[3][3], double *vec, double *res)
     /* multiplication of matrix[3][3] by vec[3] (all real)
        res=matr.vec */
{
  res[0]=matr[0][0]*vec[0]+matr[0][1]*vec[1]+matr[0][2]*vec[2];
  res[1]=matr[1][0]*vec[0]+matr[1][1]*vec[1]+matr[1][2]*vec[2];
  res[2]=matr[2][0]*vec[0]+matr[2][1]*vec[1]+matr[2][2]*vec[2];
}

/*============================================================*/

INLINE void permutate(double *vec, int *ord)
    /* permutate double vector vec using permutation ord */
{
  double buf[3];

  memcpy(buf,vec,3*sizeof(double));
  vec[0]=buf[ord[0]];
  vec[1]=buf[ord[1]];
  vec[2]=buf[ord[2]];
}

/*============================================================*/

INLINE void permutate_i(int *vec, int *ord)
    /* permutate int vector vec using permutation ord */
{
  int buf[3];

  memcpy(buf,vec,3*sizeof(int));
  vec[0]=buf[ord[0]];
  vec[1]=buf[ord[1]];
  vec[2]=buf[ord[2]];
}

/*=====================================================================*/
/* Auxillary functions */

INLINE double deg2rad(double deg)
   /* transforms angle in degrees to radians */
{
  return (PI*deg/180);
}

/*=====================================================================*/

INLINE double rad2deg(double rad)
   /* transforms angle in radians to degrees */
{
  return(180.0*rad/PI);
}

#endif /*__cmplx_h*/
