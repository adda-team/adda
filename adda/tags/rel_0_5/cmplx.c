/* FILE: cmplx.c
 * AUTH: Maxim Yurin
 * DESCR: some basic complex, and vector operations;
 *        and several auxillary functions
 * 
 *        this functions may be defined as inline
 */

#include "cmplx.h"
#include "const.h"

/*============================================================*/
/* operations on complex numbers */

FAST double cAbs2(doublecomplex a)
     /* square of absolute value of complex number; |a|^2 */
{
  return (a.r*a.r + a.i*a.i);
}

/*============================================================*/

FAST void cAdd(doublecomplex a,doublecomplex b,doublecomplex *c)
     /* add two complex numbers; c=a+b */
{
  c->r = a.r + b.r;
  c->i = a.i + b.i;
}
                                                              
/*============================================================*/

FAST void cSquare(doublecomplex a,doublecomplex *b)
     /* square of complex number; b=a^2 */
{
  b->r=a.r*a.r - a.i*a.i;
  b->i=2*a.i*a.r;
}

/*============================================================*/

FAST void cMultReal(double a,doublecomplex b,doublecomplex *c)
     /* complex multiplication by real; c=ab */
{
  c->r=a*b.r;
  c->i=a*b.i;
}

/*============================================================*/

FAST void cMult_i(doublecomplex *c)
     /* complex multiplication by i; c=i*c */
{
  double tmp;
  tmp=c->r;
  c->r=-c->i;
  c->i=tmp;
}

/*============================================================*/

FAST void cMult(doublecomplex a,doublecomplex b,doublecomplex *c)
     /* complex multiplication; c=ab */
{
  c->r=a.r*b.r - a.i*b.i;
  c->i=a.i*b.r + a.r*b.i;
}

/*============================================================*/

FAST double cMultConRe(doublecomplex a,doublecomplex b)
     /* complex multiplication; returns real(a*b_conjugated) */
{
  return (a.r*b.r + a.i*b.i);
}

/*============================================================*/

FAST double cMultConIm(doublecomplex a,doublecomplex b)
     /* complex multiplication; returns imag(a*b_conjugated) */
{
  return (a.i*b.r - a.r*b.i);
}

/*============================================================*/

FAST void cLinComb(doublecomplex a,doublecomplex b, double  c1, double c2, doublecomplex *c)
    /* linear combination of two complex numbers; c=c1*a+c2*b */
{
    c->r=c1*a.r+c2*b.r;
    c->i=c1*a.i+c2*b.i;
}

/*============================================================*/

FAST void cInv(doublecomplex a,doublecomplex *b)
     /* complex inversion; b=1/a */
{
  double tmp;
  tmp=1/(a.r*a.r + a.i*a.i);
  b->r= a.r * tmp;
  b->i= - a.i * tmp;
}

/*============================================================*/

FAST void cDiv(doublecomplex a,doublecomplex b,doublecomplex *c)
     /* complex division; c=a/b */
{
  double tmp;
  tmp=1/(b.r*b.r + b.i*b.i);
  c->r=(a.r*b.r + a.i*b.i)*tmp;
  c->i=(a.i*b.r - a.r*b.i)*tmp;
}

/*============================================================*/

FAST void cExp(double arg, doublecomplex *c)
     /* complex exponent c=Exp(i*arg)*/
{
  while (arg>TWO_PI) arg-=TWO_PI;
  c->r = cos(arg);
  c->i = sqrt(1-c->r*c->r);
  if (arg>PI) c->i*=-1;
}

/*============================================================*/

FAST void MultScal(double a,double *b,double *c)
     /* complex multiplication by real; c=ab */
{
  c[0]=a*b[0];
  c[1]=a*b[1];
  c[2]=a*b[2];
}

/*============================================================*/
/* operations on vectors (real and complex) */

FAST void vMult(double *a,double *b,double *c)
     /* multiplication of two vectors (by elements); c[i]=a[i]*b[i] */
{
  c[0]=a[0]*b[0];
  c[1]=a[1]*b[1];
  c[2]=a[2]*b[2];
}

/*============================================================*/

FAST void cScalMultRVec(double *a, doublecomplex b, doublecomplex *c)
    /* complex scalar- real vector[3] multiplication; c=b*a */
{
    c[0].r=b.r*a[0];
    c[0].i=b.i*a[0];
    c[1].r=b.r*a[1];
    c[1].i=b.i*a[1];
    c[2].r=b.r*a[2];
    c[2].i=b.i*a[2];
}

/*============================================================*/

FAST double cvNorm2(doublecomplex *a)  
    /* square of the norm of a complex vector[3] */
{
    return (a[0].r*a[0].r + a[0].i*a[0].i +
            a[1].r*a[1].r + a[1].i*a[1].i +
            a[2].r*a[2].r + a[2].i*a[2].i);
}

/*============================================================*/

FAST double DotProd(double *a, double *b)  
    /* dot product of two real vectors[3] */
{
    return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

/*============================================================*/

FAST void crDotProd(doublecomplex *a, double *b, doublecomplex *c)  
    /* dot product of complex and real vectors[3]; c=a.b */
{
    c->r=a[0].r*b[0]+a[1].r*b[1]+a[2].r*b[2];
    c->i=a[0].i*b[0]+a[1].i*b[1]+a[2].i*b[2];
}

/*============================================================*/

FAST void LinComb(double *a,double *b, double  c1, double c2, double *c)
    /* linear combination of real vectors[3]; c=c1*a+c2*b */
{
    c[0]=c1*a[0]+c2*b[0];
    c[1]=c1*a[1]+c2*b[1];
    c[2]=c1*a[2]+c2*b[2];
}

/*============================================================*/

FAST double TrSym(double *a)
    /* trace of a symmetric matrix stored as a vector of size 6 */
{
  return (a[0]+a[2]+a[5]);
}

/*============================================================*/

FAST double QuadForm(double *matr, double *vec)
    /* value of a quadratic form matr (symmetric matrix stored as 
       a vector of size 6) over a vector vec */
{
  return (vec[0]*vec[0]*matr[0]+vec[1]*vec[1]*matr[2]+vec[2]*vec[2]*matr[5]+
    2*(vec[0]*vec[1]*matr[1]+vec[0]*vec[2]*matr[3]+vec[1]*vec[2]*matr[4]));
}

/*============================================================*/
FAST void MatrVec(double matr[3][3], double *vec, double *res)
     /* multiplication of matrix[3][3] by vec[3] (all real)
        res=matr.vec */
{
  res[0]=matr[0][0]*vec[0]+matr[0][1]*vec[1]+matr[0][2]*vec[2];
  res[1]=matr[1][0]*vec[0]+matr[1][1]*vec[1]+matr[1][2]*vec[2];
  res[2]=matr[2][0]*vec[0]+matr[2][1]*vec[1]+matr[2][2]*vec[2];
}
               
/*============================================================*/

FAST void permutate(double *vec, int *ord)
    /* permutate double vector vec using permutation ord */
{
  double buf[3];

  memcpy(buf,vec,3*sizeof(double));
  vec[0]=buf[ord[0]];
  vec[1]=buf[ord[1]];
  vec[2]=buf[ord[2]];
}

/*============================================================*/

FAST void permutate_i(int *vec, int *ord)
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

FAST double deg2rad(double deg)
   /* transforms angle in degrees to radians */
{
  return (PI*deg/180);
}

/*=====================================================================*/

FAST double rad2deg(double rad)
   /* transforms angle in radians to degrees */ 
{
  return(180.0*rad/PI);
}